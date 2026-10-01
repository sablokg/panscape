use crate::filesplitpattern::fastareturn;
use crate::nanoporepacbio::{CdsExtract, CdsExtractSeq, Sequence};
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn computeintergenic(
    pathalignment: &str,
    pathreadsfasta: &str,
) -> Result<String, Box<dyn Error>> {
    let sequenceundone: Vec<Sequence> = fastareturn(pathreadsfasta).unwrap();

    let genomescale = File::open(pathalignment).expect("file not present");
    let genomeread = BufReader::new(genomescale);

    let mut linevecstore: Vec<String> = Vec::new();
    let mut alignmentids: Vec<String> = Vec::new();

    for i in genomeread.lines() {
        let line = i.expect("line not found");
        let linevec: Vec<String> = line.split("\t").map(|x| x.to_string()).collect::<Vec<_>>();
        alignmentids.push(linevec[0].clone());
        linevecstore.push(line);
    }

    // Same simplification as cdsextract::computecds: each line already
    // carries its own id, so matching it against every unique id (former
    // O(ids * lines)) was redundant. One parallel pass over the zipped
    // (line, id) pairs replaces it.
    let extract: Vec<CdsExtract> = linevecstore
        .par_iter()
        .zip(alignmentids.par_iter())
        .filter_map(|(line, id)| {
            let linecheck: Vec<&str> = line.split('\t').collect();
            if linecheck[2] == "CDS" {
                let cootuple: (usize, usize) = (
                    linecheck[3].parse::<usize>().unwrap(),
                    linecheck[4].parse::<usize>().unwrap(),
                );
                Some(CdsExtract {
                    header: id,
                    cdsordinate: vec![cootuple],
                })
            } else {
                None
            }
        })
        .collect();

    let _intergenicslice: Vec<Vec<usize>> = extract
        .par_iter()
        .map(|i| {
            (0..i.cdsordinate.len().saturating_sub(1))
                .map(|j| i.cdsordinate[j + 1].0 - i.cdsordinate[j].1)
                .collect()
        })
        .collect();

    // NOTE: each `extract` entry is built above with exactly one coordinate
    // tuple in `cdsordinate` (`vec![cootuple]`, same as the original code's
    // `cdsall.push(cootuple)` with a fresh `cdsall` per CDS line). That
    // means `i.cdsordinate.len() - 1` is always 0, so this loop's
    // `0..i.cdsordinate.len() - 1` range was always empty and `extractseq`
    // was always produced empty - the loop iterated over every
    // extract x sequenceundone pair without ever doing any work. Speed-wise
    // that's a lot of wasted iteration for a guaranteed no-op, so it's
    // short-circuited here to keep the identical (empty) output without
    // paying for the scan. If per-gene intergenic regions spanning multiple
    // CDS coordinates were actually intended, `extract` would need to group
    // coordinates by id first (e.g. into a HashMap<&str, Vec<(usize,usize)>>)
    // rather than push a single-tuple entry per CDS line - happy to wire
    // that up if that's the intended behavior.
    let extractseq: Vec<CdsExtractSeq> = Vec::new();
    let _ = &sequenceundone; // kept for interface stability / future grouping

    let mut extractwrite = BufWriter::new(File::create("extract-bed-cds-coordinate.txt").expect("file not found"));
    let mut extractseqwrite = BufWriter::new(File::create("cds-ccordinates.txt").expect("file not found"));

    for i in extract.iter() {
        writeln!(extractwrite, "{}\t{:?}", i.header, i.cdsordinate).expect("line not found");
    }

    for i in extractseq.iter() {
        writeln!(
            extractseqwrite,
            "{}\t{:?}\t{:?}",
            i.header, i.cdsordinate, i.cdsextract
        )
        .expect("line not found");
    }

    Ok("The results have been written".to_string())
}
