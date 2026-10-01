use crate::filesplitpattern::fastareturn;
use crate::nanoporepacbio::{CdsExtract, CdsExtractSeq, Sequence};
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn computecds<'a>(pathalignment: &str, pathreadsfasta: &str) -> Result<String, Box<dyn Error>> {
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

    // The original code looped over every unique id AND every line for each
    // id (O(ids * lines)), even though each line already carries its own id
    // in linecheck[0] - the outer id loop only ever matched a single id per
    // line, so it was redoing the same O(lines) work once per unique id for
    // no benefit. Zipping the id captured at parse time with its line does
    // the same match in a single O(lines) pass, and that pass is
    // parallelized since each line's result is independent.
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

    // HashMap lookup replaces the former inner scan over sequenceundone.
    let seq_by_header: HashMap<&str, &Sequence> = sequenceundone
        .iter()
        .map(|seq| (seq.header.as_str(), seq))
        .collect();

    let extractseq: Vec<CdsExtractSeq> = linevecstore
        .par_iter()
        .zip(alignmentids.par_iter())
        .filter_map(|(line, id)| {
            let linecheck: Vec<&str> = line.split('\t').collect();
            if linecheck[2] != "CDS" {
                return None;
            }
            let seq = seq_by_header.get(id.as_str())?;
            let start = linecheck[3].parse::<usize>().unwrap();
            let end = linecheck[4].parse::<usize>().unwrap();
            Some(CdsExtractSeq {
                header: id,
                cdsordinate: vec![(start, end)],
                cdsextract: vec![seq.sequence[start..end].to_string()],
            })
        })
        .collect();

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
