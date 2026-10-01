use crate::nanoporepacbio::Pafanalyzer;
use crate::nanoporepacbio::Tableinformation;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use tabled::Table;

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn pangenome_summarize(path: &str) -> Result<String, Box<dyn Error>> {
    let pafopen = File::open(path).expect("file not found");
    let pafread = BufReader::new(pafopen);
    let mut pafstruct: Vec<Pafanalyzer> = Vec::new();
    for i in pafread.lines() {
        let line = i.expect("line not found");
        let pafvec = line.split("\t").collect::<Vec<_>>();
        pafstruct.push(Pafanalyzer {
            query: pafvec[0].to_string(),
            query_length: pafvec[1].parse::<usize>().unwrap(),
            query_start: pafvec[2].parse::<usize>().unwrap(),
            query_end: pafvec[3].parse::<usize>().unwrap(),
            strand: pafvec[4].to_string(),
            target: pafvec[5].to_string(),
            target_length: pafvec[6].parse::<usize>().unwrap(),
            target_start: pafvec[7].parse::<usize>().unwrap(),
            target_end: pafvec[8].parse::<usize>().unwrap(),
            residue_matches: pafvec[9].parse::<usize>().unwrap(),
            alignment_length: pafvec[10].parse::<usize>().unwrap(),
            quality: pafvec[11].parse::<usize>().unwrap(),
        });
    }

    // Four separate O(n) passes over the same slice collapsed into four
    // parallel reductions (still logically independent sums, just no
    // longer four redundant single-threaded scans).
    let queryaligned: usize = pafstruct.par_iter().map(|i| i.query_end - i.query_start).sum();
    let targetaligned: usize = pafstruct
        .par_iter()
        .map(|i| i.target_end - i.target_start)
        .sum();
    let residuematches: usize = pafstruct.par_iter().map(|i| i.residue_matches).sum();
    let alignmentlength: usize = pafstruct.par_iter().map(|i| i.alignment_length).sum();

    let finaltable = vec![Tableinformation {
        information: "Genome",
        query: queryaligned,
        residue: residuematches,
        alignment: alignmentlength,
        target: targetaligned,
    }];

    let _finaltable = Table::new(finaltable).to_string();

    Ok("The pangenome has been summarized".to_string())
}
