use crate::nanoporepacbio::Sequence;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn fastaconvertall(path: &str) -> Result<String, Box<dyn Error>> {
    let fileopen = File::open(path).expect("path not found");
    let fileread = BufReader::new(fileopen);
    let lines: Vec<String> = fileread
        .lines()
        .map(|i| i.expect("line not present"))
        .collect();

    // Same record-based parser as filesplitpattern::fastareturn - reads one
    // sequence line per '@' header, and skips a following '+'-line's
    // quality line if present, instead of classifying every line
    // independently (which previously mispaired headers with quality
    // strings whenever quality lines were present).
    let mut sequencehead: Vec<String> = Vec::new();
    let mut sequenceseq: Vec<String> = Vec::new();
    let mut idx = 0usize;
    while idx < lines.len() {
        if let Some(header) = lines[idx].strip_prefix('@') {
            sequencehead.push(header.to_string());
            idx += 1;
            if idx < lines.len() {
                sequenceseq.push(lines[idx].clone());
                idx += 1;
            }
            if idx < lines.len() && lines[idx].starts_with('+') {
                idx += 1;
                if idx < lines.len() {
                    idx += 1;
                }
            }
        } else {
            idx += 1;
        }
    }

    let combinedseq: Vec<Sequence> = sequencehead
        .into_par_iter()
        .zip(sequenceseq.into_par_iter())
        .map(|(header, sequence)| Sequence { header, sequence })
        .collect();

    let mut filewrite = BufWriter::new(File::create("fastaconvert.fasta").expect("file not present"));
    for i in combinedseq.iter() {
        writeln!(filewrite, ">{}\n{}", i.header, i.sequence).expect("file not found");
    }

    Ok("The files have been converted".to_string())
}
