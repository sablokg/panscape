use crate::nanoporepacbio::Sequence;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn fastareturn(path: &str) -> Result<Vec<Sequence>, Box<dyn Error>> {
    let fileopen = File::open(path).expect("path not found");
    let fileread = BufReader::new(fileopen);
    let lines: Vec<String> = fileread
        .lines()
        .map(|i| i.expect("line not present"))
        .collect();

    // Record-based parser: each record starts at a '@' header line and is
    // followed by exactly one sequence line. If the line after that starts
    // with '+' (a FASTQ separator), the following quality line is skipped
    // too, so quality strings are never mistaken for sequence data. This
    // handles both plain 2-line FASTA-style records and 4-line FASTQ
    // records in the same pass, without assuming a fixed line-count for
    // the whole file (the previous version classified every line
    // independently, which silently mispaired headers with quality
    // strings once quality lines were present - see prior discussion).
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
                idx += 1; // skip the '+' separator line
                if idx < lines.len() {
                    idx += 1; // skip the quality line that follows it
                }
            }
        } else {
            // Not a header where one was expected (blank line, stray '-'
            // line, malformed record, etc.) - skip it defensively rather
            // than misreading it as sequence data.
            idx += 1;
        }
    }

    let mut combinedseq: Vec<Sequence> = Vec::with_capacity(sequencehead.len());
    for (header, sequence) in sequencehead.into_iter().zip(sequenceseq.into_iter()) {
        combinedseq.push(Sequence { header, sequence });
    }

    Ok(combinedseq)
}
