use crate::filesplitpattern::fastareturn;
use crate::nanoporepacbio::Sequence;
use rayon::prelude::*;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn selected(path: &str, ids: &str) -> Result<String, Box<dyn Error>> {
    let unpackreads: Vec<Sequence> = fastareturn(path).unwrap();
    let mut idscheck: Vec<String> = Vec::new();

    let idsopen = File::open(ids).expect("file not present");
    let idsread = BufReader::new(idsopen);
    for i in idsread.lines() {
        let line = i.expect("line not present");
        idscheck.push(line.trim().to_string());
    }

    // HashSet lookup turns the previous O(reads * ids) scan into O(reads),
    // and the per-read work is independent so it's split across threads too.
    let idsset: HashSet<&str> = idscheck.iter().map(|s| s.as_str()).collect();
    let selectedones: Vec<Sequence> = unpackreads
        .par_iter()
        .filter(|i| idsset.contains(i.header.as_str()))
        .map(|i| Sequence {
            header: i.header.clone(),
            sequence: i.sequence.clone(),
        })
        .collect();

    let mut filewrite = BufWriter::new(File::create("selectedones.fasta").expect("file not present"));
    for i in selectedones.iter() {
        writeln!(filewrite, ">{}\n{}", i.header, i.sequence).expect("file not present");
    }

    Ok("The selected reads have been written".to_string())
}
