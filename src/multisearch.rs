use crate::filesplitpattern::fastareturn;
use crate::nanoporepacbio::Multiplepattern;
use crate::nanoporepacbio::Sequence;
use rayon::prelude::*;
use regex::Regex;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn multisearchregex(path: &str, searchstr: &str) -> Result<String, Box<dyn Error>> {
    let sequencevec: Vec<Sequence> = fastareturn(path).unwrap();

    let searchregex = Regex::new(searchstr).unwrap();

    // Regex is Send + Sync, so the compiled pattern can be shared across
    // threads: each sequence's captures_iter scan runs independently.
    let collectionvec: Vec<Multiplepattern> = sequencevec
        .par_iter()
        .map(|i| {
            let vecor: Vec<(usize, usize)> = searchregex
                .captures_iter(i.sequence.as_str())
                .map(|j| (j.get(0).unwrap().start(), j.get(0).unwrap().end()))
                .collect();

            Multiplepattern {
                name: &i.header,
                collectionvec: vecor,
            }
        })
        .collect();

    let mut filewrite = BufWriter::new(File::create("multi-searchpattern.txt").expect("file not found"));

    for i in collectionvec.iter() {
        for j in i.collectionvec.iter() {
            writeln!(filewrite, "{}\t{}\t{}", i.name, j.0, j.1).expect("line not found");
        }
    }

    Ok("The search results for the multisearch have been written".to_string())
}
