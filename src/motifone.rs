use crate::filesplitpattern::fastareturn;
use crate::nanoporepacbio::Search;
use crate::nanoporepacbio::Sequence;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn motifsearch(path: &str, motif: &str) -> Result<String, Box<dyn Error>> {
    let unpack: Vec<Sequence> = fastareturn(path).unwrap();

    let sequencesearch: Vec<Search> = unpack
        .par_iter()
        .map(|i| {
            let start = i.sequence.find(motif).unwrap();
            Search {
                name: i.header.clone(),
                sequence: i.sequence.clone(),
                start,
                end: start + motif.len(),
            }
        })
        .collect();

    let mut filewrite = BufWriter::new(File::create("motipattern.txt").expect("file not present"));
    for i in sequencesearch.iter() {
        writeln!(filewrite, "{:?}\t{:?}\t{:?}", i.name, i.start, i.end).expect("file not found");
    }

    Ok("The results have been written".to_string())
}
