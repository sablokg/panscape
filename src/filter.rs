use crate::filesplitpattern::fastareturn;
use crate::nanoporepacbio::Sequence;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn readlength(path: &str, length: &str) -> Result<String, Box<dyn Error>> {
    let readiter: Vec<Sequence> = fastareturn(path).unwrap();
    let lengththreshold: usize = length.parse::<usize>().unwrap();

    let filterseq: Vec<Sequence> = readiter
        .par_iter()
        .filter(|i| i.sequence.len() >= lengththreshold)
        .map(|i| Sequence {
            header: i.header.clone(),
            sequence: i.sequence.clone(),
        })
        .collect();

    let mut filewrite = BufWriter::new(File::create("filtered-ones.fasta").expect("file not found"));
    for i in filterseq.iter() {
        writeln!(filewrite, ">{}\n{}", i.header, i.sequence).expect("line not found");
    }

    Ok("The filtered reads have been written".to_string())
}
