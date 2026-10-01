use crate::filesplitpattern::fastareturn;
use crate::nanoporepacbio::Motifcatcher;
use crate::nanoporepacbio::Sequence;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn motifcatcherupdown(
    path: &str,
    motif: &str,
    upstream: &str,
    downstream: &str,
) -> Result<String, Box<dyn Error>> {
    let unpack: Vec<Sequence> = fastareturn(path).unwrap();
    let upstreamadd: usize = upstream.parse::<usize>().unwrap();
    let downstreamadd: usize = downstream.parse::<usize>().unwrap();

    let sequencesearch: Vec<Motifcatcher> = unpack
        .par_iter()
        .map(|i| {
            let start = i.sequence.find(motif).unwrap();
            let end = start + motif.len();
            Motifcatcher {
                name: i.header.clone(),
                sequence: i.sequence.clone(),
                start,
                end,
                upstream: i.sequence[start - upstreamadd..start].to_string(),
                downstream: i.sequence[end..downstreamadd].to_string(),
            }
        })
        .collect();

    let mut filewrite = BufWriter::new(File::create("motipattern-up-downstream.txt").expect("file not present"));
    for i in sequencesearch.iter() {
        writeln!(
            filewrite,
            "{:?}\t{:?}\t{:?}\t{:?}\t{:?}",
            i.name, i.start, i.end, i.upstream, i.downstream
        )
        .expect("file not found");
    }

    Ok("The results have been written".to_string())
}
