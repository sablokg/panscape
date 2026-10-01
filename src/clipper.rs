use crate::filesplitpattern::fastareturn;
use crate::nanoporepacbio::Filesnatch;
use crate::nanoporepacbio::Sequence;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn clipperpattern(path: &str, regionfile: &str) -> Result<String, Box<dyn Error>> {
    let unpack: Vec<Sequence> = fastareturn(path).unwrap();
    let filecontent = File::open(regionfile).expect("file not present");
    let fileread = BufReader::new(filecontent);

    let mut fileanalyze: Vec<Filesnatch> = Vec::new();
    for i in fileread.lines() {
        let line = i.expect("line not present");
        let vecline = line.split(",").collect::<Vec<_>>();
        fileanalyze.push(Filesnatch {
            name: vecline[0].to_string(),
            start: vecline[1].parse::<usize>().unwrap(),
            end: vecline[2].parse::<usize>().unwrap(),
        })
    }

    // Outer loop parallelized: each read is matched against the region list
    // independently, so the O(n*m) comparisons split across threads.
    let sequencesearch: Vec<Sequence> = unpack
        .par_iter()
        .flat_map(|i| {
            fileanalyze
                .iter()
                .filter(|j| i.header == j.name)
                .map(|j| Sequence {
                    header: i.header.clone(),
                    sequence: i.sequence[j.start..j.end].to_string(),
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let mut filewrite = BufWriter::new(File::create("cliipedpattern.fasta").expect("file not present"));
    for i in sequencesearch.iter() {
        writeln!(filewrite, "{:?}\t{:?}", i.header, i.sequence).expect("file not found");
    }

    Ok("The results have been written".to_string())
}
