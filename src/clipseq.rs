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

pub fn clipseqa(path: &str, clipseq: &str) -> Result<String, Box<dyn Error>> {
    let sequenceunpack: Vec<Sequence> = fastareturn(path).unwrap();

    let clipseqa: Vec<Sequence> = sequenceunpack
        .par_iter()
        .filter_map(|i| {
            let seqlinestart = i.sequence.find(clipseq).unwrap();
            let seqlineend = seqlinestart + clipseq.len();
            let clippedregion = i.sequence[seqlinestart..seqlineend].to_string();
            if clippedregion.is_empty() {
                None
            } else {
                Some(Sequence {
                    header: i.header.clone(),
                    sequence: clippedregion,
                })
            }
        })
        .collect();

    let mut filewrite = BufWriter::new(File::create("clippedseq-file.fasta").expect("file not found"));
    for i in clipseqa.iter() {
        writeln!(filewrite, ">{}\n{}", i.header, i.sequence).expect("file not present");
    }

    Ok("The clipseq has been finished".to_string())
}
