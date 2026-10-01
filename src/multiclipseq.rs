use crate::filesplitpattern::fastareturn;
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

pub fn multiclipseqa(path: &str, clipseq: &str) -> Result<String, Box<dyn Error>> {
    let sequenceunpack: Vec<Sequence> = fastareturn(path).unwrap();

    let fileopen = File::open(clipseq).expect("file not present");
    let fileread = BufReader::new(fileopen);
    let mut clipregion: Vec<String> = Vec::new();

    for i in fileread.lines() {
        let line = i.expect("file not present");
        clipregion.push(line.to_string());
    }

    // The inner chain (each clip step searches within the previous step's
    // result) is inherently sequential per-read, but different reads are
    // fully independent, so the outer loop over sequenceunpack is what gets
    // parallelized.
    let clipseqa: Vec<Sequence> = sequenceunpack
        .par_iter()
        .flat_map(|i| {
            let mut perread: Vec<Sequence> = Vec::new();
            for j in 0..clipregion.len() {
                let seqlinestart = i.sequence.find(clipregion[j].as_str()).unwrap();
                let seqlineend = seqlinestart + clipregion[j].len();
                let clippedregion = i.sequence[seqlinestart..seqlineend].to_string();
                let nextclippedregion_start =
                    clippedregion.find(clipregion[j + 1].as_str()).unwrap();
                let nextclippedreg_end = nextclippedregion_start + clipregion[j + 1].len();
                let finalnextclipped =
                    clippedregion[nextclippedregion_start..nextclippedreg_end].to_string();
                let finalclippedstart = finalnextclipped.find(clipregion[j + 2].as_str()).unwrap();
                let finalclippedend = finalclippedstart + clipregion[j + 2].len();
                let finalclip = finalnextclipped[finalclippedstart..finalclippedend].to_string();
                let writeclipstart = finalclip.find(clipregion[j + 3].as_str()).unwrap();
                let writeclipend = writeclipstart + clipregion[j + 3].len();
                let writeseq = finalclip[writeclipstart..writeclipend].to_string();
                if !clippedregion.is_empty() {
                    perread.push(Sequence {
                        header: i.header.clone(),
                        sequence: writeseq,
                    });
                }
            }
            perread
        })
        .collect();

    let mut filewrite = BufWriter::new(File::create("clippedseq-file.fasta").expect("file not found"));
    for i in clipseqa.iter() {
        writeln!(filewrite, ">{}\n{}", i.header, i.sequence).expect("file not present");
    }

    Ok("The clipseq has been finished".to_string())
}
