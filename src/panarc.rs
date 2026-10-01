use crate::nanoporepacbio::{
    AlignmentGFF, CdsAnnotate, MrnaAnnotate, Negative, NegativeAnnotate, Positive,
    PositiveAnnotate, Sequenceadd,
};
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn metagenome_annotate(path: &str, fasta: &str) -> Result<String, Box<dyn Error>> {
    let mut vectorhold = Vec::new();
    let mut vectorstring: Vec<AlignmentGFF> = Vec::new();
    let f = File::open(path).expect("file not present");
    let read = BufReader::new(f);
    for gffreadline in read.lines() {
        let gffline = gffreadline.expect("line not present");
        if gffline.starts_with("#") {
            continue;
        } else {
            vectorhold.push(gffline)
        }
    }

    for i in vectorhold.iter() {
        let addline = i.split("\t").collect::<Vec<&str>>();
        let idhold = addline[0];
        let genomefeaturehold = addline[2];
        let starthold = addline[3]
            .to_string()
            .parse::<usize>()
            .expect("number not present");
        let endhold = addline[4]
            .to_string()
            .parse::<usize>()
            .expect("number not present");
        let strandhold = addline[6].to_string();
        vectorstring.push(AlignmentGFF {
            id: idhold.to_string(),
            genomefeature: genomefeaturehold.to_string(),
            start: starthold,
            end: endhold,
            strand: strandhold.to_string(),
        })
    }

    let mut positive = Vec::new();
    let mut negative = Vec::new();
    let new_positive = vectorstring.clone();
    let new_negative = vectorstring.clone();
    for i in new_positive.into_iter() {
        if i.strand == "+" {
            positive.push(Positive {
                id: i.id,
                genomefeature: i.genomefeature,
                start: i.start,
                end: i.end,
                strand: i.strand,
            })
        }
    }
    for i in new_negative.into_iter() {
        if i.strand == "-" {
            negative.push(Negative {
                id: i.id,
                genomefeature: i.genomefeature,
                start: i.start,
                end: i.end,
                strand: i.strand,
            })
        }
    }

    let mut header = vec![];
    let mut sequence = vec![];
    let f = File::open(fasta).expect("file not present");
    let read = BufReader::new(f);
    for i in read.lines() {
        let line = i.expect("line not present");
        if line.starts_with(">") {
            header.push(line)
        } else {
            sequence.push(line)
        }
    }

    let mut final_seq: Vec<Sequenceadd> = Vec::new();
    for i in 0..header.len() {
        final_seq.push(Sequenceadd {
            id: header[i].to_string(),
            sequence: sequence[i].to_string(),
        })
    }

    // NOTE: as in the original code, these four passes match every
    // final_seq record against every gff feature of the right type/strand -
    // there's no `i.id == j.id` check tying a sequence to its own feature,
    // so this is a genuine full cross product, not something a HashMap
    // join can shrink algorithmically (there's no matching key used).
    // Preserved exactly; what's optimized here is doing a single filter
    // pass per feature type (was being re-scanned from scratch for each of
    // the four output lists) and parallelizing the outer loop over
    // final_seq, since each final_seq record's work is independent.
    let mrna_features: Vec<&AlignmentGFF> = vectorstring
        .iter()
        .filter(|j| j.genomefeature == "mRNA")
        .collect();
    let cds_features: Vec<&AlignmentGFF> = vectorstring
        .iter()
        .filter(|j| j.genomefeature == "CDS")
        .collect();
    let cds_pos_features: Vec<&AlignmentGFF> = cds_features
        .iter()
        .copied()
        .filter(|j| j.strand == "+")
        .collect();
    let cds_neg_features: Vec<&AlignmentGFF> = cds_features
        .iter()
        .copied()
        .filter(|j| j.strand == "-")
        .collect();

    let m_rna_struct: Vec<MrnaAnnotate> = final_seq
        .par_iter()
        .flat_map(|i| {
            mrna_features
                .iter()
                .map(|j| MrnaAnnotate {
                    id: j.id.to_string(),
                    seq: i.sequence[j.start - 1..j.end].to_string(),
                    strand: j.strand.to_string(),
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let cdsstruct: Vec<CdsAnnotate> = final_seq
        .par_iter()
        .flat_map(|i| {
            cds_features
                .iter()
                .map(|j| CdsAnnotate {
                    id: j.id.to_string(),
                    seq: i.sequence[j.start - 1..j.end].to_string(),
                    strand: j.strand.to_string(),
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let positive: Vec<PositiveAnnotate> = final_seq
        .par_iter()
        .flat_map(|i| {
            cds_pos_features
                .iter()
                .map(|j| PositiveAnnotate {
                    id: j.id.to_string(),
                    seq: i.sequence[j.start - 1..j.end].to_string(),
                    strand: j.strand.to_string(),
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let negative: Vec<NegativeAnnotate> = final_seq
        .par_iter()
        .flat_map(|i| {
            cds_neg_features
                .iter()
                .map(|j| NegativeAnnotate {
                    id: j.id.to_string(),
                    seq: i.sequence[j.start - 1..j.end].to_string(),
                    strand: j.strand.to_string(),
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let _mrna_length: Vec<usize> = m_rna_struct.par_iter().map(|i| i.seq.len()).collect();
    let _cds_length: Vec<usize> = cdsstruct.par_iter().map(|i| i.seq.len()).collect();
    let _cds_length_positive: Vec<usize> = positive.par_iter().map(|i| i.seq.len()).collect();
    let _cds_length_negative: Vec<usize> = negative.par_iter().map(|i| i.seq.len()).collect();

    let mut mrna_file = BufWriter::new(File::create("mRNA.fasta").expect("file not present"));
    for i in m_rna_struct.iter() {
        writeln!(mrna_file, ">{}\n{}", i.id, i.seq).expect("line not found");
    }

    let mut cds_file = BufWriter::new(File::create("cds.fasta").expect("file not present"));
    for i in cdsstruct.iter() {
        writeln!(cds_file, ">{}\n{}", i.id, i.seq).expect("line not found");
    }

    let mut cds_positive = BufWriter::new(File::create("cds-positive.fasta").expect("file not present"));
    for i in positive.iter() {
        writeln!(cds_positive, ">{}\n{}", i.id, i.seq).expect("line not found");
    }

    let mut cds_negative = BufWriter::new(File::create("cds-negative.fasta").expect("file not present"));
    for i in negative.iter() {
        writeln!(cds_negative, ">{}\n{}", i.id, i.seq).expect("line not found");
    }
    Ok("The panarc has been completed".to_string())
}
