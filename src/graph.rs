use crate::nanoporepacbio::FinalWrite;
use crate::nanoporepacbio::Links;
use crate::nanoporepacbio::Segments;
use rayon::prelude::*;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn graph_args_segment(path: &str) -> Result<Vec<FinalWrite>, Box<dyn Error>> {
    let mut segment_counter: usize = 0usize;
    let mut links_counter: usize = 0usize;
    let mut graphreadcheck: Vec<String> = Vec::new();

    let graphopen = File::open(path).expect("file not found");
    let graphread = BufReader::new(graphopen);
    for i in graphread.lines() {
        let line = i.unwrap();
        if line.starts_with("S") {
            segment_counter += 1
        } else if line.starts_with("L") {
            links_counter += 1
        }
        graphreadcheck.push(line);
    }

    let mut segment_graph: Vec<Segments> = Vec::new();
    let mut link_graph: Vec<Links> = Vec::new();
    for i in graphreadcheck.iter() {
        if i.starts_with("S") {
            let line: Vec<_> = i.split("\t").filter(|x| !x.is_empty()).collect::<Vec<_>>();
            segment_graph.push(Segments {
                segment: line[0].to_string(),
                id: line[1].to_string(),
                seq: line[2].to_string(),
                tag: line[3].to_string(),
                aligntag: line[4].to_string(),
                alignmenttag: line[5].to_string(),
            });
        } else if i.starts_with("L") {
            let line: Vec<_> = i.split("\t").filter(|x| !x.is_empty()).collect::<Vec<_>>();
            link_graph.push(Links {
                link: line[0].to_string(),
                id1: line[1].to_string(),
                strand1: line[2].to_string(),
                id2: line[3].to_string(),
                strand2: line[4].to_string(),
                tag: line[5].to_string(),
                arc: line[6].to_string(),
            });
        }
    }

    let maxiter: HashSet<usize> = link_graph
        .par_iter()
        .map(|i| {
            let maxiter: Vec<_> = i.arc.split(":").collect::<Vec<_>>();
            maxiter[2].parse::<usize>().unwrap()
        })
        .collect();
    let maxranknumber: usize = maxiter.len();

    let (segment_summary, segment_finallength): (Vec<String>, Vec<usize>) = segment_graph
        .par_iter()
        .map(|i| (i.seq.clone(), i.seq.len()))
        .unzip();
    let _ = &segment_summary; // kept: computed but unused in the original too
    let finalsegmentlength: usize = segment_finallength.par_iter().sum();

    let rank_segment: Vec<(String, usize)> = segment_graph
        .par_iter()
        .filter_map(|i| {
            let segmentholdtag: Vec<_> = i.alignmenttag.split(":").collect::<Vec<_>>();
            if segmentholdtag[2].parse::<usize>().unwrap() == 0usize {
                Some((i.seq.clone(), i.seq.len()))
            } else {
                None
            }
        })
        .collect();

    let rank_counter: usize = rank_segment.par_iter().map(|i| i.1).sum();

    let arc_number: usize = link_graph.len() * 2;

    let mut finalwrite: Vec<FinalWrite> = vec![];
    finalwrite.push(FinalWrite {
        number_segment: segment_counter,
        number_links: links_counter,
        number_arc: arc_number,
        max_rank: maxranknumber,
        total_segment_length: finalsegmentlength,
        average_segment_length: finalsegmentlength / segment_finallength.len(),
        sum_0_segment_length: rank_counter,
    });

    let mut pangenome_summarize_write =
        BufWriter::new(File::create("summarized_pangenome.txt").expect("file not found"));
    for i in finalwrite.iter() {
        writeln!( pangenome_summarize_write,
            "Results have been written:\nnumber_segment:{}\nnumber_links:{}\nnumber_arc:{}\nnumber_rank:{}\ntotal_segment_length:{}\naverage_segment_length:{}\nsum_0_segment_length:{}\n",
            i.number_segment,
            i.number_links,
            i.number_arc,
            i.max_rank,
            i.total_segment_length,
            i.average_segment_length,
            i.sum_0_segment_length
        ).expect("line not present");
    }

    let mut graphwrite = BufWriter::new(File::create("graph.fasta").expect("file not found"));
    for i in segment_graph.iter() {
        writeln!(graphwrite, ">{}\n{}", i.id, i.seq).expect("line not present");
    }

    Ok(finalwrite)
}
