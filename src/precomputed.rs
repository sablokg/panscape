use crate::nanoporepacbio::{QuerypafView, ReferencepafView};
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn precomputealignments(pathpaf: &str) -> Result<String, Box<dyn Error>> {
    let pafaligner = File::open(pathpaf).expect("file not present");
    let pafreader = BufReader::new(pafaligner);
    let mut pafquery: Vec<QuerypafView> = Vec::new();
    let mut pafreference: Vec<ReferencepafView> = Vec::new();
    for i in pafreader.lines() {
        let line = i.expect("line not present");
        let linecheck = line.split("\t").collect::<Vec<_>>();
        pafquery.push(QuerypafView {
            query: linecheck[0].to_string(),
            length: linecheck[1].parse::<usize>().unwrap(),
            start: linecheck[2].parse::<usize>().unwrap(),
            end: linecheck[3].parse::<usize>().unwrap(),
            strand: linecheck[4].to_string(),
        });
        pafreference.push(ReferencepafView {
            reference: linecheck[5].to_string(),
            length: linecheck[6].parse::<usize>().unwrap(),
            start: linecheck[7].parse::<usize>().unwrap(),
            end: linecheck[8].parse::<usize>().unwrap(),
            residuematch: linecheck[9].parse::<usize>().unwrap(),
            alignmentblock: linecheck[10].parse::<usize>().unwrap(),
        })
    }

    let total_length_query: Vec<usize> = pafquery.par_iter().map(|i| i.end - i.start).collect();
    let total_length_reference: Vec<usize> =
        pafreference.par_iter().map(|i| i.end - i.start).collect();

    let pafreference_sum = sum_methodcall(total_length_reference.clone()).unwrap();
    let average_reference = mode_methodcall(total_length_reference.clone()).unwrap();
    let average_query = mode_methodcall(total_length_query.clone()).unwrap();
    let pafquery_sum = sum_methodcall(total_length_query.clone()).unwrap();
    println!(
        "The total reference aligned sum for the pangenome is {:?}",
        pafreference_sum
    );
    println!(
        "The total reference portion of the aligned query is {:?}",
        pafquery_sum
    );
    println!(
        "The mode for the reference pangenome alignment is: {:?}",
        average_reference
    );
    println!(
        "The model for the query pangenome alignment is {:?}",
        average_query
    );

    Ok(
        "The paf analyzer has finished and all the pangenome alignment statistics are given below:"
            .to_string(),
    )
}

fn sum_methodcall(inputarray: Vec<usize>) -> Result<usize, ()> {
    Ok(inputarray.par_iter().sum())
}

fn mode_methodcall(inputarray: Vec<usize>) -> Result<usize, ()> {
    // NOTE: despite the name, this computes an arithmetic mean (sum / count),
    // not a statistical mode - kept exactly as-is, just parallelized.
    let sum: usize = inputarray.par_iter().sum();
    Ok(sum / inputarray.len())
}
