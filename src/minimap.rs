use crate::filesplitpattern::fastareturn;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;

/*
Gaurav Sablok,
gsablok@proton.me
*/
pub fn minimapalignment(
    pathreads: &str,
    pathminimap: &str,
    pathproteins: &str,
    thread: &str,
) -> Result<String, Box<dyn Error>> {
    let sequencevector = fastareturn(pathreads).unwrap();
    let proteinfileopen = File::open(pathproteins).expect("file not present");
    let proteinfileread = BufReader::new(proteinfileopen);
    let mut proteinheader: Vec<String> = Vec::new();
    let mut proteinsequence: Vec<String> = Vec::new();

    for i in proteinfileread.lines() {
        let line = i.expect("line not found");
        if line.starts_with(">") {
            proteinheader.push(line.replace(">", ""));
        } else if !line.starts_with(">") {
            proteinsequence.push(line)
        }
    }

    let _command = Command::new(pathminimap)
        .arg("-t")
        .arg(thread)
        .arg(pathreads)
        .arg(pathproteins)
        .arg(">")
        .arg("proteinaligned.fasta")
        .output()
        .expect("command failed");

    let pafopen = File::open("proteinaligned.fasta").expect("file not present");

    let pafread = BufReader::new(pafopen);

    let mut m_rna_vec: Vec<(String, usize, usize)> = Vec::new();

    for i in pafread.lines() {
        let line = i.expect("line not present");
        if line.starts_with("#") {
            continue;
        } else if !line.starts_with("#") {
            let pafline = line.split("\t").collect::<Vec<_>>();
            if pafline[2].to_string() == "mRNA" {
                let paftuple = (
                    pafline[0].to_string(),
                    pafline[3].parse::<usize>().unwrap(),
                    pafline[4].parse::<usize>().unwrap(),
                );
                m_rna_vec.push(paftuple);
            }
        }
    }

    // HashMap lookup replaces the O(m_rna_vec * sequencevector) scan; the
    // per-match slice work is parallelized, and BTreeMap insertion (which
    // must stay sequential/ordered) happens afterward from the collected
    // results.
    let seq_by_header: HashMap<&str, &crate::nanoporepacbio::Sequence> = sequencevector
        .iter()
        .map(|j| (j.header.as_str(), j))
        .collect();

    let matched: Vec<(String, (usize, usize, String))> = m_rna_vec
        .par_iter()
        .filter_map(|i| {
            seq_by_header.get(i.0.as_str()).map(|j| {
                (
                    i.0.clone(),
                    (i.1, i.2, j.sequence[i.1..i.2].to_string()),
                )
            })
        })
        .collect();

    let mut storemap: BTreeMap<String, (usize, usize, String)> = BTreeMap::new();
    for (k, v) in matched {
        storemap.insert(k, v);
    }

    let mut filewrite =
        BufWriter::new(File::create("mrna-predicted-reads-alignment.fasta").expect("file not found"));

    for i in storemap.iter() {
        writeln!(filewrite, ">{}-{}-{}\n{}", i.0, i.1 .0, i.1 .2, i.1 .2).expect("file not found");
    }

    Ok("The results have been sent to the files: {:?}".to_string())
}
