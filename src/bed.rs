use crate::nanoporepacbio::Analyze;
use crate::nanoporepacbio::Segment;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn graph_bed_segment(path: &str) -> Result<Vec<Analyze>, Box<dyn Error>> {
    let graphopen = File::open(path).expect("file not found");
    let graphread = BufReader::new(graphopen);
    let mut segmentopen: Vec<Segment> = Vec::new();

    for i in graphread.lines() {
        let line = i.expect("line not found");
        if line.starts_with("S") {
            let linevec = line
                .split("\t")
                .filter(|x| !x.is_empty())
                .collect::<Vec<_>>();
            segmentopen.push(Segment {
                segment: linevec[0].to_string(),
                id: linevec[1].to_string(),
                seq: linevec[2].to_string(),
                connection: linevec[3].to_string(),
                tag: linevec[5].to_string(),
            });
        }
    }

    let segmentanalyzer: Vec<Analyze> = segmentopen
        .par_iter()
        .map(|i| {
            let idnumber: usize =
                i.id.split(|c: char| !c.is_numeric())
                    .collect::<Vec<_>>()
                    .into_iter()
                    .filter(|x| !x.is_empty())
                    .collect::<Vec<_>>()
                    .join("")
                    .parse::<usize>()
                    .unwrap();
            let idalignmenttag: String = i.connection.split(":").collect::<Vec<_>>()[2].to_string();
            let rankdefine: usize = i.tag.split(":").collect::<Vec<_>>()[2]
                .parse::<usize>()
                .unwrap();
            Analyze {
                tag: idalignmenttag,
                start: idnumber,
                end: idnumber + i.seq.len(),
                oritag: i.id.clone(),
                rank: rankdefine,
            }
        })
        .collect();

    let mut samfile = BufWriter::new(File::create("gfa-bed.txt").expect("file not found"));

    for i in segmentanalyzer.iter() {
        writeln!(
            samfile,
            "{}\t{}\t{}\t{}\t{}",
            i.tag, i.start, i.end, i.oritag, i.rank
        )
        .expect("line not present");
    }

    Ok(segmentanalyzer)
}
