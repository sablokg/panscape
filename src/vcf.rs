use crate::nanoporepacbio::Common;
use crate::nanoporepacbio::Commonsnatcher;
use crate::nanoporepacbio::Endpointcompare;
use crate::nanoporepacbio::Fasta;
use crate::nanoporepacbio::Fastasnatcher;
use crate::nanoporepacbio::Startpointcompare;
use crate::nanoporepacbio::VCFRange;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/*
Gaurav Sablok,
codeprog@icloud.com
*/
pub fn vcf_compare(
    path1: &str,
    path2: &str,
    pathfastaunwrap: &str,
) -> Result<String, Box<dyn Error>> {
    let vcf1_open = File::open(path1).expect("file not found");
    let vcf2_open = File::open(path2).expect("file not found");
    let vcf1_read = BufReader::new(vcf1_open);
    let vcf2_read = BufReader::new(vcf2_open);
    let mut vcf1: Vec<VCFRange> = Vec::new();
    let mut vcf2: Vec<VCFRange> = Vec::new();

    for i in vcf1_read.lines() {
        let line = i.expect("line not present");
        if line.starts_with("#") {
            continue;
        } else if !line.starts_with("#") {
            let linevec = line
                .split("\t")
                .filter(|x| !x.is_empty())
                .collect::<Vec<_>>();
            vcf1.push(VCFRange {
                name: linevec[0].to_string(),
                start: linevec[1].parse::<usize>().unwrap(),
                end: linevec[2].parse::<usize>().unwrap(),
                delorig: linevec[3].to_string(),
                deltype: linevec[4].to_string().replace(">", "").replace("<", ""),
                threshold: Box::new(linevec[5].parse::<f64>().unwrap()),
            });
        }
    }

    for i in vcf2_read.lines() {
        let line = i.expect("line not present");
        if line.starts_with("#") {
            continue;
        } else if !line.starts_with("#") {
            let linevec = line
                .split("\t")
                .filter(|x| !x.is_empty())
                .collect::<Vec<_>>();
            vcf2.push(VCFRange {
                name: linevec[0].to_string(),
                start: linevec[1].parse::<usize>().unwrap(),
                end: linevec[2].parse::<usize>().unwrap(),
                delorig: linevec[3].to_string(),
                deltype: linevec[4].to_string().replace(">", "").replace("<", ""),
                threshold: Box::new(linevec[5].parse::<f64>().unwrap()),
            });
        }
    }

    // Each of the three vcf1 x vcf2 comparisons below is gated on an
    // equality (`end==end`, `start==start`, or both) before the inequality
    // that actually decides which branch to push. Grouping vcf2 by that
    // equality key turns the O(vcf1 * vcf2) scan into O(vcf1 + vcf2) hash
    // lookups, and the remaining per-match work is independent across vcf1
    // entries so it's parallelized too.
    let vcf2_by_end: HashMap<usize, Vec<&VCFRange>> = {
        let mut map: HashMap<usize, Vec<&VCFRange>> = HashMap::new();
        for value in vcf2.iter() {
            map.entry(value.end).or_default().push(value);
        }
        map
    };
    let vcf2_by_start: HashMap<usize, Vec<&VCFRange>> = {
        let mut map: HashMap<usize, Vec<&VCFRange>> = HashMap::new();
        for value in vcf2.iter() {
            map.entry(value.start).or_default().push(value);
        }
        map
    };
    let vcf2_by_start_end: HashMap<(usize, usize), Vec<&VCFRange>> = {
        let mut map: HashMap<(usize, usize), Vec<&VCFRange>> = HashMap::new();
        for value in vcf2.iter() {
            map.entry((value.start, value.end)).or_default().push(value);
        }
        map
    };

    let end_point_compare: Vec<Endpointcompare> = vcf1
        .par_iter()
        .flat_map(|j| {
            let mut rows = Vec::new();
            if let Some(candidates) = vcf2_by_end.get(&j.end) {
                for value in candidates.iter() {
                    if j.start != value.start {
                        let (start2, difference) = if j.start > value.start {
                            (value.start, j.start - value.start)
                        } else {
                            (value.start, value.start - j.start)
                        };
                        rows.push(Endpointcompare {
                            name: j.name.clone(),
                            start1: j.start,
                            end1: j.end,
                            start2,
                            end2: value.end,
                            delorig1: j.delorig.clone(),
                            deltype1: j.deltype.clone(),
                            threshold1: Box::new(*j.threshold),
                            delorig2: value.delorig.clone(),
                            deltype2: value.deltype.clone(),
                            threshold2: Box::new(*value.threshold),
                            difference: Box::new(difference),
                        });
                    }
                }
            }
            rows
        })
        .collect();

    let start_point_compare: Vec<Startpointcompare> = vcf1
        .par_iter()
        .flat_map(|j| {
            let mut rows = Vec::new();
            if let Some(candidates) = vcf2_by_start.get(&j.start) {
                for value in candidates.iter() {
                    if j.end != value.end {
                        let difference = if j.end > value.end {
                            j.end - value.end
                        } else {
                            value.end - j.end
                        };
                        rows.push(Startpointcompare {
                            name: j.name.clone(),
                            start1: j.start,
                            start2: value.start,
                            end1: j.end,
                            end2: value.end,
                            delorig1: j.delorig.clone(),
                            deltype1: j.deltype.clone(),
                            threshold1: Box::new(*j.threshold),
                            delorig2: value.delorig.clone(),
                            deltype2: value.deltype.clone(),
                            threshold2: Box::new(*value.threshold),
                            difference: Box::new(difference),
                        });
                    }
                }
            }
            rows
        })
        .collect();

    let common_pangenome_variants: Vec<Common> = vcf1
        .par_iter()
        .flat_map(|j| {
            vcf2_by_start_end
                .get(&(j.start, j.end))
                .into_iter()
                .flatten()
                .map(|value| Common {
                    name: j.name.clone(),
                    start1: j.start,
                    start2: value.start,
                    end1: j.end,
                    end2: value.end,
                    delorig1: j.delorig.clone(),
                    deltype1: j.deltype.clone(),
                    delorig2: value.delorig.clone(),
                    deltype2: value.deltype.clone(),
                    threshold1: *j.threshold,
                    threshold2: *value.threshold,
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let fasta_unload: Vec<Fasta> = fasta_estimate(pathfastaunwrap).unwrap();

    let mut endpointcompare = BufWriter::new(File::create("endpoint_variants.txt").expect("file not present"));
    let mut startpointcompare = BufWriter::new(File::create("startpoint_variants.txt").expect("file not present"));
    let mut common_point = BufWriter::new(File::create("common_variant.txt").expect("file not present"));

    // Every remaining join below is keyed on `name == header`, so one
    // HashMap built once replaces the repeated linear scans over
    // fasta_unload, and each independent per-row lookup+slice is
    // parallelized.
    let fasta_by_header: HashMap<&str, &Fasta> = fasta_unload
        .iter()
        .map(|j| (j.header.as_str(), j))
        .collect();

    let commonpointvariant: Vec<Commonsnatcher> = common_pangenome_variants
        .par_iter()
        .filter_map(|i| {
            fasta_by_header.get(i.name.as_str()).map(|j| Commonsnatcher {
                name: i.name.clone(),
                start1: i.start1,
                end1: i.end1,
                start2: i.start2,
                end2: i.end2,
                deltype1: i.deltype1.clone(),
                delorig1: i.delorig1.clone(),
                deltype2: i.deltype2.clone(),
                delorig2: i.delorig2.clone(),
                threshold1: Box::new(i.threshold1),
                threshold2: Box::new(i.threshold2),
                sequenceregion1: j.sequence[i.start1..i.end1].to_string(),
                sequenceregion2: j.sequence[i.start2..i.end2].to_string(),
            })
        })
        .collect();

    for i in commonpointvariant.iter() {
        writeln!(
            common_point,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            i.name,
            i.start1,
            i.end1,
            i.start2,
            i.end2,
            i.delorig1,
            i.deltype1,
            i.delorig2,
            i.deltype2,
            i.threshold1,
            i.threshold2,
            i.sequenceregion1,
            i.sequenceregion2
        )
        .expect("file not present");
    }

    // implemented reverse slicing when the coorindate are not range compatible.

    // Same HashMap join as commonpointvariant above, replacing the O(n*m)
    // scan against fasta_unload; each row's reverse-slice work is
    // independent so it's parallelized.
    let endpointsnatcher: Vec<Fastasnatcher> = end_point_compare
        .par_iter()
        .filter_map(|i| {
            let j = fasta_by_header.get(i.name.as_str())?;
            if i.end1 == i.end2 && i.start1 > i.start2 {
                Some(Fastasnatcher {
                    name: i.name.clone(),
                    start1: i.start1,
                    end1: i.end1,
                    start2: i.start2,
                    end2: i.end2,
                    delorig1: i.delorig1.clone(),
                    deltype1: i.deltype1.clone(),
                    delorig2: i.delorig2.clone(),
                    deltype2: i.deltype2.clone(),
                    threshold1: i.threshold1.clone(),
                    threshold2: i.threshold2.clone(),
                    sequenceadd: j.sequence[i.start2..i.start1].chars().rev().collect(),
                    sequenceregion1: j.sequence[i.start1..i.end1].to_string(),
                    sequenceregion2: j.sequence[i.start2..i.end2].to_string(),
                })
            } else if i.end1 == i.end2 && i.start1 < i.start2 {
                Some(Fastasnatcher {
                    name: i.name.clone(),
                    start1: i.start1,
                    end1: i.end1,
                    start2: i.start2,
                    end2: i.end2,
                    delorig1: i.delorig1.clone(),
                    deltype1: i.deltype1.clone(),
                    delorig2: i.delorig2.clone(),
                    deltype2: i.deltype2.clone(),
                    threshold1: i.threshold1.clone(),
                    threshold2: i.threshold2.clone(),
                    sequenceadd: j.sequence[i.start1..i.start2].to_string(),
                    sequenceregion1: j.sequence[i.start1..i.end1].to_string(),
                    sequenceregion2: j.sequence[i.start2..i.end2].to_string(),
                })
            } else {
                None
            }
        })
        .collect();

    // implemented the reverse slicing when the range parameters are not in the range.

    let startpointsnatcher: Vec<Fastasnatcher> = start_point_compare
        .par_iter()
        .filter_map(|i| {
            let j = fasta_by_header.get(i.name.as_str())?;
            if i.start1 == i.start2 && i.end1 > i.end2 {
                Some(Fastasnatcher {
                    name: i.name.clone(),
                    start1: i.start1,
                    end1: i.end1,
                    start2: i.start2,
                    end2: i.end2,
                    deltype1: i.deltype1.clone(),
                    delorig1: i.delorig1.clone(),
                    deltype2: i.deltype2.clone(),
                    delorig2: i.delorig2.clone(),
                    threshold1: i.threshold1.clone(),
                    threshold2: i.threshold2.clone(),
                    sequenceadd: j.sequence[i.end2..i.end1].chars().rev().collect(),
                    sequenceregion1: j.sequence[i.start1..i.end1].to_string(),
                    sequenceregion2: j.sequence[i.start2..i.end2].to_string(),
                })
            } else if i.start1 == i.start2 && i.end1 < i.end2 {
                Some(Fastasnatcher {
                    name: i.name.clone(),
                    start1: i.start1,
                    end1: i.end1,
                    start2: i.start2,
                    end2: i.end2,
                    delorig1: i.delorig1.clone(),
                    deltype1: i.deltype1.clone(),
                    delorig2: i.delorig2.clone(),
                    deltype2: i.deltype2.clone(),
                    threshold1: i.threshold1.clone(),
                    threshold2: i.threshold2.clone(),
                    sequenceadd: j.sequence[i.end1..i.end2].to_string(),
                    sequenceregion1: j.sequence[i.start1..i.end1].to_string(),
                    sequenceregion2: j.sequence[i.start2..i.end2].to_string(),
                })
            } else {
                None
            }
        })
        .collect();

    for i in endpointsnatcher.iter() {
        writeln!(
            endpointcompare,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            i.name,
            i.start1,
            i.end1,
            i.start2,
            i.end2,
            i.delorig1,
            i.deltype1,
            i.delorig2,
            i.deltype2,
            i.sequenceadd,
            i.sequenceregion1,
            i.sequenceregion2
        )
        .expect("line not present");
    }

    for i in startpointsnatcher.iter() {
        writeln!(
            startpointcompare,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            i.name,
            i.start1,
            i.end1,
            i.start2,
            i.end2,
            i.delorig1,
            i.deltype1,
            i.delorig2,
            i.deltype2,
            i.sequenceadd,
            i.sequenceregion1,
            i.sequenceregion2
        )
        .expect("line not present");
    }

    Ok("VF file have been compared for the pangenomes".to_string())
}

pub fn fasta_estimate(pathfasta: &str) -> Result<Vec<Fasta>, Box<dyn Error>> {
    let fastaopen = File::open(pathfasta).expect("file not present");
    let fastaread = BufReader::new(fastaopen);
    let mut fastaholder: Vec<Fasta> = Vec::new();
    let mut fastaheader: Vec<String> = Vec::new();
    let mut fastasequence: Vec<String> = Vec::new();
    for i in fastaread.lines() {
        let line = i.expect("line not present");
        if line.starts_with(">") {
            fastaheader.push(line.replace(">", ""));
        } else {
            fastasequence.push(line);
        }
    }

    for i in 0..fastaheader.len() {
        fastaholder.push(Fasta {
            header: fastaheader[i].clone(),
            sequence: fastasequence[i].clone(),
        })
    }

    Ok(fastaholder)
}
