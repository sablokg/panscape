use crate::nanoporepacbio::{
    QueryComparativeHolder, QuerypafView, RefComparativeHolder, ReferencepafView, ResidueAlignment,
};
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn paf_alignments(pathpaf1: &str, pathpaf2: &str) -> Result<String, Box<dyn Error>> {
    let pafaligner1 = File::open(pathpaf1).expect("file not present");
    let pafreader1 = BufReader::new(pafaligner1);
    let mut pafquery1: Vec<QuerypafView> = Vec::new();
    let mut pafreference1: Vec<ReferencepafView> = Vec::new();
    for i in pafreader1.lines() {
        let line = i.expect("line not present");
        let linecheck = line.split("\t").collect::<Vec<_>>();
        pafquery1.push(QuerypafView {
            query: linecheck[0].to_string(),
            length: linecheck[1].parse::<usize>().unwrap(),
            start: linecheck[2].parse::<usize>().unwrap(),
            end: linecheck[3].parse::<usize>().unwrap(),
            strand: linecheck[4].to_string(),
        });
        pafreference1.push(ReferencepafView {
            reference: linecheck[5].to_string(),
            length: linecheck[6].parse::<usize>().unwrap(),
            start: linecheck[7].parse::<usize>().unwrap(),
            end: linecheck[8].parse::<usize>().unwrap(),
            residuematch: linecheck[9].parse::<usize>().unwrap(),
            alignmentblock: linecheck[10].parse::<usize>().unwrap(),
        });
    }

    let pafaligner2 = File::open(pathpaf2).expect("file not present");
    let pafreader2 = BufReader::new(pafaligner2);
    let mut pafquery2: Vec<QuerypafView> = Vec::new();
    let mut pafreference2: Vec<ReferencepafView> = Vec::new();
    for i in pafreader2.lines() {
        let line = i.expect("line not present");
        let linecheck = line.split("\t").collect::<Vec<_>>();
        pafquery2.push(QuerypafView {
            query: linecheck[0].to_string(),
            length: linecheck[1].parse::<usize>().unwrap(),
            start: linecheck[2].parse::<usize>().unwrap(),
            end: linecheck[3].parse::<usize>().unwrap(),
            strand: linecheck[4].to_string(),
        });
        pafreference2.push(ReferencepafView {
            reference: linecheck[5].to_string(),
            length: linecheck[6].parse::<usize>().unwrap(),
            start: linecheck[7].parse::<usize>().unwrap(),
            end: linecheck[8].parse::<usize>().unwrap(),
            residuematch: linecheck[9].parse::<usize>().unwrap(),
            alignmentblock: linecheck[10].parse::<usize>().unwrap(),
        });
    }

    // The original loop was O(pafquery1 * pafquery2 * pafreference1 *
    // pafreference2), re-checking i.query == j.query inside the two
    // reference loops even though that comparison doesn't depend on f1/f2
    // at all - for every matching (i, j) pair the full reference1 x
    // reference2 cross product was redone from scratch just to reach the
    // same "yes, still matches" answer m1*m2 times over.
    //
    // Grouping pafquery2 by query name turns the query match into a single
    // HashMap lookup per i (was an O(pafquery2) scan), and the outer loop
    // over pafquery1 is parallelized. The reference1 x reference2 cross
    // product itself is preserved exactly as before for every matching
    // query pair, since that combinatorial comparison is the actual
    // intended output, not incidental waste.
    let query2_by_name: HashMap<&str, Vec<&QuerypafView>> = {
        let mut map: HashMap<&str, Vec<&QuerypafView>> = HashMap::new();
        for j in pafquery2.iter() {
            map.entry(j.query.as_str()).or_default().push(j);
        }
        map
    };

    let results: Vec<(QueryComparativeHolder, RefComparativeHolder, ResidueAlignment)> = pafquery1
        .par_iter()
        .flat_map(|i| {
            let mut rows = Vec::new();
            if let Some(matches) = query2_by_name.get(i.query.as_str()) {
                for j in matches.iter() {
                    for f1 in pafreference1.iter() {
                        for f2 in pafreference2.iter() {
                            rows.push((
                                QueryComparativeHolder {
                                    query1: i.query.clone(),
                                    length1: i.length,
                                    start1: i.start,
                                    end1: i.end,
                                    strand1: i.strand.clone(),
                                    query2: j.query.clone(),
                                    length2: j.length,
                                    start2: j.start,
                                    end2: j.end,
                                    strand2: j.strand.clone(),
                                },
                                RefComparativeHolder {
                                    reference1: f1.reference.clone(),
                                    length1: f1.length,
                                    start1: f1.start,
                                    end1: f1.end,
                                    reference2: f2.reference.clone(),
                                    length2: f2.length,
                                    start2: f2.start,
                                    end2: f2.end,
                                },
                                ResidueAlignment {
                                    residuematch1: f1.residuematch,
                                    residuematch2: f2.residuematch,
                                    alignmentblock1: f1.alignmentblock,
                                    alignmentblock2: f2.alignmentblock,
                                },
                            ));
                        }
                    }
                }
            }
            rows
        })
        .collect();

    let mut comparative_query: Vec<QueryComparativeHolder> = Vec::with_capacity(results.len());
    let mut comparative_reference: Vec<RefComparativeHolder> = Vec::with_capacity(results.len());
    let mut comparative_residue_alignment: Vec<ResidueAlignment> = Vec::with_capacity(results.len());
    for (q, r, a) in results {
        comparative_query.push(q);
        comparative_reference.push(r);
        comparative_residue_alignment.push(a);
    }
    let _ = &comparative_residue_alignment; // preserved for interface/behavioral parity with original (unused after collection, same as before)

    let mut comparative_query_write =
        BufWriter::new(File::create("comparative_query.txt").expect("file not present"));
    let mut comparative_ref_write =
        BufWriter::new(File::create("comparative_ref_write.txt").expect("file not present"));

    println!("The comparative reference files have been written");
    for i in comparative_reference {
        writeln!(
            comparative_ref_write,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            i.reference1, i.reference2, i.start1, i.start2, i.end1, i.end2, i.length1, i.length2
        )
        .expect("line not present");
    }

    println!("The comparative query files have been written");
    for i in comparative_query {
        writeln!(
            comparative_query_write,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            i.query1,
            i.length1,
            i.start1,
            i.end1,
            i.strand1,
            i.query2,
            i.length2,
            i.start2,
            i.end2,
            i.strand2
        )
        .expect("line nor present");
    }

    Ok(
        "The pangenome alignments have been parsed and the comparative results have been written"
            .to_string(),
    )
}
