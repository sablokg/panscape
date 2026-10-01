use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me
*/
pub fn harmonicestimate(harmonicfile: &str) -> Result<String, Box<dyn Error>> {
    let fileopen1 = File::open(harmonicfile).expect("file not present");
    let fileread1 = BufReader::new(fileopen1);
    let lines: Vec<String> = fileread1
        .lines()
        .map(|i| i.expect("line not present"))
        .collect();

    // `harmonic_capture_inverse` was computed (integer 1/x, always 0 for
    // x > 1) but never used - dropped. Each line's harmonic value is
    // independent of the others, so the per-line computation is
    // parallelized; only the file read itself stays sequential.
    let harmonic_capture_final: Vec<usize> = lines
        .par_iter()
        .map(|line| {
            let linevector: Vec<usize> = line
                .split(" ")
                .map(|x| x.parse::<usize>().unwrap())
                .collect();
            let len = linevector.len();
            len / sumharmonic(linevector).unwrap()
        })
        .collect();

    let mut harmonic_file = BufWriter::new(File::create("harmonic-analysis.txt").expect("file not present"));
    for i in 0..harmonic_capture_final.len() {
        write!(harmonic_file, "{}\n", harmonic_capture_final[i]).expect("line not present");
    }

    Ok("capture harmonic mean for each same query coordinate have ben writen".to_string())
}

pub fn sumharmonic(inputvec: Vec<usize>) -> Result<usize, Box<dyn Error>> {
    Ok(inputvec.par_iter().sum())
}
