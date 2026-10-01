use crate::filesplitpattern::fastareturn;
use crate::nanoporepacbio::Readlength;
use crate::nanoporepacbio::Sequence;
use rayon::prelude::*;
use std::error::Error;

/*
Gaurav Sablok,
gsablok@proton.me
*/

pub fn stats(path: &str) -> Result<String, Box<dyn Error>> {
    let sequencevector: Vec<Sequence> = fastareturn(path).unwrap();

    // `i.sequence.to_string().len()` was cloning every read's full sequence
    // (up to 5 times per read, once per branch already checked) just to
    // read its length - `.sequence` is already a String, so `.len()` needs
    // no clone at all. The five-way bucket count is also now a single
    // parallel reduce instead of a sequential branch-counting loop.
    let (fivelength, hundredlength, hundredfiftylength, twohundredlength, morethanlength) =
        sequencevector
            .par_iter()
            .map(|i| {
                let len = i.sequence.len();
                if len <= 50000usize {
                    (1usize, 0usize, 0usize, 0usize, 0usize)
                } else if len <= 100000usize {
                    (0, 1, 0, 0, 0)
                } else if len <= 150000usize {
                    (0, 0, 1, 0, 0)
                } else if len <= 200000usize {
                    (0, 0, 0, 1, 0)
                } else {
                    (0, 0, 0, 0, 1)
                }
            })
            .reduce(
                || (0, 0, 0, 0, 0),
                |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2, a.3 + b.3, a.4 + b.4),
            );

    let readvec = vec![Readlength {
        five: fivelength,
        hundred: hundredlength,
        hundredfifty: hundredfiftylength,
        twohundred: twohundredlength,
        morethan: morethanlength,
    }];

    for i in readvec.iter() {
        println!("Read length greater than 50KB:{}\tRead length greater than 100KB:{}\tRead length greater than 150KB:{}\tRead length greater than 200KB{}\tReads length more than that:{}", i.five, i.hundred, i.hundredfifty, i.twohundred, i.morethan);
    }

    Ok("The result have been printed as follows".to_string())
}
