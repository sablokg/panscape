mod args;
mod bed;
mod bedpangenome;
mod cdsextract;
mod clipper;
mod clipseq;
mod estimate;
mod fastaconvert;
mod filesplitpattern;
mod filter;
mod graph;
mod intergenic;
mod matcher;
mod minimap;
mod motifcatcher;
mod motifone;
mod multiclipseq;
mod multimerge;
mod multisearch;
mod nanoporepacbio;
mod pafnet;
mod pafpangenome;
mod panarc;
mod pangenome;
mod pangenomereadsdatabase;
mod precomputed;
mod selectedreads;
mod singlemerge;
mod snatcher;
mod stat;
mod vcf;
use crate::args::CommandParse;
use crate::args::Commands;
use crate::bed::graph_bed_segment;
use crate::bedpangenome::pangenome_longest_alignment;
use crate::cdsextract::computecds;
use crate::clipper::clipperpattern;
use crate::clipseq::clipseqa;
use crate::estimate::harmonicestimate;
use crate::fastaconvert::fastaconvertall;
use crate::filter::readlength;
use crate::graph::graph_args_segment;
use crate::intergenic::computeintergenic;
use crate::matcher::paf_alignments;
use crate::minimap::minimapalignment;
use crate::motifcatcher::motifcatcherupdown;
use crate::motifone::motifsearch;
use crate::multiclipseq::multiclipseqa;
use crate::multimerge::pangenomemultialignment;
use crate::multisearch::multisearchregex;
use crate::pafnet::pafnet_convert;
use crate::pafpangenome::pangenome_summarize;
use crate::panarc::metagenome_annotate;
use crate::pangenome::pangenome_hifiasm;
use crate::precomputed::precomputealignments;
use crate::selectedreads::selected;
use crate::singlemerge::pangenome_merge;
use crate::snatcher::snatcherextract;
use crate::stat::stats;
use crate::vcf::vcf_compare;
use clap::Parser;
use figlet_rs::FIGfont;
use pangenomereadsdatabase::readsdatabase;

/*
Gaurav Sablok,
gsablok@proton.me
*/

fn main() {
    let standard_font = FIGfont::standard().unwrap();
    let figure = standard_font.convert("panscape");
    assert!(figure.is_some());
    println!("{}", figure.unwrap());
    let argsparse = CommandParse::parse();
    match &argsparse.command {
        Commands::FastaConvert { fastqfile, thread } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = fastaconvertall(fastqfile).unwrap();
                println!("The file has been converted: {:?}", command);
            });
        }
        Commands::ClipperAlign {
            fastqfile,
            regionfile,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = clipperpattern(fastqfile, regionfile).unwrap();
                println!("The regions have been clipped: {:?}", command);
            });
        }
        Commands::Scanner {
            fastqfile,
            motifscan,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = motifsearch(fastqfile, motifscan).unwrap();
                println!("The scanned motif files have been written:{:?}", command);
            });
        }
        Commands::Motifcatcher {
            fastqfile,
            motifscan,
            upstream,
            downstream,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command =
                    motifcatcherupdown(fastqfile, motifscan, upstream, downstream).unwrap();
                println!(
                    "The results of the motifcatcher have been written:{:?}",
                    command
                );
            });
        }
        Commands::Filterreads {
            fastqfile,
            length,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = readlength(fastqfile, length).unwrap();
                println!("The reads have been filtered: {:?}", command);
            });
        }
        Commands::Selectedreads {
            fastqfile,
            ids,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = selected(fastqfile, ids).unwrap();
                println!(
                    "The file with the selected ids have been written:{:?}",
                    command
                );
            });
        }
        Commands::ClipSeq {
            fastqfile,
            clippedseq,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = clipseqa(fastqfile, clippedseq).unwrap();
                println!("The regions have been clipped:{:?}", command);
            });
        }
        Commands::MultiClipSeq {
            fastqfile,
            clipseqfile,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = multiclipseqa(fastqfile, clipseqfile).unwrap();
                println!("The file with the clipseq have been clipped:{:?}", command);
            });
        }
        Commands::Minimap {
            fastqfile,
            minimap,
            proteins,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = minimapalignment(fastqfile, minimap, proteins, thread).unwrap();
                println!("The results of the same are as follows:{:?}", command);
            });
        }
        Commands::Pangenome {
            fastqfile,
            thread,
            proteinfasta,
        } => {
            let command = pangenome_hifiasm(fastqfile, thread, proteinfasta).unwrap();
            println!(
                "The pangenome has been assembled using the pacbiohifi reads:{:?}",
                command
            );
        }
        Commands::Stat { fastqfile, thread } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = stats(fastqfile).unwrap();
                println!("The stats for your file is as follows:{:?}", command);
            });
        }
        Commands::PangenomeSummarize { pangenome, thread } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = pangenome_summarize(pangenome).unwrap();
                println!("The pangenome has been summarized:{:?}", command);
            });
        }
        Commands::ReadMultisearch {
            readsgenome,
            multisearch,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = multisearchregex(readsgenome, multisearch).unwrap();
                println!("The read search have been written:{:?}", command);
            });
        }
        Commands::Harmonicmean {
            pafalignment,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = harmonicestimate(pafalignment).unwrap();
                println!("The harmonic analysis has been written: {:?}", command);
            });
        }
        Commands::PangenomeMatcher {
            pafgenome1,
            pafgenome2,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = paf_alignments(pafgenome1, pafgenome2).unwrap();
                println!(
                    "The comparative pangenome files have been written:{:?}",
                    command
                );
            });
        }
        Commands::PanArc {
            pafalignment,
            fastafile,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = metagenome_annotate(pafalignment, fastafile).unwrap();
                println!("The panarc has been completed: {:?}", command);
            });
        }
        Commands::Snatcher {
            pafalignment,
            queryfasta,
            referencefasta,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = snatcherextract(pafalignment, queryfasta, referencefasta).unwrap();
                println!("The reference snatcher has been completed:{:?}", command);
            });
        }
        Commands::PrecomputedPaf {
            graphalignment,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = precomputealignments(graphalignment).unwrap();
                println!("The precomputed alignment have been written:{:?}", command);
            });
        }
        Commands::PrecomputeCDS {
            pathfileminiprot,
            readsfasta,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = computecds(pathfileminiprot, readsfasta).unwrap();
                println!(
                    "The cds regions from the precomputed alignment have been written:{:?}",
                    command
                );
            });
        }
        Commands::Graph { graph, thread } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(||{ let command = graph_args_segment(graph).unwrap();
            for i in command.iter() {
                println!(
            "Results have been written:\nnumber_segment:{}\nnumber_links:{}\nnumber_arc:{}\nnumber_rank:{}\ntotal_segment_length:{}\naverage_segment_length:{}\nsum_0_segment_length:{}\n",
            i.number_segment,
            i.number_links,
            i.number_arc,
            i.max_rank,
            i.total_segment_length,
            i.average_segment_length,
            i.sum_0_segment_length
        );
            }
           });
        }
        Commands::PangenomeBed { graph, thread } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = graph_bed_segment(graph).unwrap();
                for i in command.iter() {
                    println!(
                        "{:?}\t{:?}\t{:?}\t{:?}\t{:?}",
                        i.tag, i.start, i.end, i.oritag, i.rank
                    );
                }
            });
        }
        Commands::IntergenicNoncoding {
            graph,
            readsfasta,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = computeintergenic(graph, readsfasta).unwrap();
                println!(
                    "The intergenic regions have been extracted for non-coding annotations:{:?}",
                    command
                );
            });
        }
        Commands::PanReadsDatabase { fastafile, thread } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let command = readsdatabase(fastafile).unwrap();
                println!("The reads database has been generated: {:?}", command);
            });
        }
        Commands::BedtoolAncestral {
            alignment,
            fastafile,
            threshold,
            pathprank,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let output =
                    pangenome_longest_alignment(alignment, fastafile, *threshold, pathprank)
                        .unwrap();
                println!("Results have been written:{}", output);
            });
        }
        Commands::VcfAanalyze {
            vcf1,
            vcf2,
            fasta,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let output = vcf_compare(vcf1, vcf2, fasta).unwrap();
                println!("Results have been written:{}", output);
            });
        }
        Commands::PangenomeSingleMerge {
            bed1,
            fasta,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let output = pangenome_merge(bed1, fasta).unwrap();
                println!("Results have been written:{}", output);
            });
        }
        Commands::MultiBedtoolsAncestral {
            alignment,
            fastafile,
            pathprank,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let output = pangenomemultialignment(alignment, fastafile, pathprank).unwrap();
                println!("Results have been written:{}", output);
            });
        }
        Commands::Pafnet {
            pafalignment,
            mode,
            colors,
            prefixchar,
            thread,
        } => {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(thread.parse::<usize>().unwrap())
                .build()
                .unwrap();
            pool.install(|| {
                let output = pafnet_convert(pafalignment, mode, colors, prefixchar).unwrap();
                println!("{}", output);
            });
        }
    }
}
