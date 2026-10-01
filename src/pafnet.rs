use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;
use rustc_hash::FxHashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

/*
Gaurav Sablok,
gsablok@proton.me

Ported from pafnet-fast (a petgraph-backed reimplementation of ekg/pafnet):
projects a PAF alignment file into network/graph output formats (Pajek
.net, plain edgelist, node list, id-rewritten PAF, or GEXF for Gephi) in
a single pass over the input, keeping alignment-quality columns (matching
bases, block length, mapping quality) as edge weights.
*/

/// Edge attributes carried over from the PAF record's alignment columns.
#[derive(Clone, Copy, Debug, Default)]
struct AlnWeight {
    matches: u32,
    block_len: u32,
    mapq: u8,
}

impl AlnWeight {
    /// Fraction of the alignment block that is an exact match. Guards
    /// against divide-by-zero on malformed/zero-length records.
    fn identity(&self) -> f64 {
        if self.block_len == 0 {
            0.0
        } else {
            self.matches as f64 / self.block_len as f64
        }
    }
}

/// The graph: node weight is the sequence name, edge weight is the
/// alignment quality info pulled from the PAF record.
type PafGraph = DiGraph<String, AlnWeight>;

/// Look up (or create) the node for `name`, reusing an owned allocation
/// only on first sight of a given sequence.
#[inline]
fn intern(graph: &mut PafGraph, map: &mut FxHashMap<String, NodeIndex>, name: &str) -> NodeIndex {
    if let Some(&idx) = map.get(name) {
        idx
    } else {
        let idx = graph.add_node(name.to_string());
        map.insert(name.to_string(), idx);
        idx
    }
}

/// Parse the whole PAF file in one pass, building the graph directly.
/// Returns the graph plus the name -> NodeIndex map (kept around so
/// rewrite mode can reuse it on a second, cheap pass instead of
/// rebuilding an index).
fn build_graph(path: &str) -> Result<(PafGraph, FxHashMap<String, NodeIndex>), Box<dyn Error>> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(1 << 20, file);

    let mut graph: PafGraph = DiGraph::new();
    let mut map: FxHashMap<String, NodeIndex> = FxHashMap::default();

    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        let line = line.trim_end_matches(['\n', '\r']);
        if line.is_empty() {
            continue;
        }

        let mut f = line.split('\t');
        let q = match f.next() {
            Some(v) => v,
            None => continue,
        };
        f.next(); // qlen
        f.next(); // qstart
        f.next(); // qend
        f.next(); // strand
        let t = match f.next() {
            Some(v) => v,
            None => continue,
        };
        f.next(); // tlen
        f.next(); // tstart
        f.next(); // tend
        let matches: u32 = f.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let block_len: u32 = f.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let mapq: u8 = f.next().and_then(|s| s.parse().ok()).unwrap_or(0);

        let qi = intern(&mut graph, &mut map, q);
        let ti = intern(&mut graph, &mut map, t);
        graph.add_edge(
            qi,
            ti,
            AlnWeight {
                matches,
                block_len,
                mapq,
            },
        );
    }

    Ok((graph, map))
}

fn write_pajek_net(graph: &PafGraph, out: &mut impl Write) -> Result<(), Box<dyn Error>> {
    writeln!(out, "*Vertices {}", graph.node_count())?;
    for idx in graph.node_indices() {
        // Pajek ids are 1-based.
        writeln!(out, "{} \"{}\"", idx.index() + 1, graph[idx])?;
    }
    writeln!(out, "*arcs")?;
    for e in graph.edge_references() {
        writeln!(
            out,
            "{} {} {:.4}",
            e.source().index() + 1,
            e.target().index() + 1,
            e.weight().identity()
        )?;
    }
    Ok(())
}

fn write_edgelist(graph: &PafGraph, out: &mut impl Write) -> Result<(), Box<dyn Error>> {
    for e in graph.edge_references() {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            e.source().index() + 1,
            e.target().index() + 1,
            e.weight().matches,
            e.weight().block_len,
            e.weight().mapq
        )?;
    }
    Ok(())
}

fn write_nodelist(graph: &PafGraph, out: &mut impl Write) -> Result<(), Box<dyn Error>> {
    for idx in graph.node_indices() {
        writeln!(out, "{}\t{}", idx.index() + 1, graph[idx])?;
    }
    Ok(())
}

fn write_gexf(
    graph: &PafGraph,
    out: &mut impl Write,
    colors: &FxHashMap<String, (u16, u16, u16)>,
    prefix_char: &str,
) -> Result<(), Box<dyn Error>> {
    writeln!(out, r#"<?xml version="1.0" encoding="UTF-8"?>"#)?;
    writeln!(
        out,
        r#"<gexf xmlns="http://www.gexf.net/1.2draft" xmlns:viz="http://www.gexf.net/1.1draft/viz" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd" version="1.2">"#
    )?;
    writeln!(out, "<graph defaultedgetype=\"directed\">")?;
    writeln!(out, "<nodes>")?;
    for idx in graph.node_indices() {
        let name = &graph[idx];
        let group = if prefix_char.is_empty() {
            name.as_str()
        } else {
            name.split(prefix_char).next().unwrap_or(name.as_str())
        };
        let (r, g, b) = colors.get(group).copied().unwrap_or((127, 127, 127));
        write!(out, r#"<node id="{}" label="{}">"#, idx.index() + 1, name)?;
        write!(out, r#"<viz:color r="{}" g="{}" b="{}" a="1.0"/>"#, r, g, b)?;
        writeln!(out, "</node>")?;
    }
    writeln!(out, "</nodes>")?;
    writeln!(out, "<edges>")?;
    for (i, e) in graph.edge_references().enumerate() {
        writeln!(
            out,
            r#"<edge id="{}" source="{}" target="{}" weight="{:.4}"/>"#,
            i + 1,
            e.source().index() + 1,
            e.target().index() + 1,
            e.weight().identity()
        )?;
    }
    writeln!(out, "</edges>")?;
    writeln!(out, "</graph>")?;
    writeln!(out, "</gexf>")?;
    Ok(())
}

/// Second, cheap pass: reuse the already-built name -> id map to rewrite
/// query/target names as integer ids. Only done for rewrite mode, so the
/// common cases (net/edgelist/nodelist/gexf) stay single-pass.
fn rewrite_paf(
    path: &str,
    map: &FxHashMap<String, NodeIndex>,
    out: &mut impl Write,
) -> Result<(), Box<dyn Error>> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let mut line = String::new();
    let mut buf = String::with_capacity(256);
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        let line_trimmed = line.trim_end_matches(['\n', '\r']);
        if line_trimmed.is_empty() {
            continue;
        }
        buf.clear();
        for (i, field) in line_trimmed.split('\t').enumerate() {
            if i > 0 {
                buf.push('\t');
            }
            if i == 0 || i == 5 {
                if let Some(idx) = map.get(field) {
                    buf.push_str(&(idx.index() + 1).to_string());
                    continue;
                }
            }
            buf.push_str(field);
        }
        writeln!(out, "{}", buf)?;
    }
    Ok(())
}

fn load_colors(path: &str) -> Result<FxHashMap<String, (u16, u16, u16)>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut colors = FxHashMap::default();
    for line in reader.lines() {
        let line = line?;
        let mut it = line.split_ascii_whitespace();
        let name = match it.next() {
            Some(v) => v.to_string(),
            None => continue,
        };
        let r: u16 = it.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let g: u16 = it.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let b: u16 = it.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        colors.insert(name, (r, g, b));
    }
    Ok(colors)
}

/// Convert a PAF alignment file into a network/graph format.
///
/// `mode` selects the output: "net" (Pajek .net), "edgelist", "nodelist",
/// "rewrite" (PAF rewritten with integer ids), or "gexf" (for Gephi).
/// `colors` is the path to a "group R G B" colors file used only by gexf
/// mode (pass "" for none). `prefixchar` is the character sequence names
/// are split on to obtain a group prefix for gexf coloring (pass "" for
/// none, i.e. the whole name is the group).
pub fn pafnet_convert(
    pafalignment: &str,
    mode: &str,
    colors: &str,
    prefixchar: &str,
) -> Result<String, Box<dyn Error>> {
    let (graph, map) = build_graph(pafalignment)?;

    let outname = match mode {
        "net" => "pafnet.net",
        "edgelist" => "pafnet.edges",
        "nodelist" => "pafnet.nodes",
        "rewrite" => "pafnet.ids.paf",
        "gexf" => "pafnet.gexf",
        other => return Err(format!("unknown pafnet mode: {other}").into()),
    };

    let mut out = BufWriter::with_capacity(1 << 20, File::create(outname)?);

    match mode {
        "net" => write_pajek_net(&graph, &mut out)?,
        "edgelist" => write_edgelist(&graph, &mut out)?,
        "nodelist" => write_nodelist(&graph, &mut out)?,
        "rewrite" => rewrite_paf(pafalignment, &map, &mut out)?,
        "gexf" => {
            let colormap = if colors.is_empty() {
                FxHashMap::default()
            } else {
                load_colors(colors)?
            };
            write_gexf(&graph, &mut out, &colormap, prefixchar)?
        }
        _ => unreachable!(),
    }

    out.flush()?;

    Ok(format!(
        "The paf alignment has been converted to network format: {outname}"
    ))
}
