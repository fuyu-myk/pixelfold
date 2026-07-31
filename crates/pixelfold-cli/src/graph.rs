//! Writes the residue interaction network for another tool to read.
//!
//! Three shapes: node-link JSON that a browser or a script parses directly,
//! GraphML that Cytoscape and Gephi import, and a flat edge list for a
//! dataframe.

use std::io::Write;

use anyhow::Result;
use clap::ValueEnum;
use pixelfold_core::Network;
use pixelfold_core::rin::{Edge, Node, RinAnalysis};
use serde::Serialize;
use serde_json::{Value, json};

/// How many hub residues the analysis report lists.
const TOP_HUBS: usize = 10;
/// How many connected components the analysis report lists before summarising.
const TOP_COMPONENTS: usize = 10;

/// How to write a network.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum GraphFormat {
    /// Node-link JSON: `{ "nodes": [...], "edges": [...] }` (default)
    #[default]
    Json,
    /// GraphML, for Cytoscape and Gephi
    Graphml,
    /// The edge list as tab-separated values
    Tsv,
}

/// Write `network` to `out` in `format`.
pub fn write<W: Write + ?Sized>(out: &mut W, network: &Network, format: GraphFormat) -> Result<()> {
    match format {
        GraphFormat::Json => write_json(out, network),
        GraphFormat::Graphml => write_graphml(out, network),
        GraphFormat::Tsv => write_tsv(out, network),
    }
}

/// Write a structural analysis of `network`: its connected components,
/// the residues that bridge it (highest betweenness), and the cut residues
/// whose removal would split it.
pub fn write_analysis<W: Write + ?Sized>(
    out: &mut W,
    network: &Network,
    analysis: &RinAnalysis,
) -> Result<()> {
    writeln!(
        out,
        "Residue interaction network: {} residues, {} edges, {} components",
        network.nodes.len(),
        network.edges.len(),
        analysis.components.len()
    )?;

    writeln!(out, "\nConnected components (largest first):")?;
    if analysis.components.is_empty() {
        writeln!(out, "  (none)")?;
    } else {
        for (rank, component) in analysis.components.iter().take(TOP_COMPONENTS).enumerate() {
            writeln!(out, "  {}. {} residues", rank + 1, component.len())?;
        }

        let extra = analysis.components.len().saturating_sub(TOP_COMPONENTS);
        if extra > 0 {
            writeln!(out, "  ... and {extra} more")?;
        }
    }

    writeln!(out, "\nTop hubs by betweenness centrality:")?;
    let mut ranked: Vec<usize> = (0..network.nodes.len()).collect();
    ranked.sort_by(|&a, &b| {
        analysis.betweenness[b]
            .total_cmp(&analysis.betweenness[a])
            .then(network.nodes[b].degree.cmp(&network.nodes[a].degree))
    });
    writeln!(
        out,
        "  {:<10} {:<5} {:<3} {:>6} {:>12}",
        "residue", "resn", "ss", "degree", "between"
    )?;
    for &i in ranked.iter().take(TOP_HUBS) {
        let node = &network.nodes[i];
        writeln!(
            out,
            "  {:<10} {:<5} {:<3} {:>6} {:>12.2}",
            node.id, node.resn, node.ss, node.degree, analysis.betweenness[i]
        )?;
    }

    writeln!(
        out,
        "\nArticulation points ({} cut residues):",
        analysis.articulation_points.len()
    )?;
    if analysis.articulation_points.is_empty() {
        writeln!(out, "  (none)")?;
    } else {
        let residues: Vec<String> = analysis
            .articulation_points
            .iter()
            .map(|&i| format!("{} {}", network.nodes[i].id, network.nodes[i].resn))
            .collect();
        writeln!(out, "  {}", residues.join(", "))?;
    }

    Ok(())
}

/// Write the shortest chain of interactions between two residues, one residue
/// per line, or a note when they are not connected.
pub fn write_path<W: Write + ?Sized>(
    out: &mut W,
    network: &Network,
    from_id: &str,
    to_id: &str,
    route: Option<&[usize]>,
) -> Result<()> {
    match route {
        None => writeln!(
            out,
            "No path between {from_id} and {to_id}: they are in different components"
        )?,
        Some(path) => {
            writeln!(
                out,
                "Shortest path {from_id} \u{2192} {to_id}: {} steps, {} residues",
                path.len().saturating_sub(1),
                path.len()
            )?;
            for &node in path {
                let node = &network.nodes[node];
                writeln!(out, "  {} {} [{}]", node.id, node.resn, node.ss)?;
            }
        }
    }

    Ok(())
}

/// Write a MolViewSpec (.mvsj) scene: the structure loaded from `url`, its
/// polymer drawn as a grey cartoon, and the network's residues picked out in
/// ball-and-stick. Mol* opens the file and reproduces the scene, so the network
/// selection becomes a shareable, diffable 3D view. `format` is the parse format
/// of the structure at `url` (`mmcif` or `pdb`).
pub fn write_mvsj<W: Write + ?Sized>(
    out: &mut W,
    network: &Network,
    url: &str,
    format: &str,
    title: &str,
) -> Result<()> {
    let residues: Vec<Value> = network
        .nodes
        .iter()
        .map(|node| {
            let mut expr = serde_json::Map::new();
            expr.insert("auth_asym_id".to_string(), json!(node.chain));
            expr.insert("auth_seq_id".to_string(), json!(node.resi));
            if let Some(icode) = node.icode {
                expr.insert("pdbx_PDB_ins_code".to_string(), json!(icode.to_string()));
            }
            Value::Object(expr)
        })
        .collect();

    let mut components = vec![component(json!("polymer"), "cartoon", "#8a8a99")];
    if !residues.is_empty() {
        components.push(component(json!(residues), "ball_and_stick", "#ffb454"));
    }

    let scene = json!({
        "metadata": { "version": "1", "title": title },
        "root": { "kind": "root", "children": [
            { "kind": "download", "params": { "url": url }, "children": [
                { "kind": "parse", "params": { "format": format }, "children": [
                    { "kind": "structure", "params": { "type": "model" }, "children": components }
                ]}
            ]}
        ]}
    });

    serde_json::to_writer_pretty(&mut *out, &scene)?;
    writeln!(out)?;

    Ok(())
}

/// A `component` node over `selector` (a preset string or an array of residue
/// expressions), drawn with one representation in one colour.
fn component(selector: Value, representation: &str, color: &str) -> Value {
    json!({
        "kind": "component",
        "params": { "selector": selector },
        "children": [{
            "kind": "representation",
            "params": { "type": representation },
            "children": [{ "kind": "color", "params": { "color": color } }]
        }]
    })
}

/// A distance rounded to the precision the coordinates carry.
fn distance(value: f32) -> f32 {
    (value * 100.0).round() / 100.0
}

#[derive(Serialize)]
struct JsonNetwork<'a> {
    nodes: Vec<JsonNode<'a>>,
    edges: Vec<JsonEdge<'a>>,
}

#[derive(Serialize)]
struct JsonNode<'a> {
    id: &'a str,
    chain: &'a str,
    resi: u32,
    icode: Option<char>,
    resn: &'a str,
    ss: char,
    degree: usize,
}

#[derive(Serialize)]
struct JsonEdge<'a> {
    source: &'a str,
    target: &'a str,
    kind: &'a str,
    count: usize,
    min_distance: f32,
    interchain: bool,
}

fn json_node(node: &Node) -> JsonNode<'_> {
    JsonNode {
        id: &node.id,
        chain: &node.chain,
        resi: node.resi,
        icode: node.icode,
        resn: &node.resn,
        ss: node.ss,
        degree: node.degree,
    }
}

fn json_edge(edge: &Edge) -> JsonEdge<'_> {
    JsonEdge {
        source: &edge.source,
        target: &edge.target,
        kind: edge.kind.label(),
        count: edge.count,
        min_distance: distance(edge.min_distance),
        interchain: edge.interchain,
    }
}

fn write_json<W: Write + ?Sized>(out: &mut W, network: &Network) -> Result<()> {
    let json = JsonNetwork {
        nodes: network.nodes.iter().map(json_node).collect(),
        edges: network.edges.iter().map(json_edge).collect(),
    };
    serde_json::to_writer_pretty(&mut *out, &json)?;
    writeln!(out)?;

    Ok(())
}

/// The edge list, one row per edge. Node attributes have no place in a single
/// flat table, so this carries the edges alone; JSON or GraphML carry the whole
/// graph.
fn write_tsv<W: Write + ?Sized>(out: &mut W, network: &Network) -> Result<()> {
    writeln!(out, "source\ttarget\tkind\tcount\tmin_distance\tinterchain")?;
    for edge in &network.edges {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.2}\t{}",
            edge.source,
            edge.target,
            edge.kind.label(),
            edge.count,
            distance(edge.min_distance),
            edge.interchain,
        )?;
    }

    Ok(())
}

/// The GraphML attribute keys, one per node and edge field.
const KEYS: &str = concat!(
    r#"  <key id="v_chain" for="node" attr.name="chain" attr.type="string"/>"#,
    "\n",
    r#"  <key id="v_resi" for="node" attr.name="resi" attr.type="int"/>"#,
    "\n",
    r#"  <key id="v_icode" for="node" attr.name="icode" attr.type="string"/>"#,
    "\n",
    r#"  <key id="v_resn" for="node" attr.name="resn" attr.type="string"/>"#,
    "\n",
    r#"  <key id="v_ss" for="node" attr.name="ss" attr.type="string"/>"#,
    "\n",
    r#"  <key id="v_degree" for="node" attr.name="degree" attr.type="int"/>"#,
    "\n",
    r#"  <key id="e_kind" for="edge" attr.name="kind" attr.type="string"/>"#,
    "\n",
    r#"  <key id="e_count" for="edge" attr.name="count" attr.type="int"/>"#,
    "\n",
    r#"  <key id="e_distance" for="edge" attr.name="min_distance" attr.type="double"/>"#,
    "\n",
    r#"  <key id="e_interchain" for="edge" attr.name="interchain" attr.type="boolean"/>"#,
);

fn write_graphml<W: Write + ?Sized>(out: &mut W, network: &Network) -> Result<()> {
    writeln!(out, r#"<?xml version="1.0" encoding="UTF-8"?>"#)?;
    writeln!(
        out,
        r#"<graphml xmlns="http://graphml.graphdrawing.org/xmlns">"#
    )?;
    writeln!(out, "{KEYS}")?;
    writeln!(out, r#"  <graph edgedefault="undirected">"#)?;

    for node in &network.nodes {
        writeln!(out, r#"    <node id="{}">"#, escape(&node.id))?;
        data(out, "v_chain", &escape(&node.chain))?;
        data(out, "v_resi", &node.resi.to_string())?;
        if let Some(icode) = node.icode {
            data(out, "v_icode", &escape(&icode.to_string()))?;
        }
        data(out, "v_resn", &escape(&node.resn))?;
        data(out, "v_ss", &node.ss.to_string())?;
        data(out, "v_degree", &node.degree.to_string())?;
        writeln!(out, "    </node>")?;
    }

    for edge in &network.edges {
        writeln!(
            out,
            r#"    <edge source="{}" target="{}">"#,
            escape(&edge.source),
            escape(&edge.target),
        )?;
        data(out, "e_kind", edge.kind.label())?;
        data(out, "e_count", &edge.count.to_string())?;
        data(
            out,
            "e_distance",
            &format!("{:.2}", distance(edge.min_distance)),
        )?;
        data(out, "e_interchain", &edge.interchain.to_string())?;
        writeln!(out, "    </edge>")?;
    }

    writeln!(out, "  </graph>")?;
    writeln!(out, "</graphml>")?;

    Ok(())
}

fn data<W: Write + ?Sized>(out: &mut W, key: &str, value: &str) -> Result<()> {
    writeln!(out, r#"      <data key="{key}">{value}</data>"#)?;

    Ok(())
}

/// Escape the five XML metacharacters, so an odd chain id or residue name cannot
/// break the document.
fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for ch in text.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(ch),
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use pixelfold_core::InteractionKind;
    use pixelfold_core::rin::{Edge, Network, Node};

    fn network() -> Network {
        Network {
            nodes: vec![
                Node {
                    id: "A/1".into(),
                    chain: "A".into(),
                    resi: 1,
                    icode: None,
                    resn: "LYS".into(),
                    ss: 'H',
                    degree: 1,
                },
                Node {
                    id: "A/2".into(),
                    chain: "A".into(),
                    resi: 2,
                    icode: Some('B'),
                    resn: "GLU".into(),
                    ss: 'E',
                    degree: 1,
                },
            ],
            edges: vec![Edge {
                source: "A/1".into(),
                target: "A/2".into(),
                kind: InteractionKind::SaltBridge,
                count: 1,
                min_distance: 3.456,
                interchain: false,
            }],
        }
    }

    fn rendered(format: GraphFormat) -> String {
        let mut out = Vec::new();
        write(&mut out, &network(), format).expect("writing to a vec cannot fail");

        String::from_utf8(out).expect("utf8")
    }

    #[test]
    fn json_is_node_link_and_rounds_the_distance() {
        let parsed: serde_json::Value =
            serde_json::from_str(&rendered(GraphFormat::Json)).expect("valid json");

        assert_eq!(parsed["nodes"][1]["icode"], "B");
        assert_eq!(parsed["edges"][0]["kind"], "salt-bridge");
        assert_eq!(parsed["edges"][0]["min_distance"], 3.46);
        assert_eq!(parsed["edges"][0]["interchain"], false);
    }

    #[test]
    fn the_tsv_edge_list_has_a_header_and_a_row_per_edge() {
        let text = rendered(GraphFormat::Tsv);
        let lines: Vec<&str> = text.lines().collect();

        assert_eq!(
            lines[0],
            "source\ttarget\tkind\tcount\tmin_distance\tinterchain"
        );
        assert_eq!(lines[1], "A/1\tA/2\tsalt-bridge\t1\t3.46\tfalse");
    }

    #[test]
    fn graphml_is_well_formed_and_declares_its_keys() {
        let text = rendered(GraphFormat::Graphml);

        assert!(text.starts_with(r#"<?xml version="1.0" encoding="UTF-8"?>"#));
        assert!(
            text.contains(r#"<key id="e_kind" for="edge" attr.name="kind" attr.type="string"/>"#)
        );
        assert!(text.contains(r#"<node id="A/1">"#));
        assert!(text.contains(r#"<edge source="A/1" target="A/2">"#));
        assert!(text.contains("<data key=\"e_kind\">salt-bridge</data>"));
        // The insertion code reaches the node.
        assert!(text.contains("<data key=\"v_icode\">B</data>"));
        // Every open tag is closed.
        assert_eq!(
            text.matches("<node ").count(),
            text.matches("</node>").count()
        );
        assert_eq!(
            text.matches("<edge ").count(),
            text.matches("</edge>").count()
        );
    }

    #[test]
    fn xml_metacharacters_in_a_value_are_escaped() {
        let mut net = network();
        net.nodes[0].resn = "A<&B".into();
        let mut out = Vec::new();
        write(&mut out, &net, GraphFormat::Graphml).expect("write");
        let text = String::from_utf8(out).expect("utf8");

        assert!(text.contains("A&lt;&amp;B"));
        assert!(!text.contains("A<&B"));
    }
}
