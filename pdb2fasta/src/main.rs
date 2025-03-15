use once_cell::sync::Lazy;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use regex::Regex;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use clap::Parser;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(long)]
    input: String,

    #[arg(long)]
    output: String,

    #[arg(long, default_value_t = 4)]
    threads: usize,
}

static RES_MAP: Lazy<HashMap<&'static str, char>> = Lazy::new(|| {
    let mut m = HashMap::new();
    m.insert("ALA", 'A');
    m.insert("ARG", 'R');
    m.insert("ASN", 'N');
    m.insert("ASP", 'D');
    m.insert("CYS", 'C');
    m.insert("GLU", 'E');
    m.insert("GLN", 'Q');
    m.insert("GLY", 'G');
    m.insert("HIS", 'H');
    m.insert("ILE", 'I');
    m.insert("LEU", 'L');
    m.insert("LYS", 'K');
    m.insert("MET", 'M');
    m.insert("PHE", 'F');
    m.insert("PRO", 'P');
    m.insert("SER", 'S');
    m.insert("THR", 'T');
    m.insert("TRP", 'W');
    m.insert("TYR", 'Y');
    m.insert("VAL", 'V');
    m
});

fn process_pdb_files(pdb_paths: Vec<String>, output_fasta: &str) {
    let aa_map = &*RES_MAP;

    let sequences: Vec<(String, String)> = pdb_paths
        .par_iter()
        .filter_map(|pdb_path| extract_sequence_from_pdb(pdb_path, aa_map).ok())
        .collect();

    write_fasta_file(output_fasta, &sequences);
}

fn extract_sequence_from_pdb(
    pdb_path: &str,
    aa_map: &HashMap<&str, char>,
) -> Result<(String, String), std::io::Error> {
    let content = fs::read_to_string(pdb_path)?;
    let pdb_id = extract_pdb_id(pdb_path);
    let mut sequence = String::new();

    for line in content.lines() {
        if line.starts_with("ATOM") && line.contains(" CA ") {
            if line.len() >= 20 {
                let residue = line[17..20].trim();
                if let Some(&aa) = aa_map.get(residue) {
                    sequence.push(aa);
                }
            }
        }
    }
    Ok((pdb_id, sequence))
}

fn extract_pdb_id(pdb_path: &str) -> String {
    let re = Regex::new(r"\(AF-\S+-F1\)").unwrap();

    // Extract the PDB ID from the file path
    let path = Path::new(pdb_path);
    let file_stem = path.file_stem().unwrap().to_str().unwrap();
    let pdb_id = re.replace_all(file_stem, "");

    pdb_id.to_string()

}

fn write_fasta_file(output_path: &str, sequences: &[(String, String)]) {
    let file = File::create(output_path).expect("Could not create FASTA file");
    let mut writer = BufWriter::new(file);

    for (id, seq) in sequences {
        writeln!(writer, ">{}\n{}", id, seq).expect("Error writing FASTA");
    }
}

fn main() {
    let args = Args::parse();

    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .expect("Failed to build the global thread pool");

    let file = File::open(&args.input)
        .unwrap_or_else(|_| panic!("Failed to open input file: {}", args.input));
    let reader = BufReader::new(file);
    let pdb_paths: Vec<String> = reader
        .lines()
        .filter_map(Result::ok)
        .filter(|line| !line.trim().is_empty())
        .collect();

    process_pdb_files(pdb_paths, &args.output);
}