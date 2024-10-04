use std::str;
use clap::Parser;

mod input;
use input::read_input_file;

mod extremities;
use extremities::identify_good_kmers;

mod read_graph;
use read_graph::build_sequences;

mod filter_output;
use filter_output::filter_output_sequences;

mod utils;


#[derive(Parser, Debug)]
#[command(author = None, version, about = None, long_about = None)]
struct Args {
    /// input SKA2 file [required]
    #[arg(short = 'i', long)]
    skf_file: String,

    /// ignore indels occuring after homopolymers of size n+
    #[arg(short = 'n', long)]
    n_poly: Option<u32>,
    
    /// maximum fraction of unknown states per position 
    #[arg(short = 'm', long, default_value_t = 0.2)]
    max_missing: f32,

    /// name of output files
    #[arg(short = 'o', long, default_value_t = ("skalo").to_string())]
    output_name: String,

}



fn main() {
    println!("\n      skalo v{}     \n", env!("CARGO_PKG_VERSION"));
    
    // get command line arguments
    let args = Args::parse();
    
    // read input file
    let (len_kmer, l_sample_names, all_kmers, index_map) = read_input_file(&args.skf_file);

    // identify 'good' kmers in De Bruijn graph
    let (start_kmers, end_kmers) = identify_good_kmers(len_kmer, &all_kmers, &index_map);

    // build sequences
    let sequences = build_sequences(len_kmer, &all_kmers, &start_kmers, &end_kmers, &index_map);
    
    // filter and output sequences
    filter_output_sequences(sequences, len_kmer, l_sample_names.clone(), args.n_poly, args.max_missing, &args.output_name, &args.skf_file);
    
}
