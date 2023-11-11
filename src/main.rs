use std::str;
use clap::Parser;

mod input;
use input::read_input_file;

mod extremities;
use extremities::identify_good_kmers;

mod read_graph;
use read_graph::build_sequences;

mod filter;
use filter::filter_sequences;

mod output;
use output::create_output_files;

mod utils;


#[derive(Parser, Debug)]
#[command(author = None, version, about = None, long_about = None)]
struct Args {
    /// input SKA tabular file [required]
    #[arg(short = 'i', long)]
    skf_file: String,

    /// minimum jaccard similarity between sample sets [0.8]
    #[arg(short = 's', long, default_value_t = 0.8)]
    min_jaccard: f32,
    
    /// ignore indels occuring after homopolymers of size n+ [6]
    #[arg(short = 'n', long, default_value_t = 6)]
    n_poly: u32,
    
    /// maximum fraction of unknown states per position  [0.2]
    #[arg(short = 'm', long, default_value_t = 0.2)]
    max_missing: f32,

    /// name of output files []
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
    let good_kmers = identify_good_kmers(len_kmer, &all_kmers, &index_map);
    
    // build sequences
    let sequences = build_sequences(args.min_jaccard, len_kmer, &all_kmers, &good_kmers, &index_map);
    
    // filter sequences
    let filtered_sequences = filter_sequences(sequences, len_kmer, l_sample_names.clone(), args.n_poly, args.max_missing);
    
    // output results
    create_output_files(filtered_sequences, len_kmer, l_sample_names, &args.output_name);
    
}

