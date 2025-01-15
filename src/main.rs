use std::str;
use clap::Parser;
use std::path::PathBuf;

mod input;
use input::read_input_file;

mod extremities;
use extremities::identify_good_kmers;

mod read_graph;
use read_graph::build_variant_groups;


mod compaction;
mod process_variants;
mod positioning;
mod output;
mod utils;
use crate::utils::{DATA_INFO, CONFIG, DataInfo, Config};


#[derive(Parser, Debug)]
#[command(
    author = None,
    version,
    about = None,
    long_about = None,
    term_width = 200,
    help_template = "

 Usage: {usage}

 input:
   -i, --input-skf      input SKA2 file
   -r, --reference      reference genome for variant positioning

 output:
   -o, --output-name    prefix of output files [default: skalo]
   -m, --missing        max. fraction of missing data [default: 0.2]

 graph traversal:
   -d, --depth          max. depth of recursive paths [default: 4]

 other:
   -n, --indel-kmers    max. number of internal indel k-mers [default: 2]
   -t, --threads        number of threads [default: 1]
.
"
)]
struct Args {
    /// input SKA2 file
    #[arg(short = 'i', long, help_heading = "input")]
    input_skf: String,

    /// reference genome for SNP positioning
    #[arg(short = 'r', long, help_heading = "input")]
    reference: Option<PathBuf>,

    /// prefix of output files
    #[arg(short = 'o', long, default_value_t = ("skalo").to_string(), help_heading = "output")]
    output: String,

    /// maximum fraction of missing data
    #[arg(short = 'm', long, default_value_t = 0.2, help_heading = "output")]
    missing: f32,

    /// maximum depth of recursive paths
    #[arg(short = 'd', long, default_value_t = 4, help_heading = "graph traversal")]
    depth: usize,

    /// maximum number of internal indel k-mers
    #[arg(short = 'n', long, default_value_t = 2, help_heading = "other")]
    indel_kmers: usize,

    /// number of threads
    #[arg(short = 't', long, default_value_t = 1, help_heading = "other")]
    threads: usize,
}



fn main() {
    println!("\n      skalo v{}     \n", env!("CARGO_PKG_VERSION"));
    
    // get command line arguments
    let args = Args::parse();
    
    // initialise the global CONFIG structure
    CONFIG.set(Config {
        input_file: args.input_skf.clone(),
        output_name: args.output.clone(),
        max_missing: args.missing,
        max_depth: args.depth,
        max_indel_kmers: args.indel_kmers,
        nb_threads: args.threads,
        reference_genome: args.reference.clone(),
    }).expect("failed to initialise CONFIG");
    
    // read input file
    let (len_kmer, sample_names, all_kmers, index_map) = read_input_file();
    
    // initialise the global DataInfo structure
    DATA_INFO.set(DataInfo {
        k_graph: len_kmer -1,
        sample_names: sample_names.clone(),
    }).expect("failed to initialise DATA_INFO");    
    
    // identify 'good' kmers in De Bruijn graph
    let (start_kmers, end_kmers) = identify_good_kmers(&all_kmers, &index_map);

    // identify variant groups
    build_variant_groups(all_kmers, start_kmers, end_kmers, index_map);


}
