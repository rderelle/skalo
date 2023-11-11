use hashbrown::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;


pub fn create_output_files(filtered_sequences: HashMap<String, Vec<(String, Vec<String>)>>, len_kmer: usize, l_sample_names: Vec<String>, output_name: &str) {
    // Prepare output files
    let mut out_1 = File::create(format!("{}_seq_groups.fas", output_name)).expect("Failed to create file");
    let mut out_2 = File::create(format!("{}_variants.tsv", output_name)).expect("Failed to create file");
    
    // prepare header of variant file
    let index_name: Vec<String> = l_sample_names.iter().enumerate().map(|(index, value)| format!("{}:{}", index, value)).collect();
    let formatted_string = index_name.join(", ");
    out_2.write_all(format!("#samples: {}\n", formatted_string).as_bytes())
        .expect("Failed to write to file");
    out_2.write_all(b"#pos_ali\tnb_states\ttype\tfirst_kmer\tvariants\tlast_kmer\tsamples\n")
        .expect("Failed to write to file");

    // Prepare variables for binary alignment
    let mut d_matrix_seq: HashMap<String, Vec<String>> = HashMap::new();
    let mut d_samples: HashMap<String, String> = HashMap::new();
    for (i, sample) in l_sample_names.iter().enumerate() {
        d_samples.insert(i.to_string(), sample.clone()); // map sample ID to sample full name
        d_matrix_seq.insert(sample.clone(), vec!["-".to_string(); filtered_sequences.len()]); // alignment contains only "-", to be replaced by 0 and 1
    }

    // Browse final sequences
    let mut position = 0;
    for ll in filtered_sequences.values() {
        // Output seq_groups
        for l in ll {
            let str_samples = l.1.join("|");
            out_1.write_all(format!(">{}_{}\n{}\n", position, str_samples, l.0).as_bytes())
                .expect("Failed to write to file");
        }

        // Update binary alignment
        // '-' = missing data
        // '?' = both states (= 'unknown')
        let mut state = 0;
        //  we get all sequences 1 by 1 with a maximum of 2 sequences
        //  (i.e., the two sequences corresponding to the most samples)
        //  -> if third or more sequences present, these will be ignored and corresponding samples will have an unknown ('-') state
        for l in ll.iter().take(2) {
            // collect sample names for this sequence and update its vector in d_matrix_seq
            for sample_id in &l.1 {
                let full_name = d_samples.get(sample_id).unwrap();
                let sample_vec = d_matrix_seq.get_mut(full_name).unwrap();
                if sample_vec[position] == "-" {
                    sample_vec[position] = state.to_string();
                } else {
                    sample_vec[position] = "?".to_string();
                }
            }
            // Update state
            state += 1;
        }
        
        // Extract first kmer (using first sequence)
        let first_kmer = ll[0].0[..len_kmer - 1].to_string();
        
        // Collect all sequences without the first kmer
        let l_reduced_seq: Vec<&str> = ll.iter().map(|l| &l.0[len_kmer - 1..]).collect();

        // Extract last kmer
        let mut identical = true;
        let mut n_nucl = 0;
        while identical {
            n_nucl += 1;
            
            let mut all_ends: HashSet<String> = HashSet::new();
            // extract last n nucleotide from each seq
            for seq in &l_reduced_seq {
                if n_nucl > seq.len() {
                    identical = false;
                } else {                    
                    let last_n_chars: Vec<String> = seq.chars().rev().take(n_nucl).map(|c| c.to_string()).collect();
                    let concatenated_last_chars: String = last_n_chars.into_iter().rev().collect();
                    all_ends.insert(concatenated_last_chars.clone());
                }
            }

            if all_ends.len() > 1 {
                identical = false;
            }
        }
        n_nucl -= 1;
        
        // extract last kmer from first sequence
        let pos_end = l_reduced_seq[0].len() - n_nucl;
        let last_kmer = l_reduced_seq[0][pos_end..].to_string();
        
        // remove last kmer from reduced sequences -> only middle base left
        let mut type_variant = "multi";
        let mut l_middles: Vec<String> = Vec::new();
        for seq in &l_reduced_seq {
            let pos_end = seq.len() - n_nucl;
            let mut middle_bases = &seq[..pos_end];
            if middle_bases == "" {
                middle_bases = ".";
                type_variant = "indel";
            }
            l_middles.push(middle_bases.to_string());
        }
        
        // get lists of samples
        let mut l_samples = Vec::new();
        for (_, vect) in ll {
            l_samples.push(vect.join(","));
        }

        // Save variants information
        out_2.write_all(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            position, l_middles.len(), type_variant, first_kmer,
            l_middles.join(" / "), last_kmer, l_samples.join(" / ")).as_bytes())
            .expect("Failed to write to file");

        // Update position
        position += 1;

    }
    
    // Close output files
    out_1.flush().expect("Failed to flush file");
    out_2.flush().expect("Failed to flush file");

    // Output alignment
    let mut out_3 = File::create(format!("{}_binary_ali.fas", output_name)).expect("Failed to create file");
    for (sample, l_seq) in &d_matrix_seq {
        out_3.write_all(format!(">{}\n{}\n", sample, l_seq.join("")).as_bytes())
            .expect("Failed to write to file");
    }
    out_3.flush().expect("Failed to flush file");
}


