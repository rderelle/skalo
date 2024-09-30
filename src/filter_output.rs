use hashbrown::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use crate::utils::rev_compl;


pub fn filter_output_sequences(sequences: HashMap<String, HashMap<String, Vec<String>>>, len_kmer: usize, all_samples: Vec<String>, n_homopolymer: Option<u32>, max_missing: f32, output_name: &str, input_name: &str) {
    println!(" # filter sequences");
    
    let len_kmer_graph = len_kmer -1;
    
    // prepare output files
    let mut out_1 = File::create(format!("{}_seq_groups.fas", output_name)).expect("Failed to create file");
    let mut out_2 = File::create(format!("{}_variants.tsv", output_name)).expect("Failed to create file");
    
    // write header of TSV variant file
    out_2.write_all(format!("#skalo version {}\n",env!("CARGO_PKG_VERSION")).as_bytes())
        .expect("Failed to write to file");
    out_2.write_all(format!("#parameters: m={}{}\n",max_missing,
        match n_homopolymer {
            Some(n) => format!(" n={}", n),
            None => String::new(),
        }
    ).as_bytes())
        .expect("Failed to write to file");
    out_2.write_all(format!("#input file: {} (k={})\n", input_name, len_kmer).as_bytes())
        .expect("Failed to write to file");

    let index_name: Vec<String> = all_samples.iter().enumerate().map(|(index, value)| format!("{}:{}", index, value)).collect();
    let formatted_string = index_name.join(", ");
    out_2.write_all(format!("#samples: {}\n", formatted_string).as_bytes())
        .expect("Failed to write to file");
    out_2.write_all(b"#pos_ali\tnb_states\ttype\tunclear_insert\tfirst_kmer\tvariants\tlast_kmer\tratio_missing\tsamples\n")
        .expect("Failed to write to file");
    
    // Prepare variables for binary alignment
    let mut binary_seq: HashMap<String, Vec<String>> = HashMap::new();
    let mut d_samples: HashMap<String, String> = HashMap::new();
    for (i, sample) in all_samples.iter().enumerate() {
        binary_seq.insert(sample.clone(), Vec::new());
        d_samples.insert(i.to_string(), sample.clone()); // map sample ID to sample full name
    }
    
    let mut nb_missing = 0;
    let mut nb_homopolymer = 0;
 
    // sequence groups 1 by 1
    let mut position = 0;
    for d_seq in sequences.values() {
        
        // sort variants by number of samples
        let ll_sorted_seq = sort_by_nb_samples(d_seq.clone());
        
        // check the ratio of missing samples
        let ratio_missing = calculate_ratio_missing(all_samples.len(), ll_sorted_seq.clone());
        
        if ratio_missing > max_missing {
            nb_missing += 1;
        } else {
            
            /*
            CHECK IF UNCLEAR INSERT
            */ 
            
            // extract first kmer (using first sequence)
            let first_kmer = ll_sorted_seq[0].0[..len_kmer_graph].to_string();
            
            // collect 'middle-bases' and last k-mer
            let (l_middle_bases, last_kmer) = collect_middle_bases(ll_sorted_seq.clone(), len_kmer_graph);
            
            // collect 'middle-bases' and last k-mer from Rev-Compl sequences
            let ll_rc_seq: Vec<(String, Vec<String>)> = ll_sorted_seq
                .iter()
                .map(|(s, v)| (rev_compl(s), v.clone()))
                .collect();
            let (rc_l_middle_bases, rc_last_kmer) = collect_middle_bases(ll_rc_seq.clone(), len_kmer_graph);
            
            // test if multiple positions for last k-mer
            let multiple_positions = test_multiple_positions(l_middle_bases.clone(), last_kmer.clone());
            let multiple_positions_rc = test_multiple_positions(rc_l_middle_bases, rc_last_kmer);
            
            let mut unclear_insert = "no";
            if multiple_positions || multiple_positions_rc {
                unclear_insert = "yes";
            }
            
            /*
            CHECK IK HOMOPOLYMER
            */ 

            // check if indel in homopolymer
            let mut is_homopolymer = false;
            if let Some(mut n_max) = n_homopolymer {
                // adjust n_homopolymer if higher than k-mer length to avoid the program crashing
                if n_max > len_kmer_graph.try_into().unwrap() {
                    n_max = len_kmer_graph as u32;
                    println!("     . n_homopolymer was reduced to fit the k-mer length");
                }
                // test presence homopolymer
                is_homopolymer = check_homopolymer(n_max, first_kmer.clone(), l_middle_bases.clone());
            }
                            
            if is_homopolymer {
                nb_homopolymer += 1;
            } else {
                // Output seq_groups
                for l in &ll_sorted_seq {
                    let str_samples = l.1.join("|");
                    out_1.write_all(format!(">{}_{}\n{}\n", position, str_samples, l.0).as_bytes())
                    .expect("Failed to write to file");
                }
                
                // update binary alignment ('-' = missing data; '?' = both states (= 'unknown'))
                let mut sample_done: HashSet<String> = HashSet::new();
                let mut state = 0;
                for l in ll_sorted_seq.iter().take(2) {
                    // collect sample names for this sequence and update its vector in binary_seq
                    for sample_id in &l.1 {
                        let full_name = d_samples.get(sample_id).unwrap();
                        sample_done.insert(full_name.clone());
                        let sample_vec = binary_seq.get_mut(full_name).unwrap();
                    
                        if sample_vec.len() <= position {
                            sample_vec.resize(position + 1, state.to_string());
                        } else {
                            sample_vec[position] = "?".to_string();
                        }
                    }
                    // Update state
                    state += 1;
                }
            
                // update binary alignment with missing samples
                for sample in &all_samples {
                    if !sample_done.contains(sample) {
                        binary_seq.entry(sample.to_string()).and_modify(|vec| vec.push("-".to_string()));
                    }
                }                
            
                // get type of variant
                let mut type_variant = "complex";
                for xx in &l_middle_bases {
                    if xx == "." {type_variant = "_indel_";}
                }
               
                // get lists of samples
                let mut l_samples = Vec::new();
                for (_, vect) in &ll_sorted_seq {
                
                    // sort list of samples by equivalent integers
                    let mut vect2 = vect.clone();
                    vect2.sort_by(|a, b| {
                        let a_int = a.parse::<i32>().unwrap_or(i32::MAX);
                        let b_int = b.parse::<i32>().unwrap_or(i32::MAX);
                        a_int.cmp(&b_int)
                    });
                    // save it
                    l_samples.push(vect2.join(","));
                }
            
                // round ratio_missing to the 2nd decimal
                let rounded_missing = (ratio_missing * 100.0).round() / 100.0;

                // save variants information
                out_2.write_all(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    position, l_middle_bases.len(), type_variant, unclear_insert, first_kmer,
                    l_middle_bases.join(" / "), last_kmer, rounded_missing, l_samples.join(" / ")).as_bytes())
                    .expect("Failed to write to file");
            
                // Update position
                position += 1;
            }            
                   
        }
    }
    
    // Close output files
    out_1.flush().expect("Failed to flush file");
    out_2.flush().expect("Failed to flush file");

    // Output alignment
    let mut out_3 = File::create(format!("{}_binary_ali.fas", output_name)).expect("Failed to create file");
    for (sample, l_seq) in &binary_seq {
        out_3.write_all(format!(">{}\n{}\n", sample, l_seq.join("")).as_bytes())
            .expect("Failed to write to file");
    }
    out_3.flush().expect("Failed to flush file");
    
    // final message
    println!("     . {} removed because of missing data", nb_missing);
    if n_homopolymer.is_some() {
        println!("     . {} removed because in homopolymers", nb_homopolymer);
    }   
    println!("     . {} FINAL variant groups", position);
    println!("done.");
}




fn sort_by_nb_samples(mut map_seq: HashMap<String, Vec<String>>) -> Vec<(String, Vec<String>)> {
    let mut vec: Vec<_> = map_seq.drain().collect();
    vec.sort_by(|(_, v1), (_, v2)| v2.len().cmp(&v1.len())); // sort by decreasing order
    vec
}


fn calculate_ratio_missing(nb_total: usize, ll_seq: Vec<(String, Vec<String>)>) -> f32 {
    // initialise vector of 0
    let mut ref_samples: Vec<u32> = vec![0; nb_total];
    
    // iterate through ll_seq to count present samples
    for (_, samples) in ll_seq {
        for sample in samples {
            let idx = sample.parse::<usize>().unwrap();
            ref_samples[idx] += 1;
        }
    }

    // Count 'missing' samples (i.e., != 1)
    let missing_samples = ref_samples.iter().filter(|&&count| count != 1).count() as f32;

    // Calculate ratio of missing samples
    let ratio_missing = missing_samples / nb_total as f32;
    ratio_missing
}  


fn collect_middle_bases(ll_seq: Vec<(String, Vec<String>)>, len_k_graph: usize) -> (Vec<String>, String) {
    // collect all sequences without the first kmer
    let l_reduced_seq: Vec<&str> = ll_seq.iter().map(|l| &l.0[len_k_graph..]).collect();
    
    // get start position of last k-mer (i.e. find last position for which sequences differ (from the end))
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
    
    // extract last kmer (using first sequence)
    let pos_end = l_reduced_seq[0].len() - n_nucl;
    let mut last_kmer = l_reduced_seq[0][pos_end..].to_string();
    
    // the length of last k-mer might be in some very rare cases longer than expected (only observed in variants with lot of missing samples) -> truncate it
    if last_kmer.len() > len_k_graph {
        last_kmer = last_kmer[..len_k_graph].to_string();
    }

    // extract 'middle-bases' (remove last kmer from reduced sequences -> only middle base left)
    let mut l_middles: Vec<String> = Vec::new();
    for seq in &l_reduced_seq {
        let pos2_end = seq.len() - n_nucl;
        let mut middle_bases = &seq[..pos2_end];
        if middle_bases == "" {middle_bases = ".";}
        l_middles.push(middle_bases.to_string());
    }
    
    (l_middles, last_kmer)
}



fn test_multiple_positions(l_middles: Vec<String>, last_kmer: String) -> bool {
    
    for middle in l_middles {
        
        let full_string = middle.to_string() + &last_kmer;
    
        let mut occurrences_count = 0;
        let mut start = 0;

        while let Some(index) = full_string[start..].find(&last_kmer) {
            let full_index = start + index;
            occurrences_count += 1;
            if occurrences_count > 1 {
                return true;
            }
            start = full_index + 1;
        }
    }
    
    false // only one occurrence found
}



fn check_homopolymer(limit_n: u32, first_kmer: String, l_middle_bases: Vec<String>) -> bool {
    let mut is_homopolymer = false;
    
    // extract last limit_n nucleotides of first k-mer
    let last_n_chars = &first_kmer[first_kmer.len() - limit_n as usize..];
    
    // check if these last limit_n nucleotides are the same (ie, homopolymer)
    let char_set: HashSet<_> = last_n_chars.chars().collect();
    if char_set.len() == 1 {
        // check if the insert correspond to the homopolymer repeat
        'loop_middle: for middle in l_middle_bases {
            let middle_char_set: HashSet<_> = middle.chars().collect();
            if middle_char_set == char_set {
                is_homopolymer = true;
                break 'loop_middle;
            }
        }
    }
    is_homopolymer
}
