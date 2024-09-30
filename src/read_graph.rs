use hashbrown::{HashMap, HashSet};

use crate::utils::{rev_compl, decode_kmer, get_last_nucleotide};


pub fn build_sequences(len_kmer: usize, all_kmers: &HashMap<u128, HashMap<u128, u32>>, start_kmers: &HashSet<u128>, end_kmers: &HashSet<u128>, index_map: &HashMap<u32, String>) -> HashMap<String, HashMap<String, Vec<String>>> {
    println!(" # explore graph");

    let len_kmer_graph = len_kmer -1;

    let mut built_sequences: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
    let mut seq_done: HashSet<String> = HashSet::new();
    
    // build sequences from entry nodes
    for kmer in start_kmers {
        
        // initialise container for sequences built
        let mut seq_found: HashMap<String, HashMap<String, Vec<String>>> = HashMap::new();
        
        // get k-mers to test
        let mut kmers_to_test = Vec::<&u128>::new();
        for (next_kmer, _) in all_kmers.get(kmer).unwrap().iter() {
            kmers_to_test.push(next_kmer);
        }
        kmers_to_test.sort();
        
        // consider all nucleotides and build sequence for each of them
        for next_kmer in kmers_to_test {
            // create new sequence by adding last nucleotide of the next kmer to the main kmer
            let mut sequence = decode_kmer(kmer.clone(), len_kmer_graph);            
            let next_nucl = get_last_nucleotide(next_kmer.clone());
            sequence += &next_nucl.to_string();

            // create sample variables to be used during the walk            
            let tmp_index_samples = all_kmers[kmer][next_kmer];
            let sample_set_ref: HashSet<&str> = index_map.get(&tmp_index_samples).unwrap().split('|').collect();
            
            let mut d_nb_samples: HashMap<&str, i32> = HashMap::new();
            for sample in sample_set_ref.iter() {
                *d_nb_samples.entry(sample).or_insert(0) += 1;
            }
            
            // create set of visited nodes
            let mut visited: HashSet<u128> = HashSet::new();
            visited.insert(kmer.clone());

            let mut previous_kmer: u128 = next_kmer.clone();
            let mut walking_along_path = true;
            
            let mut tmp_l_index = Vec::new();
            tmp_l_index.push(tmp_index_samples);
            
            while walking_along_path {   
                
                // compare samples of all next kmers with s_ref_samples
                let mut good_next = Vec::new();
                
                if let Some(next_kmer_data) = all_kmers.get(&previous_kmer) {
                    // get among next kmer(s) those that were not already visited
                    let mut kmer2_to_test = Vec::new();
                    for kmer2 in next_kmer_data.keys() {
                        // we ignore already visited kmers to avoid endless loops
                        if !visited.contains(kmer2) {
                            kmer2_to_test.push(kmer2.clone());
                        }
                    }
                    
                    // if only 1 kmer -> save it; if not -> compute Jaccard similarity
                    if kmer2_to_test.len() == 1 {
                        good_next.extend_from_slice(&kmer2_to_test);
                    } else {
                        // /*  uncomment here and below to remove the path selection based on Jaccard similarity
                        let mut max_sim = 0.0;
                        for kmer2 in kmer2_to_test {
                            let index_samples = next_kmer_data.get(&kmer2).unwrap();
                            let s_samples: HashSet<&str> = index_map.get(index_samples).unwrap().split('|').collect();
                            let sim = jaccard_similarity(&sample_set_ref, &s_samples);
                            if sim > max_sim {
                                good_next.clear();
                                good_next.push(kmer2.clone());
                                max_sim = sim;
                            } else if sim == max_sim {
                                good_next.push(kmer2.clone());
                            }
                        }
                        // */
                    }
                }
                
                // case only 1 next kmer
                if good_next.len() == 1 {   

                    // update sequence
                    let next_nucl = get_last_nucleotide(good_next[0]);
                    sequence += &next_nucl.to_string();

                    // update set visited
                    visited.insert(good_next[0]);
                    
                    tmp_l_index.push(all_kmers[&previous_kmer][&good_next[0]]);
                                        
                    // update previous_kmer
                    previous_kmer = good_next[0].clone();
                    
                    // save sequence if the kmer was an end kmer
                    if end_kmers.contains(&good_next[0]) {

                        // get consensus limit to rebuild s_ref_samples
                        let limit_consensus = (visited.len() as f32 * 0.5) as i32;
                        
                        // build s_ref_samples with majority rule consensus
                        for (tmp_index, tmp_count) in count_occurrences(&tmp_l_index) {                        
                            let tmp_samples: HashSet<&str> = index_map.get(&tmp_index).unwrap().split('|').collect();
                            for sample in tmp_samples.iter() {
                                d_nb_samples.entry(sample).and_modify(|count| *count += tmp_count).or_insert(tmp_count);
                            }
                        }
                        tmp_l_index.clear();
                        
                        let mut s_ref_samples: HashSet<&str> = HashSet::new();
                        for (sample, nb) in &d_nb_samples {
                            if nb >= &limit_consensus {
                                s_ref_samples.insert(&sample);
                            }
                        }
                        
                        // save sequence                                    
                        let combined_ends = format!("{}@{}", kmer, good_next[0]);
                        seq_found
                            .entry(combined_ends.clone())
                            .or_insert_with(HashMap::new)
                            .insert(sequence.clone(), s_ref_samples.iter().cloned().map(|s| s.to_string()).collect());
                        
                        // stop looking if another branch already ended on this end kmer
                        if seq_found[&combined_ends].len() > 1 {
                            walking_along_path = false;
                        }
                    }
                // case no next kmer or several good next kmers
                } else {
                    walking_along_path = false;
                }
            }            
            
        }

        // check built sequences
        for (combined_ends, d_seq) in seq_found.iter() {
            if d_seq.len() > 1 {
                // check that the pair of sequences isn't a snp (i.e., same length AND levenshtein distance = 1)
                let mut not_a_snp = false;
                let l_seq: Vec<&String> = d_seq.keys().collect();
                for (i, seq1) in l_seq.iter().enumerate() {
                    for seq2 in l_seq.iter().skip(i + 1) {
                        if seq1.len() == seq2.len() {
                            if levenshtein_distance(seq1, seq2) != 1 {
                                not_a_snp = true;
                                break;
                            }
                        } else {
                            not_a_snp = true;
                            break;
                        }
                    }
                }

                if not_a_snp {
                    // check if reverse complement already done
                    let mut new_d: HashMap<String, Vec<String>> = HashMap::new();
                    for (seq, l_samples) in d_seq.iter() {
                        let rc_seq = rev_compl(seq);
                        if !seq_done.contains(&rc_seq) {
                            new_d.insert(seq.clone(), l_samples.clone());
                            seq_done.insert(seq.clone());
                        }
                    }

                    // save sequences if still more than one sequence
                    if new_d.len() > 1 {
                        for (seq, l_samples) in new_d.iter() {
                            seq_done.insert(seq.clone());
                            built_sequences
                                .entry(combined_ends.clone())
                                .or_insert_with(HashMap::new)
                                .insert(seq.clone(), l_samples.clone());
                        }
                    }
                }
            }
        }
    }
    
    println!("     . {} variant groups", built_sequences.len());
    built_sequences
}


fn jaccard_similarity(set1: &HashSet<&str>, set2: &HashSet<&str>) -> f32 {
    // returns the Jaccard similarity value between 2 sets
    let intersection_size = set1.intersection(set2).count() as f32;
    let union_size = (set1.len() as f32 + set2.len() as f32 - intersection_size) as f32;
    intersection_size / union_size
}


fn levenshtein_distance(str1: &str, str2: &str) -> usize {
    // returns the levenshtein distance between 2 strings
    let l1: Vec<char> = str1.chars().collect();
    let l2: Vec<char> = str2.chars().collect();
    let mut dist = 0;

    for (i, nucl) in l1.iter().enumerate() {
        if nucl != &l2[i] {
            dist += 1;
        }
    }
    dist
}


fn count_occurrences(vec: &Vec<u32>) -> HashMap<u32, i32> {
    let mut count_map = HashMap::new();
    for num in vec {
        *count_map.entry(*num).or_insert(0) += 1;
    }

    count_map
}
