use hashbrown::HashMap;


pub fn filter_sequences(sequences: HashMap<String, HashMap<String, Vec<String>>>, len_kmer: usize, all_samples: Vec<String>, n_homopolymer: u32, max_missing: f32) -> HashMap<String, Vec<(String, Vec<String>)>> {
    println!(" . filter sequences");
    
    // initialise output hashmap
    let mut filtered_sequences: HashMap<String, Vec<(String, Vec<String>)>> = HashMap::new();
    
    for (combined_ends, d_seq) in sequences.iter() {
        let mut homopolymer_found = false;
        
        for seq in d_seq.keys() {
            
            // extract last n nucleotides of left kmer (or the entire left kmer if n_homopolymer too big)
            // in theory we could just do it for the 1st sequence since they all have the same left kmer
            if n_homopolymer > (len_kmer as u32 -1) {
                let test_segment = &seq[0 .. len_kmer as usize -1];
                if test_segment.chars().all(|c| c == test_segment.chars().next().unwrap()) {
                    homopolymer_found = true;
                    //println!("{}\n{}   homo_long\n", &seq, &test_segment);
                }
            
            } else {
                let test_segment = &seq[len_kmer as usize - 1 - n_homopolymer as usize .. len_kmer as usize -1];
                if test_segment.chars().all(|c| c == test_segment.chars().next().unwrap()) {
                    homopolymer_found = true;
                    //println!("{}\n                        {}   homo\n", &seq, &test_segment);            
                }
            } 
        }
        
        if !homopolymer_found {
            let mut tmp_samples: HashMap<String, u32> = all_samples.iter().enumerate().map(|(i, _)| (i.to_string(), 0)).collect();
            for l_samples in d_seq.values() {
                for sample in l_samples {
                    if let Some(count) = tmp_samples.get_mut(sample) {
                        *count += 1;
                    }
                }
            }
            
            // unknown samples correspond to those absent in all sequences or present in more than 1 sequence
            let nb_unknown = tmp_samples.values().filter(|&&count| count != 1).count() as f32;
            let missing_ratio = nb_unknown / all_samples.len() as f32;

            if missing_ratio <= max_missing {
                let mut sorted_d_seq: Vec<(&String, &Vec<String>)> = d_seq.iter().collect();
                sorted_d_seq.sort_by(|(_, l_samples1), (_, l_samples2)| l_samples2.len().cmp(&l_samples1.len()));

                let filtered_entries: Vec<(String, Vec<String>)> = sorted_d_seq
                    .into_iter()
                    .map(|(seq, l_samples)| {
                        let mut l_samples_clone = l_samples.clone();
                        l_samples_clone.sort();
                        (seq.clone(), l_samples_clone)
                    })
                    .collect();

                filtered_sequences.insert(combined_ends.clone(), filtered_entries);
            }
        }
    }

    filtered_sequences
}
