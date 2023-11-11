use hashbrown::{HashMap, HashSet};

use crate::utils::{rev_compl, encode_kmer, decode_kmer};



pub fn identify_good_kmers(len_kmer: usize, all_kmers: &HashMap<u128, HashMap<u128, u64>>, index_map: &HashMap<u64, String>) -> HashMap<u8, HashSet<u128>> {

    println!(" . identify good kmers");
    
    let len_kmer_graph = len_kmer -1;
    let mut good_kmers: HashMap<u8, HashSet<u128>> = HashMap::new();

    for (kmer, next_kmers_map) in all_kmers.iter() {
        if next_kmers_map.len() > 1 {
            let all_next_kmer: Vec<&u128> = next_kmers_map.keys().collect();

            for (i, &next_kmer) in all_next_kmer.iter().enumerate() {
                for &next_kmer2 in all_next_kmer.iter().skip(i + 1) {
                    let samples_1 = &index_map[&next_kmers_map[next_kmer]];
                    let samples_2 = &index_map[&next_kmers_map[next_kmer2]];

                    if compare_samples(samples_1, samples_2) {
                        good_kmers
                            .entry(1)
                            .or_insert_with(HashSet::new)
                            .insert(kmer.clone());
                        
                        let dna = decode_kmer(kmer.clone(), len_kmer_graph);
                        let rc = rev_compl(&dna);
                        
                        good_kmers
                            .entry(2)
                            .or_insert_with(HashSet::new)
                            .insert(encode_kmer(&rc));
                        break;
                    }
                }
            }
        }
    }
    
    good_kmers
}


fn compare_samples(str_1: &str, str_2: &str) -> bool {
    let set1: HashSet<&str> = str_1.split('|').collect();
    let set2: HashSet<&str> = str_2.split('|').collect();
    let diff_1 = set1.difference(&set2).count();
    let diff_2 = set2.difference(&set1).count();

    diff_1 != 0 && diff_2 != 0
}



