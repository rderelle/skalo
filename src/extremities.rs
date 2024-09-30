use hashbrown::{HashMap, HashSet};

use crate::utils::{rev_compl, encode_kmer, decode_kmer};



pub fn identify_good_kmers(len_kmer: usize, all_kmers: &HashMap<u128, HashMap<u128, u32>>, index_map: &HashMap<u32, String>) -> (HashSet<u128>, HashSet<u128>) {
    println!(" # identify bubble extremities");
    
    let len_kmer_graph = len_kmer -1;
    
    let mut start_kmers: HashSet<u128> = HashSet::new();
    let mut end_kmers: HashSet<u128> = HashSet::new();
    
    for (kmer, next_kmers_map) in all_kmers.iter() {
        if next_kmers_map.len() > 1 {
            let all_next_kmer: Vec<&u128> = next_kmers_map.keys().collect();

            'i_loop: for (i, &next_kmer) in all_next_kmer.iter().enumerate() {
                for &next_kmer2 in all_next_kmer.iter().skip(i + 1) {
                    let samples_1 = &index_map[&next_kmers_map[next_kmer]];
                    let samples_2 = &index_map[&next_kmers_map[next_kmer2]];

                    if compare_samples(samples_1, samples_2) {
                        start_kmers.insert(kmer.clone());
                                               
                        let dna = decode_kmer(kmer.clone(), len_kmer_graph);
                        let rc = rev_compl(&dna);

                        //uncomment to print network
                        //println!("{}	{}	red", &dna, &dna);
                        
                        end_kmers.insert(encode_kmer(&rc));

                        //uncomment to print network
                        //println!("{}	{}	red", &rc, &rc);
                        break 'i_loop;
                    }
                }
            }
        }
    }
    
    // exit program if no extremity found (e.g. cases of weeded skf files)
    if start_kmers.is_empty() {
        eprintln!("\n      Error: there is no entry node in this graph, hence no variant.\n");
        std::process::exit(1);
    }
    
    println!("     . {} entry nodes", start_kmers.len());
    (start_kmers, end_kmers)
}


fn compare_samples(str_1: &str, str_2: &str) -> bool {
    let set1: HashSet<&str> = str_1.split('|').collect();
    let set2: HashSet<&str> = str_2.split('|').collect();
    let diff_1 = set1.difference(&set2).count();
    let diff_2 = set2.difference(&set1).count();

    diff_1 != 0 && diff_2 != 0
}
