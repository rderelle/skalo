use hashbrown::{HashMap, HashSet};



pub fn rev_compl(seq: &str) -> String {
    let reverse_code: HashMap<char, char> = [
        ('A', 'T'),
        ('C', 'G'),
        ('G', 'C'),
        ('T', 'A'),
        ('N', 'N'),
    ]
    .iter()
    .cloned()
    .collect();

    let reversed_seq: String = seq.chars().rev().collect();
    let complemented_seq: String = reversed_seq.chars().map(|n| reverse_code[&n]).collect();

    complemented_seq
}


pub fn rev_compl_u128(kmer: u128, k: usize) -> u128 {
    // mask for the last 2 bits (representing one nucleotide)
    let mask = 0b11u128;

    // initialize reverse complement result
    let mut rc_u128 = 0u128;

    // for each nucleotide in the k-mer
    for i in 0..k {
        // get the last 2 bits (current nucleotide)
        let nucleotide = (kmer >> (2 * i)) & mask;

        // complement the nucleotide:
        let complement = match nucleotide {
            0b00 => 0b11, // A -> T
            0b01 => 0b10, // C -> G
            0b10 => 0b01, // G -> C
            0b11 => 0b00, // T -> A
            _ => unreachable!(),
        };

        // shift the complemented nucleotide to its reverse position
        rc_u128 |= complement << (2 * (k - i - 1));
    }
    rc_u128
}



pub fn encode_kmer(kmer: &str) -> u128 {
    let nucleotide_to_bits: [u8; 4] = [
        0b00, // A
        0b01, // C
        0b10, // G
        0b11, // T
    ];

    let mut result: u128 = 0;

    for nucleotide in kmer.chars() {
        let index = match nucleotide {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => panic!("Invalid nucleotide"),
        };
        result = (result << 2) | (nucleotide_to_bits[index] as u128);
    }
    result
}



pub fn decode_kmer(encoded: u128, k: usize) -> String {
    let bits_to_nucleotide: [char; 4] = ['A', 'C', 'G', 'T'];
    let mut kmer = String::with_capacity(k);

    let mask = (1u128 << (2 * k)) - 1;
    let mut value = encoded & mask;

    for _ in 0..k {
        let index = (value & 0b11) as usize;
        let nucleotide = bits_to_nucleotide[index];
        kmer.insert(0, nucleotide);
        value >>= 2;
    }
    kmer
}


pub fn get_last_nucleotide(encoded_kmer: u128) -> char {
    // mask the last 2 bits to get the encoded nucleotide
    let last_bits = (encoded_kmer & 0b11) as u8;
    // decode the nucleotide based on the 2-bit pattern
    match last_bits {
        0b00 => 'A',
        0b01 => 'C',
        0b10 => 'G',
        0b11 => 'T',
        _ => unreachable!(),
    }
}


/// structure to save variant information

#[derive(Clone)]
pub struct VariantInfo<'a> {
    //pub entry_kmer: u128,
    //pub exit_kmer: u128,
    pub sequence: String,
    pub is_snp: bool,
    pub visited_nodes: HashSet<u128>,
    pub count_samples: HashMap<&'a str, i32>,
    pub maj_samples: HashSet<&'a str>,
}

impl<'a> VariantInfo<'a> {
    pub fn new(
        //entry_kmer: u128,
        //exit_kmer: u128,
        sequence: String,
        is_snp: bool,
        visited_nodes: HashSet<u128>,
        count_samples: HashMap<&'a str, i32>,
        maj_samples: HashSet<&'a str>,
    ) -> Self {
        VariantInfo {
            //entry_kmer,
            //exit_kmer,
            sequence,
            is_snp,
            visited_nodes,
            count_samples,
            maj_samples,
        }
    }
    
    /*
    pub fn entry_kmer(&self) -> &u128 {
        &self.entry_kmer
    }

    pub fn exit_kmer(&self) -> &u128 {
        &self.exit_kmer
    }
    
    pub fn sequence(&self) -> &String {
        &self.sequence
    }

    pub fn is_snp(&self) -> bool {
        self.is_snp
    }

    pub fn visited_nodes(&self) -> &HashSet<u128> {
        &self.visited_nodes
    }

    pub fn count_samples(&self) -> &HashMap<&'a str, i32> {
        &self.count_samples
    }

    pub fn maj_samples(&self) -> &HashSet<&'a str> {
        &self.maj_samples
    }
    */
}


