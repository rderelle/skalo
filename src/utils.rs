use hashbrown::HashMap;



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
