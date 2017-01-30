extern crate rand;

use rand::Rng;
use std::collections::HashMap;
use std::env;
use std::io;

/* ======= CONSTANTS ======= */

const _ASCII_LEN: usize = 255;



/* ======= UTILITY FUNCTIONS ======= */

fn read_strings(a: &mut String, b: &mut String) {
    io::stdin().read_line(a);
    io::stdin().read_line(b);
    
    *a = a.trim().to_string();
    *b = b.trim().to_string();
}

fn extract_sigma(a: &String) -> String {
    let mut array: [bool; _ASCII_LEN] = [false; _ASCII_LEN];
    for c in a.chars() {
        array[c as usize] = true;
    }
    
    let mut sigma = String::new();
    for x in 0.._ASCII_LEN {
        if array[x] == true {
           sigma.push((x as u8) as char);
        }
    }
    
    sigma
}

fn define_mapping(sigma: &String, map: &mut HashMap<char, usize>) {
    let mut i: usize = 0;
    for c in sigma.chars() {
        map.insert(c, i);
        i += 1;
    }
}

fn map_vec(map: &HashMap<char, usize>, v: &mut Vec<usize>, s: &String) {
    let mut i: usize = 0;
    for c in s.chars() {
        v[i] = *map.get(&c).unwrap();
        i += 1;
    }
}

fn is_valid(s: &Vec<usize>, p: usize) -> bool{
	let mut current_min: u32 = u32::max_value()-1;
	let mut last_min: u32 = u32::max_value();
    
    for i in 0..s.len() {
        if current_min > (s[i] as u32) {
            current_min = s[i] as u32;
        }
        
        if i >= p {
            if last_min > current_min {
                return false
            }
        }
        
        if (i % p) == (p-1) {
            last_min = current_min;
            current_min = u32::max_value();
        }
    }
	true
}

fn shuffle(a: &mut Vec<usize>) {
    let mut rng = rand::thread_rng();
    
    let mut j: usize = 0;
    let mut t: usize = 0;
    
    for i in (0..a.len()).rev() {
        j = rng.gen::<usize>() % (i+1);
        t = a[j]; a[j] = a[i]; a[i] = t;
    }
}


/* ======= MATCHING SCHEMA (struct & functions) ======= */

struct MatchingSchema {
    // Can be initialized such as
    // let mut m: MatchingSchema = MatchingSchema { n: X, m: Y, ms: vec![true; X * Y] };
    n: usize,
    m: usize,
    ms: Vec<bool>
}

impl MatchingSchema {
    fn init_general(&mut self, p1: &usize, p2: &usize) {
        let p1_u: usize = *p1 as usize;
        let p2_u: usize = *p2 as usize;
        
        for i in 0..self.n {
            for j in 0..self.m {
                self.set(i, j, if *p1 != 0 && *p2 != 0 && (i/p2_u == j/p1_u) {false} else {true})
            }
        }
    }
    
    fn set(&mut self, i: usize, j: usize, value: bool) {
        self.ms[(i * self.m)+ j] = value;
    }

    fn get(&self, i: usize, j: usize) -> bool {
        self.ms[(i * self.m) + j]
    }
}


/* ======= SOLVERS et co. ======= */

fn edit_distance_enhanced(a: &Vec<usize>, b: &Vec<usize>, sig1: &Vec<usize>, sig2: &Vec<usize>,  matrix: &mut Vec<Vec<u32>>, m: &MatchingSchema) -> u32 {
    let len_a = a.len();
    let len_b = b.len();
    
    // We build the new index for the symbols
    let mut sig1_index: Vec<usize> = vec![0; sig1.len() as usize];
    let mut sig2_index: Vec<usize> = vec![0; sig2.len() as usize];
    for i in 0..sig1.len() { sig1_index[sig1[i]] = i; }
    for i in 0..sig2.len() { sig2_index[sig2[i]] = i; }
    
    for i in 0..len_a+1 { matrix[i][0] = i as u32; }
    for j in 0..len_b+1 { matrix[0][j] = j as u32; }
    
    for i in 1..len_a+1 {
        for j in 1..len_b+1 {
            let w: u32 = if m.get(sig1_index[a[i-1]], sig2_index[b[j-1]]) {1} else {0};
            
            let ops = [
                matrix[i-1][j] + 1,
                matrix[i][j-1] + 1,
                matrix[i-1][j-1] + w
            ];
            
            matrix[i][j] = *ops.iter().min().unwrap();
        }
    }
    
    matrix[len_a][len_b]
}

fn hill_climbing(a: &Vec<usize>, b: &Vec<usize>, sig1: &Vec<usize>, sig2: &Vec<usize>, p: &usize, m: &MatchingSchema) -> u32 {
    // matrix (external)
    let row: Vec<u32> = vec![0; b.len()+1];
    let mut matrix: Vec<Vec<u32>> = vec![row; a.len()+1];
    
    // compute first edit
    let mut d: u32 = edit_distance_enhanced(&a, &b, &sig1, &sig2, &mut matrix, &m);
    let mut min_dist: u32 = d;
    let mut min_min_dist: u32 = d;
    
    // vectors for the permutations
    let mut isig1_o: Vec<usize> = vec![0; sig1.len()]; for i in 0..sig1.len() { isig1_o[i] = i; }
    let mut isig2_o: Vec<usize> = vec![0; sig2.len()]; for i in 0..sig2.len() { isig2_o[i] = i; }
    let mut isig1_t: Vec<usize> = vec![0; sig1.len()]; for i in 0..sig1.len() { isig1_t[i] = i; }
    let mut isig2_t: Vec<usize> = vec![0; sig2.len()]; for i in 0..sig2.len() { isig2_t[i] = i; }
    // vectors for the fixpoint
    let mut isig1_min: Vec<usize> = vec![0; sig1.len()]; for i in 0..sig1.len() { isig1_min[i] = i; }
    let mut isig2_min: Vec<usize> = vec![0; sig2.len()]; for i in 0..sig2.len() { isig2_min[i] = i; }
    let mut isig1_min_min: Vec<usize> = vec![0; sig1.len()]; for i in 0..sig1.len() { isig1_min_min[i] = i; }
    let mut isig2_min_min: Vec<usize> = vec![0; sig2.len()]; for i in 0..sig2.len() { isig2_min_min[i] = i; }

    let mut temp: usize = 0;
    let attempts: usize = 1;
    let shuffle_tries: usize = 2;
    let mut k_shuffle: usize = 0;
    let mut tries: usize = 0;
    let mut improved: bool = true;
    
    while improved {
        improved = false;
        
        for ip in 0..sig1.len() {
            for jp in ip..sig1.len() {
                
                isig1_o = isig1_t.clone(); 
                
                temp = isig1_o[ip];
                isig1_o[ip] = isig1_o[jp];
                isig1_o[jp] = temp;
                
                if is_valid(&isig1_o, *p) {
                    
                    for ipp in 0..sig2.len() {
                        for jpp in ipp..sig2.len() {
                        
                            isig2_o = isig2_t.clone();
                        
                            temp = isig2_o[ipp];
                            isig2_o[ipp] = isig2_o[jpp];
                            isig2_o[jpp] = temp;
                        
                            d = edit_distance_enhanced(&a, &b, &isig1_o, &isig2_o, &mut matrix, &m);
                            if d < min_dist {
                                improved = true;
                                min_dist = d;
                            
                                isig1_min = isig1_o.clone();
                                isig2_min = isig2_o.clone();
                            }
                        }                    
                        isig2_o = isig2_t.clone();
                    }
                    
                }
            }
        }
        
        if improved {
            isig1_o = isig1_min.clone();
            isig2_o = isig2_min.clone();
            
            isig1_t = isig1_o.clone();
            isig2_t = isig2_o.clone();
        } else {
            if min_dist < min_min_dist {
                min_min_dist = min_dist;
                
                isig1_min_min = isig1_min.clone();
                isig2_min_min = isig2_min.clone();
                
                improved = true;
                tries = 0;
            }
            
            if tries < attempts {
                improved = true;
                tries += 1;
                
                k_shuffle = 0;
                while k_shuffle < shuffle_tries && !is_valid(&isig1_t, *p) {
                    shuffle(&mut isig1_t);
                }
                if k_shuffle == shuffle_tries {
                    isig1_t = isig1_o.clone();
                }
                shuffle(&mut isig2_t);
                
                isig1_o = isig1_t.clone();
                isig2_o = isig2_t.clone();
                
                min_dist = edit_distance_enhanced(&a, &b, &isig1_o, &isig2_o, &mut matrix, &m);
            }
        }
    }
    
    min_min_dist
}

fn main() {
    let mut s1 = String::new();
    let mut s2 = String::new();
    let mut p1: usize = 1;
    let mut p2: usize = 1;
    
    // I create a vector of string in order to parse arguments
    let args: Vec<String> = env::args().collect();
    // Update p1 and p2 if used
    if args.len() > 2 {
        p1 = args[1].parse().unwrap();
        p2 = args[2].parse().unwrap();
    }
    
    // I read the strings
    read_strings(&mut s1, &mut s2);
    
    // I define sigmas
    let sigma1 = extract_sigma(&s1);
    let sigma2 = extract_sigma(&s2);
    
    /*println!("π1: {}, π2: {}.", p1, p2);
    println!("∑1: {} (len = {}), ∑2: {} (len = {}).", sigma1, sigma1.len(), sigma2, sigma2.len());*/
    
    // Mappings
    let mut map_sigma1: HashMap<char, usize> = HashMap::new();
    let mut map_sigma2: HashMap<char, usize> = HashMap::new();
    define_mapping(&sigma1, &mut map_sigma1);
    define_mapping(&sigma2, &mut map_sigma2);
    
    // Vectors
    let mut isigma1: Vec<usize> = vec![0; sigma1.len()];
    let mut isigma2: Vec<usize> = vec![0; sigma2.len()];
    let mut is1: Vec<usize> = vec![0; s1.len()];
    let mut is2: Vec<usize> = vec![0; s2.len()];
    // Populate vectors
    for i in 0..sigma1.len() { isigma1[i] = i; }
    for i in 0..sigma2.len() { isigma2[i] = i; }
    map_vec(&map_sigma1, &mut is1, &s1);
    map_vec(&map_sigma2, &mut is2, &s2);
    
    //println!("{:?}, {:?}", is1, is2);
    //println!("{:?}, {:?}", isigma1, isigma2);
    
    // Matching schema
    let mut m: MatchingSchema = MatchingSchema { n: sigma1.len(), m: sigma2.len(), 
        ms: vec![true; sigma1.len() * sigma2.len()] };
    m.init_general(&p1, &p2);
    
    /*println!("M:");
    for i in 0..m.n {
        for j in 0..m.m {
            print!("{} ", if m.get(i, j) {1} else {0})
        }
        println!();
    }*/
    
    //println!("{}", edit_distance(&s1, &s2));
    //println!("{}", edit_distance_enhanced(&is1, &is2, &isigma1, &isigma2, &m));
    println!("{}", hill_climbing(&is1, &is2, &isigma1, &isigma2, &p1, &m));
}
