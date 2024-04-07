use std::collections::HashMap;
// use csv::Writer;
use std::env;
use std::fs::{self, read_dir, File};
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
fn main() {
    let test_top = PathBuf::from("/home/david/Documents/Academic/UVM/Harris_Lab/Mt_sequences");
    let genotype = env::args()
        .nth(1)
        .expect("please enter a valid genome name");
    // let top_dir = env::current_dir().expect("bad top dir");
    let genomes = create_full_path(test_top.clone(), String::from("genomes"));
    let annotation = create_full_path(test_top.clone(), String::from("annotations"));
    let dir_geno = get_name(genotype.clone(), genomes.clone());
    let dir_anno = get_name(genotype.clone(), annotation.clone());
    let full_geno = create_full_path(genomes.clone(), dir_geno[1].clone());
    let full_anno = create_full_path(genomes.clone(), dir_anno[1].clone());
    println!("{:?}", full_geno);
    println!("{:?}", full_anno);
}

fn create_full_path(tdir: PathBuf, dir: String) -> PathBuf {
    let mut s_dir = tdir.clone().into_os_string().into_string().expect("Nope");
    s_dir.push_str("/");
    s_dir.push_str(&dir);
    PathBuf::from(s_dir)
}

fn get_name(pat: String, search_dir: PathBuf) -> Vec<String> {
    let mut patmatch = String::new();
    let mut mdir = String::new();
    if let Ok(files) = read_dir(search_dir) {
        for file in files {
            let sdir = file.unwrap().path().into_os_string().into_string().unwrap();
            let sp = sdir.split('/').last().unwrap().split('.').next().unwrap();
            if pat == sp {
                patmatch.push_str(&sp);
                let msp = sdir.split('/').last().unwrap();
                mdir.push_str(msp);
            }
        }
    } else {
        println!("Path not found");
    }
    let mut rvec: Vec<String> = Vec::new();
    rvec.push(patmatch);
    rvec.push(mdir);
    rvec
}

fn read_fasta<P>(filename: P) -> HashMap<String, String>
where
    P: AsRef<Path>,
{
    let file = File::open(filename).expect("Could not open file");
    let buf = BufReader::new(file);
    let mut fasta = HashMap::new();
    let mut curid = String::new();
    let mut curseq = String::new();

    for line in buf.lines() {
        let line = line.expect("Could not read line");
        if line.starts_with('>') {
            if !curid.is_empty() {
                fasta.insert(curid.clone(), curseq.clone());
                curseq.clear();
            }
            curid = line[..].trim().to_string();
        } else {
            curseq.push_str(line.trim());
        }
    }
    fasta
}
