use csv::Writer;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::f32::NAN;
use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::Instant;
use std::{env, usize};

fn main() {
    env::set_var("RUST_BACKTRACE", "0");

    let start = Instant::now();
    let test_top1 = PathBuf::from("/home/david/Documents/Academic/UVM/Harris_Lab/Mt_sequences");
    // let test_top2 = PathBuf::from("/media/david/WorkDrive/Documents/UVM/HarrisLab/Mt_Data");
    let genotype = env::args()
        .nth(1)
        .expect("please enter a valid genome name (arg 1)");
    let size_string = env::args()
        .nth(2)
        .expect("please enter the size for the search window (arg 2)");
    let size = size_string.trim().parse::<i32>().unwrap();
    let pattern = env::args()
        .nth(3)
        .expect("please enter the pattern to search for (arg 3)");
    // let top_dir = env::current_dir().expect("bad top dir");
    let genomes = create_full_path(test_top1.clone(), String::from("genomes"));
    let annotation = create_full_path(test_top1.clone(), String::from("annotations"));
    let dir_geno = get_name(genotype.clone(), genomes.clone());
    let dir_anno = get_name(genotype.clone(), annotation.clone());
    let full_geno = create_full_path(genomes.clone(), dir_geno.clone());
    let full_anno = create_full_path(annotation.clone(), dir_anno.clone());
    assert!(Path::new(&full_anno).exists());
    assert!(Path::new(&full_geno).exists());

    let decogeno = read_fasta(full_geno);
    let decoanno = read_fasta(full_anno);

    let full_genome = chromosomes(decogeno);

    let akeys = decoanno.keys();
    for key in akeys {
        let info = get_info(parse_header(key.to_string()));
        let id = &info[0];
        let strand = get_strand(info.clone());
        println!("{:?}", id);
        let begend = get_begin_end(info.clone());
        let window = build_window(begend[0], begend[1], size);
        search_seq(full_genome.clone(), window, pattern.clone(), strand);
    }

    let duration = start.elapsed();
    println!("Total time to completion: {:?}", duration);
}

fn chromosomes(hash: HashMap<String, String>) -> String {
    let chrom = String::from(
"MtrunA17Chr1,MtrunA17Chr2,MtrunA17Chr3,MtrunA17Chr4,MtrunA17Chr5,MtrunA17Chr6,MtrunA17Chr7,MtrunA17Chr8"
    );
    let mut full_genome = String::new();
    let ckeys = hash.keys();

    for chm in chrom.split(",") {
        for key in ckeys.clone() {
            if key.to_string().contains(chm) {
                let val = hash.get(key).unwrap();
                full_genome.push_str(val);
            }
        }
    }

    full_genome
}

fn get_begin_end(info: Vec<String>) -> Vec<i32> {
    let info_split = info[0].split(" ");
    let mut begend: Vec<i32> = Vec::new();
    let mut begin = String::new();
    let mut end = String::new();
    let mut misc = String::new();

    for fo in info_split {
        let fon = fo.split("=").next().unwrap();
        let fol = fo.split("=").last().unwrap();
        match fon {
            "begin" => begin.push_str(&fol),
            "end" => end.push_str(&fol),
            _ => misc.push_str(&fol),
        }
    }

    let ibegin = begin.to_string().trim().parse::<i32>().unwrap();
    let iend = end.to_string().trim().parse::<i32>().unwrap();

    begend.push(ibegin);
    begend.push(iend);
    assert!(begend.len() == 2, "{:?} Did not pass", info[0]);
    assert!(begend[0] < begend[1]);
    begend
}

fn build_window(begin: i32, end: i32, size: i32) -> Vec<i32> {
    let left_bound: i32 = begin - size;
    let right_bound: i32 = end + size;
    let mut window: Vec<i32> = Vec::new();
    window.push(left_bound);
    window.push(right_bound);
    window
}

fn get_strand(info: Vec<String>) -> String {
    let info = &info[0];
    let info_sp = info.split(" ");
    let mut strand_type = String::new();
    for element in info_sp {
        if element.contains("strand") {
            strand_type.push_str(element);
        }
    }
    strand_type
}

fn minus_strand_invsersion(seq: String, pat: String) -> Vec<String> {
    let mut inversion: Vec<String> = Vec::new();
    let mut seq_inv = String::new();
    let pat_inv = pat.chars().rev().collect();

    for nuc in seq.chars() {
        match nuc {
            'a' => seq_inv.push_str("t"),
            'c' => seq_inv.push_str("g"),
            't' => seq_inv.push_str("a"),
            'g' => seq_inv.push_str("c"),
            'n' => seq_inv.push_str("n"),
            _ => println!("This is not a nucleotide, or n"),
        }
    }

    inversion.push(seq_inv);
    inversion.push(pat_inv);

    println!("Sequence and pattern inverted");

    inversion
}

fn search_seq(seq: String, window: Vec<i32>, pattern: String, strand: String) {
    let strand = strand.split("=").last().unwrap();
    if strand == "-" {
        let seq = seq.to_lowercase();
        let pattern = pattern.to_lowercase();
        let inversion = minus_strand_invsersion(seq, pattern);
        let seq = &inversion[0];
        let pattern = &inversion[1];
        let bseq = seq.as_bytes();
        let bpat = pattern.as_bytes();
        let left_bound = window[0] as usize;
        let right_bound = window[1] as usize;
        let search_area = &bseq[left_bound..right_bound];

        for i in 0..search_area.len() - bpat.len() {
            if bpat == &search_area[i..i + bpat.len()] {
                let location = left_bound + i;
                println!("Occurance at {}", location);
            }
        }
    } else {
        let seq = seq.to_lowercase();
        let pattern = pattern.to_lowercase();
        let bseq = seq.as_bytes();
        let bpat = pattern.as_bytes();
        let left_bound = window[0] as usize;
        let right_bound = window[1] as usize;
        let search_area = &bseq[left_bound..right_bound];

        for i in 0..search_area.len() - bpat.len() {
            if bpat == &search_area[i..i + bpat.len()] {
                let location = left_bound + i;
                println!("Occurance at {}", location);
            }
        }
    }
}

fn collapse_lines(hash: HashMap<String, String>) -> String {
    let seqs = hash.values();
    let mut full = String::new();
    for idv in seqs {
        full.push_str(&idv);
    }
    full
}

fn create_csv(hash: HashMap<String, String>, file_type: String, name: String, path: PathBuf) {
    let keys = hash.keys();
    let mut csv_name = name;
    csv_name.push_str(&file_type);
    csv_name.push_str(".csv");
    let csv_path = create_full_path(path.clone(), csv_name.clone());
    if Path::new(&csv_path).exists() {
        println!(
            "{} exists, Please remove the file from the directory so it is not overwritten",
            csv_name.clone()
        );
    } else {
        let mut wrt = Writer::from_path(csv_path).expect("Did not write csv");
        wrt.write_record(&["info", "misc"]);
        for head in keys {
            let ph = parse_header(head.to_string());
            let hinfo = get_info(ph);

            wrt.write_record(&hinfo).expect("Did not write line");
        }
    }
}

fn create_full_path(tdir: PathBuf, dir: String) -> PathBuf {
    let mut s_dir = tdir.clone().into_os_string().into_string().expect("Nope");
    s_dir.push_str("/");
    s_dir.push_str(&dir);
    PathBuf::from(s_dir)
}

fn get_name(pat: String, search_dir: PathBuf) -> String {
    let mut patmatch = String::new();
    let mut mdir = String::new();
    let mut fdir = String::from("/medtr.");
    if let Ok(files) = read_dir(search_dir.clone()) {
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

    if search_dir
        .clone()
        .into_os_string()
        .into_string()
        .unwrap()
        .split('/')
        .last()
        .unwrap()
        == String::from("genomes")
    {
        fdir.push_str(&mdir);
        fdir.push_str(".genome_main.fna.gz");
    }

    if search_dir
        .clone()
        .into_os_string()
        .into_string()
        .unwrap()
        .split('/')
        .last()
        .unwrap()
        == String::from("annotations")
    {
        fdir.push_str(&mdir);
        fdir.push_str(".cds.fna.gz");
    }

    let mut rvec = String::new();
    rvec.push_str(&mdir);
    rvec.push_str(&fdir);
    rvec
}

fn read_fasta<P>(filename: P) -> HashMap<String, String>
where
    P: AsRef<Path>,
{
    let start = Instant::now();
    let file = File::open(filename).expect("Could not open file");
    let gz = MultiGzDecoder::new(file);
    let buf = BufReader::new(gz);
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
            // println!("{:?}", &line[..].trim());
            curid = line[..].trim().to_string();
        } else {
            curseq.push_str(line.trim());
        }
    }
    let duration = start.elapsed();
    println!("It took {:?} to decode and read", duration);
    fasta
}

fn parse_header(header: String) -> Vec<String> {
    let sp = header.split(' ');
    let mut definition_body: String = String::new();
    let mut svec: Vec<String> = Vec::new();

    for i in sp {
        /* finish getting the definition into one element */
        if !i.contains("=") && !i.contains(">") {
            definition_body.push_str(i);
        } else {
            svec.push(i.to_string());
        }
    }
    let mut definition: String = String::new();
    for i in &svec {
        if i.contains("def=") {
            definition.push_str(&i);
        }
    }
    definition.push_str(&definition_body);
    if let Some(ind) = svec.iter().position(|val| val.contains("def=")) {
        svec.swap_remove(ind);
    }
    svec.push(definition);
    svec
}

fn get_info(header: Vec<String>) -> Vec<String> {
    let mut info: String = String::new();
    let mut misc: String = String::new();
    let mut all = Vec::new();

    for idv in header {
        let mut hsp = idv.split("=");
        if idv.starts_with('>') {
            info.push_str("id=");
            info.push_str(&idv);
            info.push_str(" ");
        } else {
            match hsp.next() {
                Some("gn") => {
                    info.push_str(&idv);
                    info.push_str(" ");
                }
                Some("strand") => {
                    info.push_str(&idv);
                    info.push_str(" ");
                }
                Some("begin") => {
                    info.push_str(&idv);
                    info.push_str(" ");
                }
                Some("end") => {
                    info.push_str(&idv);
                    info.push_str(" ");
                }
                Some("loc") => {
                    info.push_str(&idv);
                    info.push_str(" ");
                }
                Some("len") => {
                    info.push_str(&idv);
                    info.push_str(" ");
                }
                Some(_) => {
                    misc.push_str(&idv);
                    misc.push_str(" ");
                }
                None => println!("Or this one"),
            }
        }
    }
    all.push(info.trim().to_string());
    all.push(misc.trim().to_string());
    all
}
