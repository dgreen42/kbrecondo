use csv::Writer;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::Instant;
use std::{env, usize};

fn main() {
    env::set_var("RUST_BACKTRACE", "0");

    let start = Instant::now();
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
    let name = env::args()
        .nth(4)
        .expect("please endter the name for the csv (arg 4)");
    let seq_type = env::args()
        .nth(5)
        .expect("please endter the type of sequence you are searching around (arg 5)");
    let species = env::args()
        .nth(6)
        .expect("please endter the species tag (arg 6)");
    let top_dir = env::current_dir().expect("bad top dir");
    let genomes = create_full_path(top_dir.clone(), String::from("genomes"));
    let annotation = create_full_path(top_dir.clone(), String::from("annotations"));
    let dir_geno = get_name(
        genotype.clone(),
        genomes.clone(),
        seq_type.clone(),
        species.clone(),
    );
    let dir_anno = get_name(
        genotype.clone(),
        annotation.clone(),
        seq_type.clone(),
        species.clone(),
    );
    let full_geno = create_full_path(genomes.clone(), dir_geno.clone());
    let full_anno = create_full_path(annotation.clone(), dir_anno.clone());
    println!("{:?}", full_geno);
    println!("{:?}", full_anno);
    assert!(Path::new(&full_anno).exists());
    assert!(Path::new(&full_geno).exists());

    let decogeno = read_fasta(full_geno);
    let decoanno = read_fasta(full_anno);

    let akeys = decoanno.keys();
    let gkeys = decogeno.keys();

    // std::process::exit(1);
    let mut csv_name = name;
    csv_name.push_str("_");
    csv_name.push_str(&pattern);
    csv_name.push_str("_");
    csv_name.push_str(&seq_type);
    csv_name.push_str("_");
    csv_name.push_str(&species);
    csv_name.push_str(".csv");
    let csv_path = create_full_path(top_dir.clone(), csv_name.clone());

    if Path::new(&csv_path).exists() {
        println!(
            "{} exists, Please remove the file from the directory so it is not overwritten",
            csv_name.clone()
        );
    } else {
        let mut wrt = Writer::from_path(csv_path).expect("Did not write csv");
        wrt.write_record(&[
            "id",
            "length",
            "begin",
            "end",
            "strand",
            "occurance.location",
            "info",
        ])
        .expect("Did not write fist line");

        for gk in gkeys {
            let info_geno = get_info(parse_header(gk.to_string()));
            let acc = get_element(info_geno.clone(), String::from("acc"));
            let chromosome = acc.split("=").last().unwrap();

            for ak in akeys.clone() {
                let info_anno = get_info(parse_header(ak.to_string()));
                let chrom = get_element(info_anno.clone(), String::from("chr"));
                let spchrom = chrom.split("=").last().unwrap();
                // come back to this and see if it makes sense. Not sure that the chromosomes that
                // are used in the ladder assert! are differnet, try to reason about why this test
                // makes sense
                // !! I may have done this because it is taking the length from the actual sequence
                // and comparing it to the length that we are getting from the header.
                // I'm high right now so please come back and check this
                let chrom_len = get_element(info_geno.clone(), String::from("len"));
                let spchrom_len = chrom_len.split("=").last().unwrap();
                let chrom_len_num = spchrom_len.to_string().parse::<i32>().unwrap();

                if chromosome == spchrom {
                    let chrom_seq = decogeno.get_key_value(gk).unwrap().1;
                    let chrom_seq_len = chrom_seq.len() as i32;
                    assert!(chrom_seq_len == chrom_len_num);
                    let header = &info_anno[0];
                    let strand = get_element(info_anno.clone(), String::from("strand"));
                    let begend = get_begin_end(info_anno.clone());
                    let window = build_window(begend[0], begend[1], size, chrom_seq_len);
                    let occrances =
                        search_seq(chrom_seq.to_string(), window, pattern.clone(), strand);
                    write_csv(&mut wrt, info_anno.clone(), occrances);
                }
            }
            let dur = start.elapsed();
            let dur_minutes = dur.as_secs() / 60;
            let dur_remainder = dur.as_secs() % 60;
            println!(
                "Time to complete {} search: {:?}, minutes and {:?} seconds",
                gk.to_string(),
                dur_minutes,
                dur_remainder
            );
        }
    }

    let duration = start.elapsed();
    let duration_minutes = duration.as_secs() / 60;
    let duration_remainder = duration.as_secs() % 60;
    println!(
        "Total time to completion: {:?}, minutes and {:?} seconds",
        duration_minutes, duration_remainder
    );
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

fn build_window(begin: i32, end: i32, window_size: i32, seq_len: i32) -> Vec<i32> {
    let mut window: Vec<i32> = Vec::new();
    if begin < window_size {
        let left_bound: i32 = (begin + 1) - begin;
        window.push(left_bound);
        if end + window_size > seq_len {
            let window_end_diff = (end + window_size) - seq_len;
            let right_bound: i32 = (end + window_size) - window_end_diff;
            window.push(right_bound);
        } else {
            let right_bound: i32 = end + window_size;
            window.push(right_bound);
        }
    } else {
        let left_bound: i32 = begin - window_size;
        window.push(left_bound);
        if end + window_size > seq_len {
            let window_end_diff = (end + window_size) - seq_len;
            let right_bound: i32 = (end + window_size) - window_end_diff;
            window.push(right_bound);
        } else {
            let right_bound: i32 = end + window_size;
            window.push(right_bound);
        }
    }
    window
}

fn get_element(info: Vec<String>, element: String) -> String {
    let info = &info[0];
    let info_sp = info.split(" ");
    let mut element_type = String::new();
    for i in info_sp {
        if i.contains(&element) {
            element_type.push_str(i);
        }
    }
    element_type
}

fn minus_strand_invsersion(pat: String) -> String {
    let mut inversion = String::new();
    let pat_inv: String = pat.chars().rev().collect();
    // instead of inveerting the entire sequence, just invert the and reverse it
    for nuc in pat_inv.chars() {
        match nuc {
            'a' => inversion.push_str("t"),
            'c' => inversion.push_str("g"),
            't' => inversion.push_str("a"),
            'g' => inversion.push_str("c"),
            'n' => inversion.push_str("n"),
            _ => println!("This is not a nucleotide, or n"),
        }
    }

    inversion
}

fn search_seq(seq: String, window: Vec<i32>, pattern: String, strand: String) -> Vec<String> {
    let strand = strand.split("=").last().unwrap();
    let mut occurances: Vec<String> = Vec::new();
    if strand == "-" {
        let seq = seq.to_lowercase();
        let pattern = pattern.to_lowercase();
        let inversion = minus_strand_invsersion(pattern);
        let pattern = inversion;
        let bseq = seq.as_bytes();
        let bpat = pattern.as_bytes();
        let left_bound = window[0] as usize;
        let right_bound = window[1] as usize;
        let search_area = &bseq[left_bound..right_bound];

        for i in 0..search_area.len() - bpat.len() {
            if bpat == &search_area[i..i + bpat.len()] {
                let location = left_bound + i;
                occurances.push(location.to_string());
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
                occurances.push(location.to_string());
            }
        }
    }
    occurances
}

fn write_csv(writer: &mut Writer<File>, info: Vec<String>, occurances: Vec<String>) {
    let mut id = String::new();
    let mut strand = String::new();
    let mut length = String::new();
    let mut begin = String::new();
    let mut end = String::new();
    let def_info = &info[1];
    let mut misc = String::new();
    let info = info[0].split(" ");
    for i in info {
        let sp1 = i.split("=").next().unwrap();
        let sp2 = i.split("=").last().unwrap();
        match sp1 {
            "begin" => begin.push_str(sp2),
            "end" => end.push_str(sp2),
            "len" => length.push_str(sp2),
            "strand" => strand.push_str(sp2),
            "id" => id.push_str(sp2),
            _ => misc.push_str(sp2),
        }
    }
    let mut record: Vec<String> = Vec::new();
    for i in occurances.iter() {
        record.push(id.clone());
        record.push(length.clone());
        record.push(begin.clone());
        record.push(end.clone());
        record.push(strand.clone());
        record.push(i.to_string());
        record.push(def_info.to_string());
        writer
            .write_record(record.clone())
            .expect("Did not write record");
        record.clear();
    }
}

fn create_full_path(tdir: PathBuf, dir: String) -> PathBuf {
    let mut s_dir = tdir.clone().into_os_string().into_string().expect("Nope");
    s_dir.push_str("/");
    s_dir.push_str(&dir);
    PathBuf::from(s_dir)
}

fn get_name(pat: String, search_dir: PathBuf, seq_type: String, species: String) -> String {
    let mut patmatch = String::new();
    let mut mdir = String::new();
    let mut fdir = String::from("/");
    fdir.push_str(&species);
    fdir.push_str(".");
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
        fdir.push_str(".");
        fdir.push_str(&seq_type);
        fdir.push_str(".fna.gz");
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
                Some("chr") => {
                    info.push_str(&idv);
                    info.push_str(" ");
                }
                Some("acc") => {
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
