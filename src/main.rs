use csv::Writer;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::env;
use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn main() {
    let start = Instant::now();
    // let test_top1 = PathBuf::from("/home/david/Documents/Academic/UVM/Harris_Lab/Mt_sequences");
    let test_top2 = PathBuf::from("/media/david/WorkDrive/Documents/UVM/HarrisLab/Mt_Data");
    let genotype = env::args()
        .nth(1)
        .expect("please enter a valid genome name");
    // let top_dir = env::current_dir().expect("bad top dir");
    let genomes = create_full_path(test_top2.clone(), String::from("genomes"));
    let annotation = create_full_path(test_top2.clone(), String::from("annotations"));
    let dir_geno = get_name(genotype.clone(), genomes.clone());
    let dir_anno = get_name(genotype.clone(), annotation.clone());
    let full_geno = create_full_path(genomes.clone(), dir_geno.clone());
    let full_anno = create_full_path(annotation.clone(), dir_anno.clone());
    assert!(Path::new(&full_anno).exists());
    assert!(Path::new(&full_geno).exists());

    let decogeno = read_fasta(full_geno);
    let decoanno = read_fasta(full_anno);

    let full_geno = collapse_lines(decogeno.clone());

    create_csv(
        decogeno,
        String::from("geno"),
        genotype.clone(),
        test_top2.clone(),
    );
    create_csv(
        decoanno,
        String::from("anno"),
        genotype.clone(),
        test_top2.clone(),
    );
    let duration = start.elapsed();
    println!("Total time to completion: {:?}", duration);
}

fn get_begin_end(genome: String, info: Vec<String>) {}

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
    let mut svec: Vec<String> = Vec::new();
    for i in sp {
        svec.push(i.to_string());
    }

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
    all.push(info);
    all.push(misc);
    all
}
