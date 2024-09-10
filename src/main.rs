use crate::{pathing::*, read_write::*, search::*};
use csv::Writer;
use std::collections::{hash_map::Keys, HashMap};
use std::io::stdin;
use std::path::Path;
use std::{env, usize};

struct Dirs {
    genome: String,
    annotation: String,
}

fn main() {
    env::set_var("RUST_BACKTRACE", "0");

    let genotype = env::args()
        .nth(1)
        .expect("please enter a valid genome name (arg 1) or type -help");
    if genotype == "-help" {
        println!("
Kbrecondo version v1.0.0

Description: searches fastas for sequences in a window with several search parameters. Hard coded for the Medicago truncatula 
             genome and annotations as formatted in:
             https://data.legumeinfo.org/Medicago/truncatula/

Format: kbrecondo [genotype] [window] [search] [output] [sequence] [species] [option]

Example: 

Manually input search:
    kbrecondo A17 1000 ACTAC test cds medtr -n

Fasta input search:
    kbrecondo A17 1000 search.fasta test cds medtr -f

Flags:
    -n: normal search, runs a search in the format of the manually input search
    -f: use a fasta as the sequence input (currently only takes one sequence but in the future will take mutliple for bulk search)
    -m: asks for a csv containing gene ids to use to guide searching. The first column of the csv must have the gene ids. the ids
        will be matched based on what is found in the headers of the fasta

");
        std::process::exit(3);
    }
    let size_string = env::args()
        .nth(2)
        .expect("please enter the size for the search window (arg 2)");
    let size = size_string.trim().parse::<i32>().unwrap();
    let raw_pattern = env::args()
        .nth(3)
        .expect("please enter the pattern to search for (arg 3)");
    let name = env::args()
        .nth(4)
        .expect("please enter the name for the csv (arg 4)");
    let seq_type = env::args()
        .nth(5)
        .expect("please enter the type of sequence you are searching around (arg 5)");
    let species = env::args()
        .nth(6)
        .expect("please enter the species tag (arg 6)");
    let option = env::args()
        .nth(7)
        .expect("please enter a option tag (arg 7)");

    let top_dir = env::current_dir().expect("bad top dir");

    let genomes = create_full_path(top_dir.clone(), "genomes");
    let annotation = create_full_path(top_dir.clone(), "annotations");

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

    let mut pattern = String::new();
    let mut pat_identifier = String::new();

    if option == "-f" {
        let search_path = create_full_path(top_dir.clone(), raw_pattern.clone());
        assert!(Path::new(&search_path).exists());
        let search_fasta = read_search_fasta_single(search_path);
        let search_key = search_fasta.keys().next().unwrap();
        let fasta_vals = search_fasta.get_key_value(search_key).unwrap().1;
        pat_identifier.push_str(search_key);
        pattern.push_str(fasta_vals);
    }

    if option == "-m" {
        pattern.push_str(&raw_pattern);
        pat_identifier.push_str(&raw_pattern);
    }

    if option == "-n" {
        pattern.push_str(&raw_pattern);
        pat_identifier.push_str(&raw_pattern);
    }

    let pattern_size = pattern.len();
    let test_size = size as usize;
    if test_size < pattern_size {
        println!("Window Size (arg 2) must be larger than the length of the search pattern");
        std::process::exit(3);
    }

    assert!(Path::new(&full_anno).exists());
    assert!(Path::new(&full_geno).exists());

    let decogeno = read_fasta(full_geno);
    let decoanno = read_fasta(full_anno);

    let akeys = decoanno.keys();
    let gkeys = decogeno.keys();

    // std::process::exit(1);
    let mut csv_name = name;
    csv_name.push_str("_");
    csv_name.push_str(&pat_identifier);
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

        println!("Now Searching for {}", &pattern);

        if option == "-m" {
            let mut csv_path = String::new();
            println!("\nPlease enter path to csv");
            stdin()
                .read_line(&mut csv_path)
                .expect("Did not enter a path to csv");
            let csv_path = csv_path.trim();
            let csv_path = create_full_path(top_dir.clone(), csv_path.to_string());
            println!("{:?}", csv_path);
            let search_map = read_csv_first_col(csv_path);
            search(
                decogeno.clone(),
                akeys,
                gkeys,
                search_map.keys(),
                size,
                pattern,
                wrt,
                option,
            );
        } else {
            let search_map: HashMap<String, String> = HashMap::new();
            search(
                decogeno.clone(),
                akeys,
                gkeys,
                search_map.keys(),
                size,
                pattern,
                wrt,
                option,
            );
        }
    }
}

mod pathing {

    use ::std::{fs::read_dir, path::PathBuf};

    pub fn create_full_path(top_dir: PathBuf, dir: &str) -> PathBuf {
        let mut s_dir = top_dir
            .clone()
            .into_os_string()
            .into_string()
            .expect("Nope");
        s_dir.push_str("/");
        s_dir.push_str(dir);
        PathBuf::from(s_dir)
    }

    pub fn get_name(pat: String, search_dir: PathBuf, seq_type: String, species: String) -> String {
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
            let seq_sp: Vec<_> = seq_type.split(" ").collect();
            if seq_sp[1] == "-f" {
                fdir.push_str(&mdir);
                fdir.push_str(".");
                fdir.push_str(&seq_type);
                fdir.push_str(".fna.gz");
            }
            if seq_sp[1] == "-b" {
                fdir.push_str(&mdir);
                fdir.push_str(".");
                fdir.push_str(&seq_type);
                fdir.push_str(".bed.gz");
            }
        }

        let mut rvec = String::new();
        rvec.push_str(&mdir);
        rvec.push_str(&fdir);
        println!("{:?}", rvec);
        rvec
    }
}

mod read_write {

    use csv::{Reader, Writer};
    use flate2::read::MultiGzDecoder;
    use std::{
        collections::HashMap,
        fs::File,
        io::{BufRead, BufReader},
        path::{Path, PathBuf},
        time::Instant,
    };

    struct Reads {
        filename: PathBuf,
    }

    impl Reads {
        pub fn read_csv_first_col<P>(&self) -> HashMap<String, String>
        where
            P: AsRef<Path>,
        {
            let mut id_map: HashMap<String, String> = HashMap::new();
            let mut rdr = Reader::from_path(&self.filename).expect("Csv does not exist");
            for record in rdr.records() {
                let rec = record.expect("268, could not read record");
                let id = &rec[0];
                let sudoname = &rec[1];
                id_map.insert(id.to_string(), sudoname.to_string());
            }

            id_map
        }
    }
    pub fn read_fasta<P>(&self) -> HashMap<String, String>
    where
        P: AsRef<Path>,
    {
        let start = Instant::now();
        let file = File::open(&self.filename).expect("Could not open file");
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

    pub fn write_csv(writer: &mut Writer<File>, info: &Vec<String>, occurances: Vec<String>) {
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
}

mod search {

    use crate::info::*;
    use crate::read_write::write_csv;
    use csv::Writer;
    use std::{
        collections::{hash_map::Keys, HashMap},
        fs::File,
        io::{BufRead, BufReader},
        path::Path,
        time::Instant,
    };

    pub fn search2(
        decogeno: HashMap<String, String>,
        akeys: Keys<'_, String, String>,
        gkeys: Keys<'_, String, String>,
        search_map: Keys<'_, String, String>,

        size: i32,
        pattern: String,
        mut wrt: Writer<File>,
        option: String,
    ) {
        let start = Instant::now();
        for gk in gkeys {
            let geno_header = Key {
                key: gk.to_string(),
                file_type: String::new(),
            };
            println!("\nSearching: {:?}", gk.to_string());
            let info_geno = geno_header.parse_header().get_info();
            let acc = info_geno.get_element("acc");
            let chromosome = acc.split("=").last().unwrap();

            for ak in akeys.clone() {
                let anno_header = Key {
                    key: ak.to_string(),
                    file_type: String::new(),
                };
                let info_anno = anno_header.parse_header().get_info();
                let chrom = info_anno.get_element("chr");
                let spchrom = chrom.split("=").last().unwrap();
                let chrom_len = info_geno.get_element("len");
                let spchrom_len = chrom_len.split("=").last().unwrap();
                let chrom_len_num = spchrom_len.to_string().parse::<i32>().unwrap();

                if option == "-m" {
                    let gn = info_anno.get_element("gn");
                    for search in search_map.clone() {
                        let search_test = search.split("_").last().unwrap();
                        if gn.trim().contains(search_test) && !search_test.is_empty() {
                            if chromosome == spchrom {
                                let chrom_seq = decogeno.get_key_value(gk).unwrap().1;
                                let chrom_seq_len = chrom_seq.len() as i32;
                                assert!(chrom_seq_len == chrom_len_num);
                                let strand = info_anno.get_element("strand");
                                let begend = info_anno.get_begin_end();
                                let window =
                                    build_window(begend[0], begend[1], size, chrom_seq_len);
                                let occrances = search_seq(
                                    chrom_seq.to_string(),
                                    window,
                                    pattern.clone(),
                                    strand,
                                );
                                write_csv(&mut wrt, &info_anno.info, occrances);
                            }
                        } else {
                            continue;
                        }
                    }
                } else {
                    if chromosome == spchrom {
                        let chrom_seq = decogeno.get_key_value(gk).unwrap().1;
                        let chrom_seq_len = chrom_seq.len() as i32;
                        assert!(chrom_seq_len == chrom_len_num);
                        let strand = info_anno.get_element("strand");
                        let begend = info_anno.get_begin_end();
                        let window = build_window(begend[0], begend[1], size, chrom_seq_len);
                        let occrances =
                            search_seq(chrom_seq.to_string(), window, pattern.clone(), strand);
                        write_csv(&mut wrt, &info_anno.info, occrances);
                    }
                }
            }
            let dur = start.elapsed();
            let dur_min = dur.as_secs() / 60;
            let dur_rem = dur.as_secs() % 60;
            println!(
                "Completed in {:?} minutes and {:?} seconds",
                dur_min, dur_rem
            )
        }
    }

    pub fn build_window(begin: i32, end: i32, window_size: i32, seq_len: i32) -> Vec<i32> {
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

    pub fn minus_strand_invsersion(pat: String) -> String {
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

    pub fn read_search_fasta_single<P>(filename: P) -> HashMap<String, String>
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
                curid = line[..].trim().to_string();
            } else {
                curseq.push_str(line.trim());
            }
        }
        fasta.insert(curid.clone(), curseq.clone());
        fasta
    }

    pub fn search_seq(
        seq: String,
        window: Vec<i32>,
        pattern: String,
        strand: String,
    ) -> Vec<String> {
        let strand = strand.split("=").last().unwrap();
        let mut occurances: Vec<String> = Vec::new();
        if seq.len() > pattern.len() {
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
        } else {
            occurances.push(String::from("pattern larger than sequence"));
        }

        occurances
    }
}

mod info {

    use std::{fs::read_dir, path::PathBuf};

    pub struct Key {
        pub key: String,
        pub file_type: String,
    }

    struct Info {
        pub info: Vec<String>,
    }

    impl Key {
        pub fn parse_header(&self) -> Info {
            let mut svec: Vec<String> = Vec::new();
            if self.file_type == "-f" {
                let sp = self.key.split(' ');
                let mut definition_body: String = String::new();

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
            }
            if self.file_type == "-b" {
                // MAKE BED THING HERE OR USE CONVERSION TOOL ATOB FUNCITONALITY
            }

            Info { info: svec }
        }
    }

    impl Info {
        pub fn get_begin_end(&self) -> Vec<i32> {
            let info_split = self.info[0].split(" ");
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
            assert!(begend.len() == 2, "{:?} Did not pass", self.info[0]);
            assert!(begend[0] < begend[1]);
            begend
        }

        pub fn get_info(&self) -> Info {
            let mut info: String = String::new();
            let mut misc: String = String::new();
            let mut all = Vec::new();

            for idv in &self.info {
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
            Info { info: all }
        }

        pub fn get_element(&self, element: &str) -> String {
            let info = &self.info[0];
            let info_sp = info.split(" ");
            let mut element_type = String::new();
            for i in info_sp {
                if i.contains(element) {
                    element_type.push_str(i);
                }
            }
            element_type
        }
    }
}
