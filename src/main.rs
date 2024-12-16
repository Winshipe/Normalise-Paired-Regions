use rust_htslib::{bam, bam::Read, bam::record::Record};
use std::fs::{self, File};
use std::env;
//use std::io::LineWriter;
use std::io::Write;
use std::collections::HashMap;
use std::io::BufReader;
use std::io::BufRead;
use std::path::Path;
use itertools::Itertools;
use std::vec::Vec;
use std::cmp::max;
use std::cmp::min;
use rand::seq::index::sample;
use rand::seq::index::IndexVec;
//use rand::seq::SliceRandom;
use clap::{arg, Arg, Command, ArgAction};
use std::iter::zip;

struct Params {
    min_covg: f64,
    allow_unpaired: bool,
    verbose: bool
}

impl Default for Params {
    fn default() -> Params {
        Params { //these come from instrain
            min_covg: 5.0,
            allow_unpaired: false,
            verbose: false
        }

    }
}
fn read_regions(path: &str, ) -> (Vec<Vec<String>>){
    let reader = BufReader::new(File::open(path).expect("Cannot open regions file"));
    let mut output = Vec::with_capacity(100);
    for l in reader.lines() {
        let mut line = l.unwrap();
        if (line.starts_with('#') || line.starts_with("@")){
            continue;
        }
        let mut split = line.trim().split("\t");
        let mut split_line: Vec<String> = split.map(|s| s.to_string()).collect(); //necessary to secure ownership of vec contents before returning
        output.push(split_line.to_owned());
    }
    output.sort();
    return output;
}
fn split_region_str(region: &str) -> (String, i64, i64){
    // region = scaffold:start-end
    let mut split = region.clone().split(":");
    let scaff = split.next().unwrap();
    let mut start_end = split.next().unwrap();
    split = start_end.split("-");
    let start = split.next().unwrap().parse::<i64>().unwrap();
    let end = split.next().unwrap().parse::<i64>().unwrap();
    return (scaff.to_string(), start, end);
}
fn read_records_as_pairs(bam_file: &mut bam::IndexedReader, reads: &mut HashMap<String, (Record, Option<Record>)>) -> () {
    for r in bam_file.records(){
        let read = r.unwrap();
        if read.is_proper_pair(){
            let rname =  String::from_utf8_lossy(read.qname()).into_owned();
            if reads.contains_key(&rname){
                reads.get_mut(&rname).unwrap().1 = Some(read);
            } else {
                reads.insert(rname, (read, None));
            }
        }
    }
}
fn read_records(bam_file: &mut bam::IndexedReader, reads: &mut Vec<Record>) -> () {
    for r in bam_file.records(){
        let read = r.unwrap();
        if read.is_proper_pair(){
            reads.push(read);
        }
    }
}
fn complain_and_return_0(scaff: &String, start: i64, end: i64) -> usize {
    eprintln!("\x1b[93mUnable to read {}:{}-{}!\x1b[0m",scaff, start, end);
    return 0;
}
fn unwrap_one_more_layer(result: &mut rust_htslib::errors::Result<bam::Record>, scaff: &String, start: i64, end: i64) -> usize{
    //println!("One more layer");
    let output = match result {
        Err(r) => complain_and_return_0(scaff, start, end),
        Ok(r) => r.seq_len() //friggin' finally
    };
    return output;
}
/*fn get_read_len(bam_file: &mut bam::IndexedReader, scaff: &String, start: i64, end: i64) -> usize {
    //this thing is packaged as an Option<Result<type>> so gotta do some work to unwrap it
    bam_file.fetch((&scaff, start, end));
    let possible_read = bam_file.records().next();
    let read_len = match possible_read {
        None => complain_and_return_0(scaff, start, end),
        Some(mut r) => unwrap_one_more_layer(&mut r, scaff, start, end)
    };
    return read_len;
}*/
fn sum_read_lengths (bam_file: &mut bam::IndexedReader, scaff: &String, start: i64, end: i64) -> f64 {
    bam_file.fetch((&scaff, start, end));
    return bam_file.rc_records()
    .map(|x| x.expect("Failure parsing Bam file"))
    .filter(|read|read.is_proper_pair())
    .fold(0 , |summed, read| summed + read.seq_len()) as f64;
}

fn find_coverages(bam_files: &mut Vec<bam::IndexedReader>, split_regions: & Vec<(String, i64, i64)>, params: &Params) -> Vec<f64>{
    let mut coverages = Vec::<f64>::with_capacity(split_regions.len());
    for i in 0..split_regions.len(){
        let (scaff, start, end) = split_regions[i].clone();
        bam_files[i].fetch((&scaff, start, end));
        let mut read_count = 0.0; 
        let mut read_names = bam_files[i].records()
            .map(|r| r.expect("Failure parsing BAM file"))
            .filter(|r| r.is_proper_pair() && !r.is_mate_unmapped())
            .map(|r| String::from_utf8_lossy(r.qname()).into_owned())
            .collect::<Vec<String>>();
        if !params.allow_unpaired {
            read_count = (read_names.iter().duplicates().count() * 2) as f64; 
        } else {
            read_count = read_names.len() as f64;
        }
            
        let reads_per_base = (read_count) / ((end - start) as f64);
        //let read_len = bam_files[i].records().next().unwrap().seq_len() as f64;
        let read_len = sum_read_lengths(&mut bam_files[i], &scaff, start, end) / read_count;
        coverages.push(reads_per_base * read_len);
    }
    return coverages;
}


fn normalize_region_set(bam_files: &mut Vec<bam::IndexedReader>, out_files: &mut Vec<bam::Writer>, regions: &Vec<String>, params: &Params) -> (){
    //normalize one abitrarily large set of regions, e.g. a gene found in 4 different organisms
    let mut rng = rand::thread_rng();
    
    let split_regions: Vec<(String, i64, i64)> = regions.iter().map(|r| split_region_str(&r)).collect();
    let mut coverages = find_coverages(bam_files, &split_regions, params);
    //cant use &coverages.iter().min().unwrap() b/c floats can be NaN which isn't ordered
    let mut smallest_covg = params.min_covg.max(coverages.iter().fold(f64::INFINITY,|min_so_far, &b| min_so_far.min(b))); //folds must start with initial value
    
    for rgn_idx in 0..regions.len(){
        let (scaff, start, end) = split_regions[rgn_idx].clone();  
        let covg = coverages[rgn_idx];
        let mut covg_ratio = smallest_covg / covg;
        if covg_ratio > 1.0 {
            if params.verbose {
                println!("Keeping {}:{}-{} as is because it is under the minimum coverage ({:.2} : {})",scaff, start, end, covg, params.min_covg);
            }
            covg_ratio = 1.0;
        } else if covg_ratio == 1.0 {
            if params.verbose {
                println!("Keeping {}:{}-{} as is  ({:.2} : {:.2})",scaff, start, end, covg, smallest_covg);
            }
        } else {
            if params.verbose {
                println!("Normalizing region {}:{}-{} from {:.2} to {:.2}x ({:.3} fold downscaling)",scaff, start, end, covg, smallest_covg, covg_ratio);
            }
        }
        
        bam_files[rgn_idx].fetch((&scaff, start, end));
        if !params.allow_unpaired {
            let mut pairs_map: HashMap<String, (Record, Option<Record>)> = HashMap::with_capacity(5000);
            read_records_as_pairs(&mut bam_files[rgn_idx], &mut pairs_map);
            let pairs: Vec<(Record, Record)> = pairs_map.values().filter_map(|p| Some(p.0.to_owned()).zip(p.1.to_owned())).collect();
            let nreads = (covg_ratio * (pairs.len() as f64)) as usize;
            if params.verbose {println!("Writing {} pairs (originally {})", nreads, pairs.len());}
            let mut indices: Vec<usize> = sample(&mut rng, pairs.len(), nreads).into_iter().collect();
            
        // println!("{}",nreads);
            for idx in indices{
                out_files[rgn_idx].write(&pairs[idx].0);
                out_files[rgn_idx].write(&pairs[idx].1);
            }
        } else {
            let mut reads: Vec<Record> = Vec::with_capacity(5000);
            read_records(&mut bam_files[rgn_idx], &mut reads);
            let nreads = (covg_ratio * (reads.len() as f64)) as usize;
            if params.verbose {println!("Writing {} reads (originally {})", nreads, reads.len());}
            let mut indices: Vec<usize> = sample(&mut rng, reads.len(), nreads).into_iter().collect();
        // println!("{}",nreads);
            for idx in indices{
                out_files[rgn_idx].write(&reads[idx]);
            }
        }
        
    }
}

fn normalize_given_regions(in_paths: &Vec<String>, out_paths: &Vec<String>, region_sets: &Vec<Vec<String>>, params: &Params){
    let mut in_bam_files = Vec::<bam::IndexedReader>::with_capacity(in_paths.len());
    let mut in_bam_headers = Vec::<bam::Header>::with_capacity(in_paths.len());
    let mut out_bam_files = Vec::<bam::Writer>::with_capacity(in_paths.len());
    
    for i in 0..in_paths.len(){
        println!("Opening {} for reading...",in_paths[i]);
        let bai_path =in_paths[i].clone() + ".bai";
        if !(Path::new(&bai_path).exists()){
            println!("Indexing BAM {}",in_paths[i]);
            rust_htslib::bam::index::build(in_paths[i].clone(),Some(bai_path),  rust_htslib::bam::index::Type::Bai, 1);
        }
        in_bam_files.push(bam::IndexedReader::from_path(in_paths[i].clone()).expect("Cannot open bam file!"));
        in_bam_headers.push(bam::Header::from_template(in_bam_files[i].header()));
        println!("Opening {} for writing...",out_paths[i]);
        out_bam_files.push(bam::Writer::from_path(out_paths[i].clone(), &in_bam_headers[i], bam::Format::Bam).expect("Cannot open bam file for writing!"));
    }
    println!("Normalizing Regions...");
    if params.min_covg > 0.0 {
        println!("Minimum coverage set to {}",params.min_covg);
    }
    let mut i = 0;
    for regions in region_sets{
        normalize_region_set(&mut in_bam_files, &mut out_bam_files, &regions, params);
    }
}
fn build_argparser() -> Command {
    let output = Command::new("Normalize Paired Regions") //CLAP, pretty neat but short codes are very limiting
        .about("Normalizes coverage between regions given by a tab separated file\nie chr:1-1000<tab>chr2:1000-2000<tab>...")
        .arg(Arg::new("regions")
            .help("Path to tab-separated regions file")
            .short('r')
            .long("regions")
            .value_name("PATH")
            .required(true)
        )
        .arg(Arg::new("bams")
            .help("Paths to bam files, must match # of cols in regions file! ./bam_1 ./bam_2 [... ./bam_n]")
            .short('b')
            .long("bams")
            .value_name("PATHS")
            .num_args(1..)
            .required(true)
        )
        .arg(Arg::new("unpaired")
            .help("Allows single reads through")
            .short('u')
            .long("allow_unpaired")
            .action(ArgAction::SetTrue)
            .required(false)
        )
        .arg(Arg::new("verbose")
           .help("prints debugging info")
            .short('v')
            .long("verbose")
            .action(ArgAction::SetTrue)
            .required(false)
        )
        .arg(arg!(-c --min_covg "Minimum coverage (default 0)").required(false).action(ArgAction::Set)); //can use macro but harder to debug
        
        return output;
}
fn main(){
    let argparser_obj = build_argparser();
    let mut argparser = argparser_obj.get_matches();
    let regions = read_regions(argparser.get_one::<String>("regions").expect("Regions file path not found"));
    let bams_paths: Vec<String> = argparser.remove_many("bams").expect("`bams` is required").collect();
    let bam_out_paths: Vec<String> = bams_paths.iter().map(|p| "normalized_".to_owned() + Path::new(p).file_name().unwrap().to_str().unwrap()).collect();
    let default_cvg = "0.0".to_string();
    let min_covg_str = argparser.get_one::<String>("min_covg").unwrap_or(&default_cvg);
    let params = Params {
        verbose: *argparser.get_one::<bool>("verbose").unwrap_or(&false),
        min_covg: min_covg_str.parse::<f64>().expect("Failed to parse minimum coverage"),
        allow_unpaired:  *argparser.get_one::<bool>("unpaired").unwrap_or(&false)
    };
    normalize_given_regions(&bams_paths,&bam_out_paths, &regions, &params);
}