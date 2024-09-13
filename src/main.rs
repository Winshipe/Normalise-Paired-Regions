use rust_htslib::{bam, bam::Read, bam::record::Record};
use std::fs::{self, File};
use std::env;
//use std::io::LineWriter;
use std::io::Write;
//use std::collections::HashMap;
use std::io::BufReader;
use std::io::BufRead;
use std::path::Path;
//use std::cmp;
use std::vec::Vec;
//use regex::Regex; 
//use regex::Captures; 
//use std::thread;
//use std::sync::Arc;
use rand::seq::SliceRandom;
use clap::{arg, Arg, Command, ArgAction};

fn read_regions(path: &str, ) -> (Vec<(String, String)>){
    let reader = BufReader::new(File::open(path).expect("Cannot open regions file"));
    let mut output = Vec::with_capacity(100);
    for l in reader.lines() {
        let mut line = l.unwrap();
        if (line.starts_with('#') || line.starts_with("@")){
            continue;
        }
        let mut split = line.trim().split("\t");
        let first_region = split.next().expect("Unable to read first region in line");
        let snd_region = split.next().expect("Unable to read second region in line");
        output.push((first_region.to_owned(), snd_region.to_owned()));
    }
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
fn read_records(bam_file3: &mut bam::IndexedReader, bam_file25: &mut bam::IndexedReader, reads3: &mut Vec<Record>, reads25: &mut Vec<Record>) -> () {
    for r in bam_file3.records(){
        let read = r.unwrap();
        if read.is_proper_pair(){
            reads3.push(read);
        }
    }
    for r in bam_file25.records(){
        let read = r.unwrap();
        if read.is_proper_pair(){
            reads25.push(read);
        }
    }
}
fn normalize_region(bam_file3: &mut bam::IndexedReader, bam_file25: &mut bam::IndexedReader, out3: &mut bam::Writer, out25: &mut bam::Writer, region3: &str, region25: &str, min_covg: &f64) -> (){
    let (scaff3, start3, end3) = split_region_str(&region3);
    let (scaff25, start25, end25) = split_region_str(&region25);
    bam_file3.fetch((&scaff3, start3, end3));
    bam_file25.fetch((&scaff25, start25, end25));
    let mut reads3: Vec<Record> = Vec::with_capacity(5000);
    let mut reads25: Vec<Record> = Vec::with_capacity(5000);
    let mut rng = rand::thread_rng();
    read_records(bam_file3, bam_file25, &mut reads3, &mut reads25);
    if reads3.len() > reads25.len() {
        let nreads = reads25.len();// sample 3 down to the same number of reads as 25
        let mut covg = 0.0;
        if nreads > 0 {
            covg = ((nreads * reads25[0].seq_len()) as f64)   / ((end25 - start25) as f64);
        }
        if(covg < *min_covg){
            println!("Skipping {} {} because they are under the minimum coverage ({:.2})",region3, region25, covg);
        } else {
            println!("Region 1 {} Region 2 {}; # reads in 1: {} # reads in 2: {}; downsampling region 1 to {} reads ({:.2}x)", region3, region25, reads3.len(), reads25.len(), nreads, covg);
            for record in reads3.choose_multiple(&mut rng, nreads){ //downsample and write to new bam, comes from rand::
                out3.write(&record);
            }
            for record in reads25 {
                out25.write(&record);
            }
        }
    } else { // same as above but backwards
        let nreads = reads3.len(); //(fraction * (reads25.len() as f64)) as usize; // dont need to compute a fraction when using choose_multiple
        let mut covg = 0.0;
        if nreads > 0 {
            covg = ((nreads * reads3[0].seq_len()) as f64)   / ((end25 - start25) as f64);
        }
        if(covg < *min_covg){
            println!("Skipping {} {} because they are under the minimum coverage ({:.2})",region3, region25, covg);            
        } else {
            println!("Region 1 {} Region 2 {}; # reads in 1 {} # reads in 2 {}; downsampling region 2 to {} reads ({:.2}x)", region3, region25, reads3.len(), reads25.len(), nreads, covg);
            for record in reads25.choose_multiple(&mut rng, nreads){
                out25.write(&record);
            }
            for record in reads3 {
                out3.write(&record);
            }
        }
    }
}
fn normalize_given_regions(path3: &str, path25: &str, outpath3: &str, outpath25: &str, regions: Vec<(String, String)>, min_covg: &f64){
    let mut bam_file3 = bam::IndexedReader::from_path(path3).expect("Cannot open bam 3 file!");
    let mut bam_file25 = bam::IndexedReader::from_path(path25).expect("Cannot open bam 25 file!");
    let mut bam_header3 = bam::Header::from_template(bam_file3.header());
    let mut bam_header25 = bam::Header::from_template(bam_file25.header());
    let mut out3 = bam::Writer::from_path(outpath3, &bam_header3, bam::Format::Bam).expect("Cannot open bam file 3 for writing!");
    let mut out25 = bam::Writer::from_path(outpath25, &bam_header25, bam::Format::Bam).expect("Cannot open bam file 25 for writing!");
    println!("Normalizing Regions...");
    if *min_covg > 0.0 {
        println!("Minimum coverage set to {}",min_covg);
    }
    for (region3, region25 ) in regions{
        normalize_region(&mut bam_file3, &mut bam_file25, &mut out3, &mut out25, &region3, &region25, min_covg);
    }
}
fn build_argparser() -> Command {
    let output = Command::new("Normalize Paired Regions")
        .about("Normalizes coverage between paired regions given by a tab separated file ie chr:1-1000<tab>chr2:1000-2000")
        .arg(Arg::new("bam_1")
            .help("Path to first BAM file")
            .short('1')
            .long("bam_1")
            .value_name("PATH")
            .required(true)
        )
        .arg(Arg::new("bam_2")
            .help("Path to second BAM file")
            .short('2')
            .long("bam_2")
            .value_name("PATH")
            .required(true)
        )
        .arg(Arg::new("regions")
            .help("Path to tab-separated regions file")
            .short('r')
            .long("regions")
            .value_name("PATH")
            .required(true)
        )
        .arg(Arg::new("output_1")
            .help("Output path for 1st BAM file")
            .short('o')
            .long("output_1")
            .value_name("PATH")
            .required(false)
        )
        .arg(Arg::new("output_2")
            .help("Output path for 1st BAM file")
            .short('n')
            .long("output_2")
            .value_name("PATH")
            .required(false)
        )
        .arg(arg!(-c --min_covg "Minimum coverage (default 0)").required(false).action(ArgAction::Set)); //can use macro but harder to debug
        return output;
}
fn main(){
    let argparser_obj = build_argparser();
    let argparser = argparser_obj.get_matches();
    let regions = read_regions(argparser.get_one::<String>("regions").expect("Regions file path not found"));
    let bam1 = argparser.get_one::<String>("bam_1").expect("Bam-1 file path not found");
    let bam2 = argparser.get_one::<String>("bam_2").expect("Bam-2 file path not found");
    let default_1 = "normalized_".to_owned() + Path::new(bam1).file_name().unwrap().to_str().unwrap();
    let out1 = argparser.get_one::<String>("output_1").unwrap_or(&default_1);
    let default_2 = "normalized_".to_owned() + Path::new(bam2).file_name().unwrap().to_str().unwrap();
    let out2 = argparser.get_one::<String>("output_2").unwrap_or(&default_2);
    let default_cvg = "0.0".to_string();
    let min_covg_str = argparser.get_one::<String>("min_covg").unwrap_or(&default_cvg);
    let min_covg = min_covg_str.parse::<f64>().expect("Failed to parse minimum coverage");
    normalize_given_regions(bam1,bam2,out1,out2,regions, &min_covg);
}