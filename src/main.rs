use rust_htslib::{bam, bam::Read, bam::record::Record};
use std::fs::{self, File};
use std::env;
//use std::io::LineWriter;
use std::io::Write;
//use std::collections::HashMap;
use std::io::BufReader;
use std::io::BufRead;
use std::path::Path;

use std::vec::Vec;
use std::cmp::max;
use std::cmp::min;
use rand::seq::index::sample;
use rand::seq::index::IndexVec;
//use rand::seq::SliceRandom;
use clap::{arg, Arg, Command, ArgAction};
use std::iter::zip;

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
    let output = match result {
        Err(r) => complain_and_return_0(scaff, start, end),
        Ok(r) => r.seq_len() //friggin' finally
    };
    return output;
}
fn get_read_len(bam_file: &mut bam::IndexedReader, scaff: &String, start: i64, end: i64) -> usize {
    //this thing is packaged as an Option<Result<type>> so gotta do some work to unwrap it
    let possible_read = bam_file.records().next();
    let read_len = match possible_read {
        None => complain_and_return_0(scaff, start, end),
        Some(mut r) => unwrap_one_more_layer(&mut r, scaff, start, end)
    };
    return read_len;
}
fn normalize_region_set(bam_files: &mut Vec<bam::IndexedReader>, out_files: &mut Vec<bam::Writer>, regions: &Vec<String>, min_covg: f64) -> (){
    //normalize one abitrarily large set of regions, e.g. a gene found in 4 different organisms
    let mut rng = rand::thread_rng();
    let mut coverages = Vec::<f64>::with_capacity(regions.len());
    let split_regions: Vec<(String, i64, i64)> = regions.iter().map(|r| split_region_str(&r)).collect();
    for i in 0..regions.len(){
        let (scaff, start, end) = split_regions[i].clone();
        bam_files[i].fetch((&scaff, start, end));
        let reads_per_base = (bam_files[i].records().count() as f64) / ((end - start) as f64);
        let read_len = get_read_len(&mut bam_files[i], &scaff, start, end) as f64;
        coverages.push(reads_per_base * read_len);
    }
    //cant use &coverages.iter().min().unwrap() b/c floats can be NaN which isn't ordered
    let mut smallest_covg = min_covg.max(coverages.iter().fold(f64::INFINITY,|min_so_far, &b| min_so_far.min(b))); //folds must start with initial value
    
    for i in 0..regions.len(){
        let (scaff, start, end) = split_regions[i].clone();
        let mut reads: Vec<Record> = Vec::with_capacity(5000);
        let covg = &coverages[i];
        bam_files[i].fetch((&scaff, start, end));
        read_records(&mut bam_files[i], &mut reads);
        let covg_ratio = smallest_covg / covg;
        if covg_ratio > 1.0 {
            println!("Skipping {}:{}-{} because it is under the minimum coverage ({:.2} : {})",scaff, start, end, covg, min_covg);
            continue;
        }
        let nreads = (covg_ratio * (reads.len() as f64)) as usize;
        let mut indices: Vec<usize> = sample(&mut rng, reads.len(), nreads).into_iter().collect();
        for idx in indices{

        }
    }
}
/*        if reads_per_base3 > reads_per_base25 {
        //theoretically I could move the if else contents to their own fn but that makes printing the regions
        // in order harder for not much much more readability
        let nreads = (reads3.len() as f64 * (reads_per_base25 / reads_per_base3)) as usize;// sample 3 down to the same number of reads as 25
        let mut covg = 0.0;
        if nreads > 0 {
            covg = ((nreads * reads3[0].seq_len()) as f64)   / ((end3 - start3) as f64);
        }
        if(covg < *min_covg){
            println!("Skipping {} {} because they are under the minimum coverage ({:.2})",region3, region25, covg);
        } else {
            println!("Region 1 {} Region 2 {}; # reads in 1: {} # reads in 2: {}; downsampling region 1 to {} reads ({:.2}x)", region3, region25, reads3.len(), reads25.len(), nreads, covg);
            let mut indices: Vec<usize> = sample(&mut rng, reads25.len(), nreads).into_iter().collect();
            indices.sort();            
            //for record in reads3.choose_multiple(&mut rng, nreads){ //downsample and write to new bam, comes from rand::
            for idx in indices{
                out3.write(&reads3[idx]);
                //out3.write(&record);
            }
            for record in reads25 {
                out25.write(&record);
            }
        }
    } else { // same as above but backwards
        let nreads = (reads25.len() as f64 * (reads_per_base3 / reads_per_base25)) as usize;// sample 3 down to the same number of reads as 25
        let mut covg = 0.0;
        if nreads > 0 {
            covg = ((nreads * reads25[0].seq_len()) as f64)   / ((end25 - start25) as f64);
        }
        if(covg < *min_covg){
            println!("Skipping {} {} because they are under the minimum coverage ({:.2})",region3, region25, covg);            
        } else {
            println!("Region 1 {} Region 2 {}; # reads in 1 {} # reads in 2 {}; downsampling region 2 to {} reads ({:.2}x)", region3, region25, reads3.len(), reads25.len(), nreads, covg);
            let mut indices: Vec<usize> = sample(&mut rng, reads25.len(), nreads).into_iter().collect();
            indices.sort();
            //for record in reads25.choose_multiple(&mut rng, nreads){
            for idx in indices{
                out25.write(&reads25[idx]);
            }
            for record in reads3 {
                out3.write(&record);
            }
        }
    }
}
    */
fn normalize_given_regions(in_paths: &Vec<String>, out_paths: &Vec<String>, region_sets: Vec<Vec<String>>, min_covg: f64){
    /*let mut bam_file3 = bam::IndexedReader::from_path(path3).expect("Cannot open bam 3 file!");
    let mut bam_file25 = bam::IndexedReader::from_path(path25).expect("Cannot open bam 25 file!");
    let mut bam_header3 = bam::Header::from_template(bam_file3.header());
    let mut bam_header25 = bam::Header::from_template(bam_file25.header());
    
    let mut out3 = bam::Writer::from_path(outpath3, &bam_header3, bam::Format::Bam).expect("Cannot open bam file 3 for writing!");
    let mut out25 = bam::Writer::from_path(outpath25, &bam_header25, bam::Format::Bam).expect("Cannot open bam file 25 for writing!");
    */

    let mut in_bam_files = Vec::<bam::IndexedReader>::with_capacity(in_paths.len());
    let mut in_bam_headers = Vec::<bam::Header>::with_capacity(in_paths.len());
    let mut out_bam_files = Vec::<bam::Writer>::with_capacity(in_paths.len());
    let mut i = 0;
    for (in_path, out_path) in zip(in_paths, out_paths){
        println!("Opening {} for reading...",in_path);
        in_bam_files.push(bam::IndexedReader::from_path(in_path).expect("Cannot open bam file!"));
        in_bam_headers.push(bam::Header::from_template(in_bam_files[i].header()));
        println!("Opening {} for writing...",out_path);
        out_bam_files.push(bam::Writer::from_path(out_path, &in_bam_headers[i], bam::Format::Bam).expect("Cannot open bam file for writing!"));
    }
    println!("Normalizing Regions...");
    if min_covg > 0.0 {
        println!("Minimum coverage set to {}",min_covg);
    }

    for regions in region_sets{
        normalize_region_set(&mut in_bam_files, &mut out_bam_files, &regions, min_covg);
    }
}
fn build_argparser() -> Command {
    let output = Command::new("Normalize Paired Regions") //CLAP, pretty neat but short codes are very limiting
        .about("Normalizes coverage between paired regions given by a tab separated file ie chr:1-1000<tab>chr2:1000-2000")
        /* .arg(Arg::new("bam_1")
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
        )*/
        .arg(Arg::new("regions")
            .help("Path to tab-separated regions file")
            .short('r')
            .long("regions")
            .value_name("PATH")
            .required(true)
        )/*
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
        )*/
        .arg(Arg::new("bams")
            .help("Paths to bam files, must match # of cols in regions file! ./bam_1 ./bam_2 [... ./bam_n]")
            .short('b')
            .long("bams")
            .value_name("PATHS")
            .num_args(1..)
            .required(true)
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
    let min_covg = min_covg_str.parse::<f64>().expect("Failed to parse minimum coverage");
    normalize_given_regions(&bams_paths,&bam_out_paths, regions, min_covg);
}