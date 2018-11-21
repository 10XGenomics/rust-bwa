extern crate bindgen;
extern crate skeptic;
use std::path::PathBuf;
use std::process::Command;
use std::env;


fn main() {

    Command::new("make").current_dir("bwa")
                        .arg("CFLAGS=-g -Wall -O2 -fPIC")
                        .arg("libbwa.a")
                        .status().ok().expect("Failed to build bwa");

    
    println!("cargo:rustc-link-search=bwa");
    println!("cargo:rustc-link-lib=static=bwa");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("wrapper.h")
        .whitelist_function("mem_align1_core")
        .whitelist_function("mem_sam_pe")
        .whitelist_function("mem_opt_init")
        .whitelist_function("bwa_idx_load")
        .whitelist_function("bwa_idx_destroy")
        .whitelist_function("mem_process_seq_pe")
        .whitelist_function("bwa_fill_scmat")
        .whitelist_var("BWA_IDX_*")

        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");



    // generates doc tests for `README.md`.
    // skeptic::generate_doc_tests(&["README.md"]);
}
