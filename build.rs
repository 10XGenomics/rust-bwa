#[cfg(feature = "bindgen")]
extern crate bindgen;

use std::env;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    Command::new("make")
        .current_dir("bwa")
        .arg("CFLAGS=-g -Wall -O2 -fPIC")
        .arg("libbwa.a")
        .status()
        .ok()
        .expect("Failed to build bwa");

    println!("cargo:rustc-link-search=bwa");
    println!("cargo:rustc-link-lib=static=bwa");
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    // If bindgen is enabled, use it
    #[cfg(feature = "bindgen")]
    {
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
        bindings
            .write_to_file(out_path.join("bindings.rs"))
            .expect("Couldn't write bindings!");
    }

    #[cfg(all(not(feature = "bindgen"), target_os = "linux"))]
    {
        fs::copy("linux_prebuilt_bindings.rs", out_path.join("bindings.rs"))
            .expect("couldn't copy prebuilt bindings");
        println!("cargo:rerun-if-changed=linux_prebuilt_bindings.rs");
    }
}
