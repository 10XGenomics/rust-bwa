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

    #[cfg(target_os = "linux")]
    {
        fs::copy("linux_prebuilt_bindings.rs", out_path.join("bindings.rs"))
            .expect("couldn't copy prebuilt bindings");
        println!("cargo:rerun-if-changed=linux_prebuilt_bindings.rs");
    }
}
