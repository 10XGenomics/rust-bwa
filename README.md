# Rust-bwa 

Rust-bwa is a Rust wrapper of the BWA api. Pass read-pair information in, and get Rust-htslib BAM records back.
Get started quickly with default settings & a reasonable paired-end model. See docs for more details on customing
parameters or the paired-end model.

```
extern crate bwa;
use bwa::BwaAligner;

fn main() {
    let bwa = BwaAligner::from_path(&"tests/test_ref.fa").unwrap();

    let r1 = b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGG";
    let q1 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";
    let r2 = b"TGCTGCGTAGCAGATCGACCCAGGCATTCCCTAGCGTGCTCATGCTCTGGCTGGTAAACGCACGGATGAGGGCAAAAAT";
    let q2 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";

    let (r1_alns, _r2_alns) = bwa.align_read_pair(b"read_name", r1, q1, r2, q2).unwrap();
    println!("r1 mapping -- tid: {}, pos: {}", r1_alns[0].tid(), r1_alns[0].pos());
}
```

Pre-built rust bindings were generated using `bindgen` for linux using the command:

```
~/.cargo/bin/bindgen --no-doc-comments --whitelist-function mem_align1_core --whitelist-function mem_sam_pe --whitelist-function mem_opt_init --whitelist-function bwa_idx_load --whitelist-function bwa_idx_destroy --whitelist-function mem_process_seq_pe --whitelist-function bwa_fill_scmat --whitelist-var "BWA_IDX_.*" wrapper.h -o linux_prebuilt_bindings.rs
```

`bindgen` can be installed using `cargo install bindgen`. See the documentation [here](https://rust-lang.github.io/rust-bindgen/command-line-usage.html).
