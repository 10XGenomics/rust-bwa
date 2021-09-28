// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.

//! Rust-bwa provides a simple API wrapper around the BWA aligner.
//! Pass read-pair information in, and get Rust-htslib BAM records
//! back.
//!
//! ```
//! use bwa::BwaAligner;
//!
//! let bwa = BwaAligner::from_path(&"tests/test_ref.fa").unwrap();
//!
//! let r1 = b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGG";
//! let q1 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";
//! let r2 = b"TGCTGCGTAGCAGATCGACCCAGGCATTCCCTAGCGTGCTCATGCTCTGGCTGGTAAACGCACGGATGAGGGCAAAAAT";
//! let q2 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";
//!
//! let (r1_alns, _r2_alns) = bwa.align_read_pair(b"read_name", r1, q1, r2, q2);
//! println!("r1 mapping -- tid: {}, pos: {}", r1_alns[0].tid(), r1_alns[0].pos());
//! ```

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

extern crate libc;
extern crate rust_htslib;

extern crate thiserror;

use std::ffi::{CStr, CString};
use std::path::Path;
use std::ptr;
use std::sync::{Arc, Mutex};

use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::HeaderView;

// include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

/// BWA settings object. Currently only default settings are enabled
pub struct BwaSettings {
    bwa_settings: bwa_sys::mem_opt_t,
}

impl BwaSettings {
    /// Create a `BwaSettings` object with default BWA parameters
    pub fn new() -> BwaSettings {
        let ptr = unsafe { bwa_sys::mem_opt_init() };
        let bwa_settings = unsafe { *ptr };
        unsafe { libc::free(ptr as *mut libc::c_void) };
        BwaSettings { bwa_settings }
    }

    /// Set alignment scores
    pub fn set_scores(
        mut self,
        matchp: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> BwaSettings {
        self.bwa_settings.a = matchp;
        self.bwa_settings.b = mismatch;
        self.bwa_settings.o_del = gap_open;
        self.bwa_settings.o_ins = gap_open;
        self.bwa_settings.e_del = gap_extend;
        self.bwa_settings.e_ins = gap_extend;

        unsafe {
            bwa_sys::bwa_fill_scmat(matchp, mismatch, self.bwa_settings.mat.as_mut_ptr());
        }
        self
    }

    /// Set clipping score penalties
    pub fn set_clip_scores(mut self, clip5: i32, clip3: i32) -> BwaSettings {
        self.bwa_settings.pen_clip5 = clip5;
        self.bwa_settings.pen_clip3 = clip3;
        self
    }

    /// Set unpaired read penalty
    pub fn set_unpaired(mut self, unpaired: i32) -> BwaSettings {
        self.bwa_settings.pen_unpaired = unpaired;
        self
    }

    /// Mark shorter splits as secondary
    pub fn set_no_multi(mut self) -> BwaSettings {
        self.bwa_settings.flag |= 0x10; // MEM_F_NO_MULTI
        self
    }
}

#[derive(Debug, thiserror::Error)]
#[error("{0}")]
pub struct ReferenceError(String);

/// A BWA reference object to perform alignments to.
/// Must be loaded from a BWA index created with `bwa index`
pub struct BwaReference {
    bwt_data: *const bwa_sys::bwaidx_t,
    contig_names: Vec<String>,
    contig_lengths: Vec<usize>,
}
unsafe impl Sync for BwaReference {}

impl BwaReference {
    /// Load a BWA reference from disk. Pass the fasta filename of the
    /// original reference as `path`
    pub fn open<P: AsRef<Path>>(path: P) -> Result<BwaReference, ReferenceError> {
        let idx_file = CString::new(path.as_ref().to_str().unwrap()).unwrap();
        let idx = unsafe { bwa_sys::bwa_idx_load(idx_file.as_ptr(), 0x7 as i32) }; // FIXME -- use BWA_IDX_ALL

        if idx.is_null() {
            return Err(ReferenceError(format!(
                "couldn't load reference: {:?}",
                path.as_ref()
            )));
        }

        let mut contig_names = Vec::new();
        let mut contig_lengths = Vec::new();
        let num_contigs = unsafe { (*(*idx).bns).n_seqs };

        for i in 0..num_contigs as isize {
            unsafe {
                let name = CStr::from_ptr((*(*(*idx).bns).anns.offset(i)).name);
                let sz = (*(*(*idx).bns).anns.offset(i)).len;

                let name_string = name.to_owned().into_string().unwrap();
                contig_names.push(name_string);
                contig_lengths.push(sz as usize)
            }
        }

        Ok(BwaReference {
            bwt_data: idx,
            contig_names,
            contig_lengths,
        })
    }

    pub fn create_bam_header(&self) -> Header {
        let mut header = Header::new();
        self.populate_bam_header(&mut header);
        header
    }

    pub fn populate_bam_header(&self, header: &mut Header) {
        for (ref contig_name, &len) in self.contig_names.iter().zip(self.contig_lengths.iter()) {
            add_ref_to_bam_header(header, &contig_name, len);
        }
    }
}

impl Drop for BwaReference {
    fn drop(&mut self) {
        unsafe {
            bwa_sys::bwa_idx_destroy(self.bwt_data as *mut bwa_sys::bwaidx_t);
        }
    }
}

fn add_ref_to_bam_header(header: &mut Header, seq_name: &str, seq_len: usize) {
    let mut header_rec = HeaderRecord::new(b"SQ");
    header_rec.push_tag(b"SN", &seq_name);
    header_rec.push_tag(b"LN", &seq_len);
    header.push_record(&header_rec);
}

/// Paired-end statistics structure used by BWA to score paired-end reads
pub struct PairedEndStats {
    inner: [bwa_sys::mem_pestat_t; 4],
}

impl PairedEndStats {
    /// Generate a 'simple' paired-end read structure that standard forward-reverse
    /// pairs as created by TruSeq, Nextera, or Chromium Genome sample preparations.
    pub fn simple(avg: f64, std: f64, low: i32, high: i32) -> PairedEndStats {
        let pe_stat_null = || bwa_sys::mem_pestat_t {
            failed: 1,
            low: 0,
            high: 0,
            avg: 0.0,
            std: 100.0,
        };

        let pes = [
            pe_stat_null(),
            bwa_sys::mem_pestat_t {
                failed: 0,
                low,
                high,
                avg,
                std,
            },
            pe_stat_null(),
            pe_stat_null(),
        ];

        PairedEndStats { inner: pes }
    }

    pub fn default() -> PairedEndStats {
        Self::simple(200.0, 100.0, 35, 600)
    }
}

/// A BWA aligner. Carries everything required to align
/// reads to a reference and generate BAM records.
pub struct BwaAligner {
    reference: BwaReference,
    header_view: Arc<Mutex<HeaderView>>,
    settings: BwaSettings,
    pe_stats: PairedEndStats,
}
// this is not automatically derived because of an interior
//   mutable pointer inside HeaderView. It _is_ mutated
//   by the Record::from_sam function, so guard it with a mutex
unsafe impl Sync for BwaAligner {}

impl BwaAligner {
    /// Load a BWA reference from the given path and use default BWA settings and paired-end structure.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<BwaAligner, ReferenceError> {
        let bwa_ref = BwaReference::open(path)?;
        Ok(BwaAligner::new(
            bwa_ref,
            BwaSettings::new(),
            PairedEndStats::default(),
        ))
    }

    pub fn new(
        reference: BwaReference,
        settings: BwaSettings,
        pe_stats: PairedEndStats,
    ) -> BwaAligner {
        let header = reference.create_bam_header();
        let header_view = Arc::new(Mutex::new(HeaderView::from_header(&header)));

        BwaAligner {
            reference,
            header_view,
            settings,
            pe_stats,
        }
    }

    /// Align a read-pair to the reference.
    pub fn align_read_pair(
        &self,
        name: &[u8],
        r1: &[u8],
        q1: &[u8],
        r2: &[u8],
        q2: &[u8],
    ) -> (Vec<Record>, Vec<Record>) {
        let name = CString::new(name).unwrap();
        let raw_name = name.into_raw();

        // Prep input data -- need to make copy of reads since BWA will edit the strings in-place
        // FIXME - set an id -- used for a random hash
        let mut r1 = Vec::from(r1);
        let mut q1 = Vec::from(q1);
        let mut r2 = Vec::from(r2);
        let mut q2 = Vec::from(q2);

        let read1 = bwa_sys::bseq1_t {
            l_seq: r1.len() as i32,
            name: raw_name,
            seq: r1.as_mut_ptr() as *mut i8,
            qual: q1.as_mut_ptr() as *mut i8,
            comment: ptr::null_mut(),
            id: 0,
            sam: ptr::null_mut(),
        };

        let read2 = bwa_sys::bseq1_t {
            l_seq: r2.len() as i32,
            name: raw_name,
            seq: r2.as_mut_ptr() as *mut i8,
            qual: q2.as_mut_ptr() as *mut i8,
            comment: ptr::null_mut(),
            id: 0,
            sam: ptr::null_mut(),
        };

        let mut reads = [read1, read2];

        // Align the read pair. BWA will write the SAM data back to the bwa_sys::bseq1_t.sam field
        unsafe {
            let r = *(self.reference.bwt_data);
            let settings = self.settings.bwa_settings;
            bwa_sys::mem_process_seq_pe(
                &settings,
                r.bwt,
                r.bns,
                r.pac,
                reads.as_mut_ptr(),
                self.pe_stats.inner.as_ptr(),
            );
            let _ = CString::from_raw(raw_name);
        }

        // Parse the results from the SAM output & convert the htslib Records
        let sam1 = unsafe { CStr::from_ptr(reads[0].sam) };
        let sam2 = unsafe { CStr::from_ptr(reads[1].sam) };

        let recs1 = self.parse_sam_to_records(sam1.to_bytes());
        let recs2 = self.parse_sam_to_records(sam2.to_bytes());

        unsafe {
            libc::free(reads[0].sam as *mut libc::c_void);
            libc::free(reads[1].sam as *mut libc::c_void);
        }

        (recs1, recs2)
    }

    fn parse_sam_to_records(&self, sam: &[u8]) -> Vec<Record> {
        let mut records = Vec::new();

        for slc in sam.split(|x| *x == b'\n') {
            if slc.len() > 0 {
                let record = {
                    let header_view = self.header_view.lock().unwrap();
                    Record::from_sam(&header_view, slc).unwrap()
                };
                records.push(record);
            }
        }

        records
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn load_aligner() -> BwaAligner {
        let aln = BwaAligner::from_path("tests/test_ref.fa");
        aln.unwrap()
    }

    #[test]
    fn test_load_aligner() {
        let _ = load_aligner();
    }

    fn read_simple() -> [&'static [u8]; 5] {
        let name: &[u8] = b"@chr_727436_727956_3:0:0_1:0:0_0/1";
        let r1  : &[u8] = b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGGCTGGCGCGGCTGATTAATGACATTCCTCTTCCCGGTACAACGGGCGTTGAGCGCGAACTTTTTCGCGCACT";
        let q1  : &[u8] = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        let r2  : &[u8] = b"TGCTGCGTAGCAGATCGACCCAGGCATTCCCTAGCGTGCTCATGCTCTGGCTGGTAAACGCACGGATGAGGGCAAAAATCACCGCAATCCCGCTGGCGGCAGAAAGAAAGTTTTGCACCGTTAAGCCCGCCATCTGGCTGAAATAGCTCA";
        let q2  : &[u8] = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        [name, r1, q1, r2, q2]
    }

    fn read_split() -> [&'static [u8]; 5] {
        let name = b"@chr_1561275_1561756_1:0:0_2:0:0_5c/1";
        let r1 = b"GCATCGATAAGCAGGTCAAATTCTCCCGTCATTATCACCTCTGCTACTTAAATTTCCCGCTTTATAAGCCGATTACGGCCTGGCATTACCCTATCCATAATTTAGGTGGGATGCCCGGTGCGTGGTTGGCAGATCCGCTGTTCTTTATTT";
        let q1 = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        let r2 = b"TCATCGACCCAGGTATCATCGCGACGGGTACGATTACTGGCGAAGGTGAGAATGTTTAAAATCCAGCCGCCGAGTTTTTCAGCAATGGTCACCCATGACCAACCGGTGAACAACGTGAGGGCCGCTGCCCAAACGCATAGCAGCGCAATA";
        let q2 = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        [name, r1, q1, r2, q2]
    }

    fn align_read(r: [&[u8]; 5]) -> (Vec<Record>, Vec<Record>) {
        let bwa = load_aligner();
        bwa.align_read_pair(r[0], r[1], r[2], r[3], r[4])
    }

    #[test]
    fn simple_align() {
        let (r1, r2) = align_read(read_simple());
        assert_eq!(r1[0].pos(), 727806);
        assert_eq!(r2[0].pos(), 727435);
    }

    #[test]
    fn split_align() {
        let (r1, r2) = align_read(read_split());
        assert_eq!(r1.len(), 2);
        assert_eq!(r1[0].pos(), 931375);
        assert_eq!(r1[1].pos(), 932605);
        assert_eq!(r2[0].pos(), 932937);
    }

    #[test]
    fn header() {
        let reference = BwaReference::open("tests/test_ref.fa").unwrap();
        let hdr = b"@SQ\tSN:PhiX\tLN:5386\n@SQ\tSN:chr\tLN:4639675";
        assert_eq!(
            reference.create_bam_header().to_bytes().as_slice(),
            &hdr[..]
        );
    }
}
