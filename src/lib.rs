//! Transcript Integrity Number (TIN) computation from RNA-seq BAM data.
//!
//! TIN measures coverage uniformity across the mRNA body using Shannon entropy.
//! Given per-base coverage `C_i` at sampled mRNA positions:
//! `p_i = C_i / sum(C_j)`, `H = -sum(p_i * ln(p_i))` (non-zero only), `TIN = 100 * exp(H) / l`
//! where `l` is the total number of sampled positions (including zero-coverage positions).
//!
//! Transcripts with fewer than `min_cov` UNIQUE read-start positions are assigned TIN=0.
//! The sample count `n` is halved while `n >= mRNA_len`; when `mRNA_len <= n`, all bases
//! are used. Exon boundary positions are always added to the position set (tin.py behaviour).
//!
//! ## Origin
//!
//! This crate is an independent Rust reimplementation of `tin.py` based on:
//! - Wang L, Wang S, Li W. "`RSeQC`: quality control of RNA-seq experiments."
//!   Bioinformatics. 2012;28(16):2184-5. doi:10.1093/bioinformatics/bts356
//! - Wang L, Nie J, Sicotte H, et al. "Measure transcript integrity using RNA-seq data."
//!   BMC Bioinformatics. 2016;17:58. doi:10.1186/s12859-016-0922-z
//! - The BED12 format specification
//! - Black-box behaviour testing against `RSeQC` `tin.py` 5.0.4
//!
//! No source code from the GPL-v3 upstream was used as reference during
//! implementation. Test fixtures are independently generated.
//!
//! License: MIT OR Apache-2.0.
//! Upstream credit: `RSeQC` <https://rseqc.sourceforge.net/> (GPL-v3).

#![allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]

mod bed;
mod model;
mod sample;
mod score;
mod stats;

pub use model::{SampleReport, TinOutput, TinReport, TinSummary};
pub use score::{compute_tin, compute_tin_score, effective_sample_size};
pub use stats::summarise;
