use std::collections::HashSet;
use std::fs::File;
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use rsomics_common::{Result, RsomicsError};

use crate::bed::parse_bed12;
use crate::model::TinOutput;
use crate::sample::pick_positions;

/// Shannon entropy `H = -sum(p_i * ln(p_i))` over non-zero values.
///
/// Returns 0.0 when all values are zero or the list is empty.
fn shannon_entropy(vals: &[f64]) -> f64 {
    let total: f64 = vals.iter().sum();
    if total == 0.0 {
        return 0.0;
    }
    vals.iter()
        .filter(|&&v| v > 0.0)
        .map(|&v| {
            let p = v / total;
            -p * p.ln()
        })
        .sum()
}

/// Compute TIN from a coverage vector and the total position count `l`.
///
/// Only non-zero coverage positions contribute to entropy; `l` is the denominator
/// (total sampled positions including zero-coverage ones), matching `tin.py`.
/// Returns 0.0 when total coverage is zero.
#[must_use]
pub fn compute_tin_score(mrna_cov: &[u32], n: usize) -> f64 {
    let mrna_len = mrna_cov.len();
    if mrna_len == 0 || n == 0 {
        return 0.0;
    }

    let step = mrna_len / n;
    if step == 0 {
        return 0.0;
    }

    let sample: Vec<f64> = (0..n).map(|i| f64::from(mrna_cov[i * step])).collect();
    let total: f64 = sample.iter().sum();
    if total == 0.0 {
        return 0.0;
    }
    let h = shannon_entropy(&sample);
    #[allow(clippy::cast_precision_loss)]
    let l = n as f64;
    100.0 * h.exp() / l
}

/// Choose the number of sample positions: halve `n` while `n >= mrna_len`.
///
/// From `tin.py`: "if this number is larger than the length of mRNA (L),
/// it will be halved until it's smaller than L."
#[must_use]
pub fn effective_sample_size(mut n: usize, mrna_len: usize) -> usize {
    while n >= mrna_len {
        n /= 2;
        if n == 0 {
            return 0;
        }
    }
    n
}

/// Per-transcript accumulator used during the linear BAM sweep.
struct TxAccum {
    positions: Vec<i64>,         // sorted, 1-based sampled genomic positions
    coverage: Vec<f64>,          // coverage[i] corresponds to positions[i]
    unique_starts: HashSet<i64>, // 0-based read starts within [tx_start, tx_end)
}

#[allow(clippy::too_many_lines)]
/// Compute TIN for all transcripts in a BED12 from a sorted, indexed BAM.
///
/// Uses a single sequential BAM scan rather than per-transcript indexed seeks,
/// avoiding repeated BGZF-decompress cycles when transcripts are dense.
pub fn compute_tin(
    bam_path: &Path,
    bed_path: &Path,
    min_cov: u64,
    sample_size: usize,
) -> Result<Vec<TinOutput>> {
    let transcripts = parse_bed12(bed_path)?;

    // Sort by chrom then tx_start for sweep-line ordering; keep original index for output order.
    let mut tx_by_idx: Vec<(usize, &_)> = transcripts.iter().enumerate().collect();
    tx_by_idx.sort_unstable_by(|a, b| {
        a.1.chrom
            .cmp(&b.1.chrom)
            .then(a.1.tx_start.cmp(&b.1.tx_start))
    });

    let positions_vec: Vec<Vec<i64>> = transcripts
        .iter()
        .map(|tx| pick_positions(tx, sample_size))
        .collect();

    let mut accums: Vec<TxAccum> = positions_vec
        .iter()
        .map(|pos| TxAccum {
            coverage: vec![0.0; pos.len()],
            positions: pos.clone(),
            unique_starts: HashSet::new(),
        })
        .collect();

    let file = File::open(bam_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bam_path.display())))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    // chrom → sorted (tx_start, tx_end, original_index) for sweep-line overlap queries
    let mut chrom_txs: std::collections::HashMap<String, Vec<(i64, i64, usize)>> =
        std::collections::HashMap::new();
    for (orig_idx, tx) in transcripts.iter().enumerate() {
        chrom_txs
            .entry(tx.chrom.clone())
            .or_default()
            .push((tx.tx_start, tx.tx_end, orig_idx));
    }
    for v in chrom_txs.values_mut() {
        v.sort_unstable_by_key(|&(s, _, _)| s);
    }

    let mut record = bam::Record::default();
    loop {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(RsomicsError::Io(e)),
        }

        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_qc_fail() {
            continue;
        }

        let start_1based: i64 = match record.alignment_start() {
            Some(Ok(pos)) => {
                #[allow(clippy::cast_possible_wrap)]
                let v = usize::from(pos) as i64;
                v
            }
            _ => continue,
        };
        let start_0based = start_1based - 1;

        // Determine reference name from the record's reference sequence ID.
        let Some(Ok(ref_id)) = record.reference_sequence_id() else {
            continue;
        };
        let chrom: &str = match header.reference_sequences().get_index(ref_id) {
            Some((name, _)) => match std::str::from_utf8(name.as_ref()) {
                Ok(s) => s,
                Err(_) => continue,
            },
            None => continue,
        };

        let Some(tx_list) = chrom_txs.get(chrom) else {
            continue;
        };

        // Compute read end from CIGAR (reference-consuming length).
        let cigar = record.cigar();
        let read_ref_len: i64 = cigar
            .iter()
            .filter_map(std::result::Result::ok)
            .map(|op| match op.kind() {
                Kind::Match
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
                | Kind::Deletion
                | Kind::Skip => {
                    #[allow(clippy::cast_possible_wrap)]
                    {
                        op.len() as i64
                    }
                }
                _ => 0,
            })
            .sum();
        // Find transcripts overlapping [start_0based, start_0based + read_ref_len) using
        // binary search into the chromosome-sorted transcript list.
        // A transcript overlaps the read if tx_start < search_end AND tx_end > start_0based.
        let search_end = start_0based + read_ref_len; // exclusive 0-based read end
        let first_candidate = tx_list.partition_point(|&(s, _, _)| s < start_0based);
        let upper = tx_list.partition_point(|&(s, _, _)| s < search_end);

        // Determine lower bound: look back for transcripts that started before start_0based
        // but whose tx_end still overlaps (long transcripts).
        let lower = if first_candidate > 0 {
            let mut lo = first_candidate;
            while lo > 0 && tx_list[lo - 1].1 > start_0based {
                lo -= 1;
            }
            lo
        } else {
            0
        };

        for &(tx_start, tx_end, orig_idx) in &tx_list[lower..upper] {
            if tx_end <= start_0based {
                continue;
            }
            let accum = &mut accums[orig_idx];
            let positions = &accum.positions;

            if positions.is_empty() {
                continue;
            }

            // Track unique read starts within [tx_start, tx_end) for min_cov gate.
            if start_0based >= tx_start && start_0based < tx_end {
                accum.unique_starts.insert(start_0based);
            }

            // Walk CIGAR match intervals, accumulate coverage at sampled positions.
            let last_pos = *positions.last().unwrap();
            #[allow(clippy::cast_sign_loss)]
            let mut ref_cursor = start_1based as usize;
            for op_result in cigar.iter() {
                let Ok(op) = op_result else { break };
                let len = op.len();
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        #[allow(clippy::cast_possible_wrap)]
                        let block_start = ref_cursor as i64;
                        #[allow(clippy::cast_possible_wrap)]
                        let block_end = block_start + len as i64;
                        if block_start > last_pos {
                            break;
                        }
                        let lo = positions.partition_point(|&p| p < block_start);
                        let hi = positions.partition_point(|&p| p < block_end);
                        for idx in lo..hi {
                            accum.coverage[idx] += 1.0;
                        }
                        ref_cursor += len;
                    }
                    Kind::Deletion | Kind::Skip => {
                        ref_cursor += len;
                    }
                    Kind::SoftClip | Kind::HardClip | Kind::Pad | Kind::Insertion => {}
                }
            }
        }
    }

    let mut results: Vec<Option<TinOutput>> = (0..transcripts.len()).map(|_| None).collect();
    for (orig_idx, tx) in transcripts.iter().enumerate() {
        let accum = &accums[orig_idx];
        let positions = &accum.positions;

        let tin = if positions.is_empty() || accum.unique_starts.len() as u64 <= min_cov {
            0.0
        } else {
            let h = shannon_entropy(&accum.coverage);
            if h == 0.0 {
                0.0
            } else {
                #[allow(clippy::cast_precision_loss)]
                let l_f = positions.len() as f64;
                100.0 * h.exp() / l_f
            }
        };

        eprint!(" {} transcripts finished", orig_idx + 1);
        results[orig_idx] = Some(TinOutput {
            gene_id: tx.gene_id.clone(),
            chrom: tx.chrom.clone(),
            tx_start: tx.tx_start,
            tx_end: tx.tx_end,
            tin,
        });
    }

    eprintln!();
    Ok(results.into_iter().map(|r| r.unwrap()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tin_uniform_100_positions() {
        // All 100 positions have coverage 10 → p_i = 1/100 → H = ln(100) → TIN = 100
        let cov: Vec<u32> = vec![10; 100];
        let tin = compute_tin_score(&cov, 100);
        assert!(
            (tin - 100.0).abs() < 1e-10,
            "expected TIN=100 for uniform, got {tin}"
        );
    }

    #[test]
    fn tin_single_spike() {
        // All coverage at one position → H = 0 → TIN = exp(0)/n * 100 = 1.0
        let mut cov = vec![0u32; 100];
        cov[50] = 20;
        let tin = compute_tin_score(&cov, 100);
        assert!(
            (tin - 1.0).abs() < 1e-10,
            "expected TIN≈1 for spike, got {tin}"
        );
    }

    #[test]
    fn tin_zero_coverage() {
        let cov = vec![0u32; 100];
        let tin = compute_tin_score(&cov, 100);
        assert!(tin == 0.0, "expected 0.0, got {tin}");
    }

    #[test]
    fn effective_sample_size_no_halving() {
        assert_eq!(effective_sample_size(100, 1000), 100);
    }

    #[test]
    fn effective_sample_size_halving_once() {
        // n=100, mrna_len=100: 100 >= 100 → halve to 50 → 50 < 100
        assert_eq!(effective_sample_size(100, 100), 50);
    }

    #[test]
    fn effective_sample_size_halving_twice() {
        // n=100, mrna_len=50: 100≥50 → 50≥50 → 25 < 50
        assert_eq!(effective_sample_size(100, 50), 25);
    }
}
