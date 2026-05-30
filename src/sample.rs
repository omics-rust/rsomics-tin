use std::collections::HashSet;

use crate::model::Transcript;

/// Build the sampled genomic positions for a transcript, matching `tin.py` exactly.
///
/// Positions are 1-based genomic coordinates. Exon boundary positions are always
/// included in addition to the equally-spaced sample. When `mrna_len <= sample_size`,
/// every mRNA base is returned (all-bases mode).
///
/// Returns the sorted, deduplicated position list. Avoids materialising the full
/// per-base mRNA index list by walking exons arithmetically.
pub(crate) fn pick_positions(tx: &Transcript, sample_size: usize) -> Vec<i64> {
    #[allow(clippy::cast_sign_loss)]
    let mrna_len = tx.mrna_len() as usize;
    if mrna_len == 0 {
        return Vec::new();
    }

    // Collect exon boundary positions (1-based).
    let mut exon_bounds: Vec<i64> = Vec::with_capacity(tx.exons.len() * 2);
    for &(ex_start, ex_end) in &tx.exons {
        exon_bounds.push(ex_start + 1);
        exon_bounds.push(ex_end);
    }

    // Generate the equally-spaced sample without building the full per-base list.
    // tin.py: step = int(mRNA_size / sample_size); indx = range(0, len(gene_all_base), step)
    // gene_all_base[i] = genomic position of mRNA base at index i.
    // We compute gene_all_base[i] by walking exons: find which exon contains mRNA index i.
    let chose_bases: Vec<i64> = if mrna_len <= sample_size {
        // All-bases mode: every mRNA position in genomic order.
        let mut v: Vec<i64> = Vec::with_capacity(mrna_len);
        for &(ex_start, ex_end) in &tx.exons {
            for gpos in (ex_start + 1)..=ex_end {
                v.push(gpos);
            }
        }
        v
    } else {
        let step = mrna_len / sample_size;
        // Walk exons to map mRNA index i → genomic 1-based position.
        // Collect indices 0, step, 2*step, ... < mrna_len.
        let mut sampled: Vec<i64> = Vec::with_capacity(mrna_len / step + 1);
        let mut mrna_offset: usize = 0; // mRNA index of the start of the current exon
        let mut next_idx: usize = 0; // next mRNA sample index to collect

        for &(ex_start, ex_end) in &tx.exons {
            #[allow(clippy::cast_sign_loss)]
            let exon_len = (ex_end - ex_start) as usize;
            // Indices for this exon: [mrna_offset, mrna_offset + exon_len)
            // Sample indices in this range: next_idx, next_idx+step, ...
            while next_idx < mrna_offset + exon_len {
                let local = next_idx - mrna_offset; // offset within this exon
                #[allow(clippy::cast_possible_wrap)]
                let gpos = ex_start + 1 + local as i64; // 1-based genomic
                sampled.push(gpos);
                next_idx += step;
            }
            mrna_offset += exon_len;
        }
        sampled
    };

    // Merge exon bounds + sampled bases, deduplicate (first-seen wins), then sort.
    // Mirrors tin.py: sorted(uniqify(exon_bounds + chose_bases)).
    let mut seen: HashSet<i64> = HashSet::with_capacity(exon_bounds.len() + chose_bases.len());
    let mut merged: Vec<i64> = Vec::with_capacity(exon_bounds.len() + chose_bases.len());
    for pos in exon_bounds.iter().chain(chose_bases.iter()) {
        if seen.insert(*pos) {
            merged.push(*pos);
        }
    }
    merged.sort_unstable();
    merged
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Transcript;

    #[test]
    fn pick_positions_single_exon() {
        // tx: chr1:0-1000, one exon [0,1000), mrna_len=1000, sample_size=100
        // step=10, gene_all_base=[1..=1000], indx=[0,10,...,990]
        // chose_bases=[1,11,...,991], exon_bounds=[1,1000]
        // uniqify([1,1000]+[1,11,...,991])=[1,1000,11,...,991] → sorted=[1,11,...,991,1000]
        let tx = Transcript {
            gene_id: "tx".into(),
            chrom: "chr1".into(),
            tx_start: 0,
            tx_end: 1000,
            exons: vec![(0, 1000)],
        };
        let positions = pick_positions(&tx, 100);
        assert_eq!(positions.len(), 101, "expected 101 positions");
        assert_eq!(positions[0], 1, "first position should be 1 (1-based)");
        assert_eq!(
            *positions.last().unwrap(),
            1000,
            "last position should be 1000"
        );
    }
}
