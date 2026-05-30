use crate::model::{TinOutput, TinSummary};

/// Compute sample-level summary statistics over non-zero TIN scores.
///
/// Mirrors `tin.py`: only transcripts with TIN > 0 contribute. Standard deviation
/// uses the population formula (divide by N, not N-1).
#[must_use]
pub fn summarise(outputs: &[TinOutput]) -> TinSummary {
    let nonzero: Vec<f64> = outputs
        .iter()
        .filter(|r| r.tin > 0.0)
        .map(|r| r.tin)
        .collect();

    if nonzero.is_empty() {
        return TinSummary {
            mean: 0.0,
            median: 0.0,
            stdev: 0.0,
        };
    }

    #[allow(clippy::cast_precision_loss)]
    let n = nonzero.len() as f64;
    let mean = nonzero.iter().sum::<f64>() / n;

    let mut sorted = nonzero.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if sorted.len() % 2 == 1 {
        sorted[sorted.len() / 2]
    } else {
        let mid = sorted.len() / 2;
        sorted[mid - 1] / 2.0 + sorted[mid] / 2.0
    };

    // Population standard deviation (divides by N, matching numpy default)
    let variance = nonzero.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
    let stdev = variance.sqrt();

    TinSummary {
        mean,
        median,
        stdev,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::TinOutput;

    #[test]
    fn summarise_nonzero_only() {
        let outputs = vec![
            TinOutput {
                gene_id: "a".into(),
                chrom: "chr1".into(),
                tx_start: 0,
                tx_end: 100,
                tin: 90.0,
            },
            TinOutput {
                gene_id: "b".into(),
                chrom: "chr1".into(),
                tx_start: 0,
                tx_end: 100,
                tin: 0.0,
            },
            TinOutput {
                gene_id: "c".into(),
                chrom: "chr1".into(),
                tx_start: 0,
                tx_end: 100,
                tin: 80.0,
            },
        ];
        let s = summarise(&outputs);
        assert!((s.mean - 85.0).abs() < 1e-10, "mean = {}", s.mean);
        assert!((s.median - 85.0).abs() < 1e-10, "median = {}", s.median);
        // pop stdev of [90, 80]: sqrt(((90-85)^2 + (80-85)^2) / 2) = sqrt(25) = 5
        assert!((s.stdev - 5.0).abs() < 1e-10, "stdev = {}", s.stdev);
    }
}
