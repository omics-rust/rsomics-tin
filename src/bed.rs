use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::model::Transcript;

/// Parse BED12 file into a list of transcripts.
pub(crate) fn parse_bed12(path: &Path) -> Result<Vec<Transcript>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| RsomicsError::Io(std::io::Error::other(format!("reading BED12: {e}"))))?;

    let mut transcripts = Vec::new();

    for line in content.lines() {
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        let chrom = fields[0].to_string();
        let Ok(tx_start) = fields[1].parse::<i64>() else {
            continue;
        };
        let Ok(tx_end) = fields[2].parse::<i64>() else {
            continue;
        };
        let gene_id = fields[3].to_string();

        let Ok(block_count) = fields[9].parse::<usize>() else {
            continue;
        };

        let block_sizes: Vec<i64> = fields[10]
            .trim_end_matches(',')
            .split(',')
            .filter(|s| !s.is_empty())
            .filter_map(|s| s.parse().ok())
            .collect();
        let block_starts: Vec<i64> = fields[11]
            .trim_end_matches(',')
            .split(',')
            .filter(|s| !s.is_empty())
            .filter_map(|s| s.parse().ok())
            .collect();

        if block_sizes.len() != block_count || block_starts.len() != block_count {
            continue;
        }

        let exons: Vec<(i64, i64)> = block_starts
            .iter()
            .zip(block_sizes.iter())
            .map(|(&rel_start, &size)| {
                let abs_start = tx_start + rel_start;
                (abs_start, abs_start + size)
            })
            .collect();

        if exons.is_empty() {
            continue;
        }

        transcripts.push(Transcript {
            gene_id,
            chrom,
            tx_start,
            tx_end,
            exons,
        });
    }

    Ok(transcripts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_bed12_single_exon() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "chr1\t0\t1000\ttx1\t0\t+\t0\t1000\t0\t1\t1000,\t0,").unwrap();
        let txs = parse_bed12(f.path()).unwrap();
        assert_eq!(txs.len(), 1);
        assert_eq!(txs[0].mrna_len(), 1000);
        assert_eq!(txs[0].exons, vec![(0, 1000)]);
    }

    #[test]
    fn parse_bed12_multi_exon() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        // 3 exons: offsets 0,200,300 sizes 100,150,200 → absolute: [100,200],[300,450],[400,600]
        writeln!(
            f,
            "chr1\t100\t600\ttx2\t0\t+\t100\t600\t0\t3\t100,150,200,\t0,200,300,"
        )
        .unwrap();
        let txs = parse_bed12(f.path()).unwrap();
        assert_eq!(txs[0].mrna_len(), 450);
        assert_eq!(txs[0].exons, vec![(100, 200), (300, 450), (400, 600)]);
    }
}
