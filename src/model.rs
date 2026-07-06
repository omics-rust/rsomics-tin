/// One row in the per-transcript output (`.tin.xls`).
#[derive(Debug, Clone)]
pub struct TinOutput {
    pub gene_id: String,
    pub chrom: String,
    pub tx_start: i64,
    pub tx_end: i64,
    /// TIN score in [0, 100]. 0.0 means below minCov or zero coverage.
    pub tin: f64,
}

/// Sample-level summary statistics emitted to `.summary.txt`.
#[derive(Debug, serde::Serialize)]
pub struct TinSummary {
    pub mean: f64,
    pub median: f64,
    pub stdev: f64,
}

/// Per-BAM entry in the `--json` report.
#[derive(Debug, serde::Serialize)]
pub struct SampleReport {
    pub bam_file: String,
    pub mean: f64,
    pub median: f64,
    pub stdev: f64,
    pub n_transcripts: usize,
}

/// The `--json` result object: one entry per processed BAM.
#[derive(Debug, serde::Serialize)]
pub struct TinReport {
    pub samples: Vec<SampleReport>,
}

/// One transcript parsed from BED12.
#[derive(Debug, Clone)]
pub(crate) struct Transcript {
    pub(crate) gene_id: String,
    pub(crate) chrom: String,
    /// 0-based, half-open `[tx_start, tx_end)`
    pub(crate) tx_start: i64,
    pub(crate) tx_end: i64,
    /// Exon blocks as 0-based half-open `[exon_start, exon_end)` in genomic coordinates.
    pub(crate) exons: Vec<(i64, i64)>,
}

impl Transcript {
    pub(crate) fn mrna_len(&self) -> i64 {
        self.exons.iter().map(|(s, e)| e - s).sum()
    }
}
