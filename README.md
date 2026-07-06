# rsomics-tin

Compute the **Transcript Integrity Number (TIN)** for each transcript in a BED12
gene model from a sorted, indexed BAM file. TIN is an entropy-based score
(0–100) measuring read-coverage uniformity across the mRNA body — a
transcript-level analogue of RIN that is more sensitive for low-quality RNA.

A Rust port of RSeQC's `tin.py`, byte-exact with the upstream output and
substantially faster (single linear BAM scan instead of per-transcript indexed
queries).

## Install

```sh
cargo install rsomics-tin
```

## Usage

```sh
rsomics-tin -i sample.bam -r genes.bed12
```

Writes `<bam_basename>.tin.xls` (per-transcript: geneID, chrom, tx_start,
tx_end, TIN) and `<bam_basename>.summary.txt` (sample-level mean/median/stdev
TIN).

Options:

| flag | meaning | default |
|---|---|---|
| `-i, --input` | sorted+indexed BAM(s): comma-separated, a directory, or a text file of paths | required |
| `-r, --refgene` | reference gene model (BED12) | required |
| `-c, --minCov` | minimum reads mapped to a transcript (below → TIN=0) | 10 |
| `-n, --sample-size` | equally-spaced mRNA positions sampled per transcript | 100 |

With `--json`, a single JSON envelope is written to stdout; the file outputs
are still written. The `result` object carries one entry per processed BAM:

```json
{
  "result": {
    "samples": [
      { "bam_file": "sample.bam", "mean": 82.4, "median": 85.1, "stdev": 12.3, "n_transcripts": 1200 }
    ]
  }
}
```

## Method

For each transcript, per-base coverage is computed across the concatenated mRNA
(introns skipped). With relative coverage `p_i = C_i / Σ C_j`, Shannon entropy
`H = −Σ p_i ln p_i` (over covered positions), and `TIN = 100 × exp(H) / n`,
where `n` is the number of sampled positions. Transcripts with fewer than
`--minCov` unique read-start positions are assigned TIN = 0.

## Origin

This crate is an independent Rust reimplementation of `tin.py` (RSeQC) based on:

- The published method: Wang, Wang, Kalari, et al. "Measure transcript integrity
  using RNA-seq data." *BMC Bioinformatics* 2016;17:58.
  doi:[10.1186/s12859-016-0922-z](https://doi.org/10.1186/s12859-016-0922-z)
- The public BED12 + SAM/BAM format specifications
- Black-box behavior testing against the `tin.py` 5.0.4 binary

No source code from the GPL-licensed RSeQC upstream was used as reference during
implementation.

License: MIT OR Apache-2.0.
Upstream credit: [RSeQC](https://rseqc.sourceforge.net/) (GNU GPL).
