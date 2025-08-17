# fastats

CLI to generate statistics from [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files:

- Generates **[BED](https://en.wikipedia.org/wiki/BED_(file_format)) files** for **non-masked** (`A|C|G|T`), **soft-masked** (`a|c|g|t`), and **hard-masked regions** (`n|N`), per sequence.
- Stores **overall statistics** (GC content, ratios of masked bases) to `stdout` and **JSON**.

## Details

```text
CLI to generate FASTA file statistics (masking, GC content, etc.).

Usage: fastats [OPTIONS] <FASTA_FILE>

Arguments:
  <FASTA_FILE>  

Options:
  -o, --output-dir <OUTPUT_DIR>
          The output directory for the BED and summary files. [default: .]
  -q, --quiet
          Do not print results on stdout.
      --ignore-iupac
          Enable this to avoid failing when encountering a sequence character that is not in ('A', 'C', 'T', 'G', 'N', 'a', 'c', 't', 'g', 'n').
      --no-bed-output
          Do not store masking regions into BED files.
      --match-regex <SEQUENCE_MATCH_REGEX>
          Regular expression to focus the analysis on sequences matching a specific regular expression. [default: .*]
  -h, --help
          Print help
  -V, --version
          Print version
```

## Sample output

### Bed file per sequence

For each sequence, BED files that report the non-masked, soft-masked, and hard-masked regions are define. 
They use the simple three-column BED format.
Sample output:

```text
chr9 0 10000
chr9 40529470 40529480
...
```

### Summary statistics

Summary statistics are printed out to `stdout` and into a `summary.json` file.
Sample output:

```json
[
  {
    "sequence_name": "sample_sequence",
    "non_masked_bases": 304,
    "soft_masked_bases": 36936,
    "hard_masked_bases": 0,
    "non_masked_ratio": 0.00816326530612245,
    "soft_masked_ratio": 0.9918367346938776,
    "hard_masked_ratio": 0.0,
    "gc_content": 0.4293233082706767,
    "other_iupac_bases": 0,
    "sequence_length": 37240,
    "checksum_sha256": "4b2a8b27c0f83f7d72600e33af490149d027b3e6c1e81987730a7561cde563a8"
  },
  ...
]
```

## Usage examples

### Get sorted list of sequence names

```shell
fastats hg38.fasta | jq '.[].sequence_name'
```

### Calculate the overall sequence length

```shell
fastats hg38.fasta | jq '.[].sequence_length' | paste -sd+ | bc
```

### Print stats for all sequences without a `_` in the name

```shell
fastats hg38.fasta --match-regex "[^_]*"
```

## Notes

- Note that the base `n` is _not_ considered soft-masked (so the sum of all non-masked, soft-masked, hard-masked, and non-supported IUPAC code bases equals the overall sequence length).

- Ambiguous [IUPAC codes](https://genome.ucsc.edu/goldenPath/help/iupac.html) (i.e., any code except `N`, `A`, `C`, `G`, or `T`) are not supported. To ingest sequences containing such IUPAC codes, use `--ignore-iupac`.
