# fastats

CLI to generate statistics from FASTA files:

- Generates **BED files** for **non-masked**, **soft-masked**, and **hard-masked regions**, per sequence.
- Stores **overall statistics** (GC content, ratios of masked bases) to `stdout` and **JSON**.

## Details

```text
Provides a CLI to generate statistics for FASTA files.

Usage: fastats [OPTIONS] <FASTA_FILE>

Arguments:
  <FASTA_FILE>  

Options:
  -o, --output-dir <OUTPUT_DIR>
          The output directory for the BED and summary files. [default: .]
  -q, --quiet
          Do not print results on stdout.
      --no-bed-output
          Do not store masking regions into BED files.
      --match-regex <SEQUENCE_MATCH_REGEX>
          Regular expression to focus the analysis on sequences matching a specific regular expression. [default: .*]
  -h, --help
          Print help
  -V, --version
          Print version
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
