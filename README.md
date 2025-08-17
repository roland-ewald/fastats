# fastats

Generate statistics from FASTA files.


## Usage examples

```
```


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
fastats hg38.fasta  --match-regex "[^_]*"
```