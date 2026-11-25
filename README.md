# subclonal

Allele frequencies from bacterial pileups, without genotype calls.

## The problem

A genotype caller answers one question: *which base is this sample?* For a
clonal isolate that is the right question. For a bacterial population sampled
as colonies, it is not.

A site at 60% A and 40% G is not a noisy A. It is two lineages coexisting in
that sample, and the ratio is the measurement. Every genotype caller will
collapse it to a single call, and the thing you were trying to see is gone
before you look at the output.

`subclonal` reports the frequency at each variable site in each sample and
stops there. Deciding what a frequency means is your job, not the tool's.

## What you need

- `samtools`
- Python 3.7+
- One sorted, indexed BAM per sample, all mapped to the same reference

## Install

```bash
git clone https://github.com/clairezyy603/subclonal.git
cd subclonal
chmod +x subclonal.py
```

## Use

Write a manifest with `sample` and `bam` columns:

```
sample      bam
colony_01   bams/colony_01.bam
colony_02   bams/colony_02.bam
colony_03   bams/colony_03.bam
```

```bash
./subclonal.py freq \
    --reference reference.fna \
    --manifest samples.tsv \
    --min-depth 10 \
    -o frequencies.tsv
```

One row per variable site, two columns per sample — frequency and depth:

| contig | pos | ref | alt | n_samples_variant | colony_01_freq | colony_01_depth | … |
| --- | --- | --- | --- | --- | --- | --- | --- |
| chr | 503494 | G | A | 7 | 0.982 | 46 | … |
| chr | 2950850 | C | T | 1 | 0.511 | 43 | … |

A site is written when samples disagree. Sites where every sample carries the
alt allele are invariant *between* samples and carry no information about
population structure, so they are dropped unless you pass `--keep-invariant`.

### Merging adjacent sites

```bash
./subclonal.py events --freq frequencies.tsv --max-gap 200 -o events.tsv
```

Sites within `--max-gap` of each other that are carried by exactly the same
set of samples are almost certainly one event, not several — a short indel, a
gene conversion tract, a recombination block. Counting them separately
inflates anything you compute downstream.

In the dataset this was written for, 702 variable sites merged into 605
events, and a functional enrichment that looked like 3.7-fold became
2.5-fold. The 3.7 was four adjacent positions in one gene being counted four
times. **Run this before any enrichment or burden test.**

Check that your result is not an artefact of the window:

```bash
for g in 200 500 1000 2000; do
  ./subclonal.py events --freq frequencies.tsv --max-gap $g -o /dev/null
done
```

## Options worth knowing

| flag | default | what it does |
| --- | --- | --- |
| `--min-depth` | 10 | below this a sample is `NA` at that site rather than 0% |
| `--min-freq` | 0.10 | the frequency at which a sample counts as carrying the alt allele |
| `--mapq` | 20 | minimum mapping quality — see the paralogue note below |
| `--region` | — | `contig:start-end`, for testing on one locus before running genome-wide |
| `--keep-invariant` | off | also write sites where every sample agrees |

## Validation

`validate.py` builds a synthetic case with known answers: two lineages
differing at planted positions, mixed at specified ratios, reads simulated
from the mixture, and the recovered frequencies compared against what went in.
It includes three adjacent positions carried by the same samples, so `events`
is tested too.

```bash
python3 validate.py
```

This is not decoration. The original version of this parser had an off-by-one
in the field offset — it read the wrong column of the pileup and returned zero
differences for every pair of samples. On real data that looks like a
biological result: *these colonies are identical*. It was only visible against
a dataset where the answer was known in advance. Run the validation after any
change.

## Limits

**Frequencies are not genotypes, and this tool will not make that decision for
you.** 40% could be two lineages in a mixed colony, or contamination, or
mismapping. The tool reports the number; separating those is a question about
your experiment, not about the data.

**Paralogues will produce false frequencies.** If a gene has close homologues
elsewhere in the genome, reads land on the wrong copy and you get an apparent
intermediate frequency with no biology behind it. `--mapq 20` removes the worst
of it, but it is not enough on its own. Before trusting a signal in a gene
family, check coverage against the genome mean and check the fraction of reads
at MAPQ 0 in that window. A region that looks variable and has 18% MAPQ-0 reads
is telling you about the aligner.

**It is only as good as the reference.** Every frequency is relative to one
assembly. If the reference happens to carry a minority allele at some position,
every sample will look variant there — which says nothing about the population.
In the dataset this was written for, two such sites appeared in 36 and 33 of 44
colonies, and both turned out to be the reference being the odd one out.

**Indels are skipped, not called.** The parser steps over them correctly so
they do not corrupt the base counts, but they are not reported. If indels
matter for your question, this is the wrong tool.

**Low depth is reported, not imputed.** A sample below `--min-depth` at a site
gets `NA`, not 0%. Treat those as missing.

## Licence

MIT
