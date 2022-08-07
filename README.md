# RivIndel
IN CONSTRUCTION!

## Motivation
Complex indels (cxIndels) often introduce multiple deletions, insertions and mismatches during the alignment phase.
Complex indels are prone to be found on highly unstable genomes such as Lung Adenocarcinoma (LUAD).
A correct interpretation in terms of net base change is crucial: e.g. exon 19 EGFR inframe deletions are actionable, but not frameshift deletions.

These variants are typically called as multiple individual variants by some popular somatic callers (Mutect2, Lancet, etc).
Since these individual variants are closely located, there is a high chance that they come from a single compound variant.

RivIndel detects complex indels from targeted sequencing (Illumina).
RivIndel scans every targeted region (supplied as a BED file) in a similar fashion than indelSeek.
However, RivIndel validates candidate cxIdels after a successful assembly of candidate reads (plus neighbouring soft-clipped reads).
This, in turn, leads to a better allele frequency estimation.

## Method overview
![Slide](img/rivindel_schema2.png)

## Install RivIndel locally

1. Clone RivIndel repository:
```
git clone https://github.com/bdolmo/RivIndel.git
cd RivIndel/
```

2. Install via pip
```
pip3 install -e .
```

## Usage
```
python3 RivIndel.py --bam <input.bam> --bed <test.bed> --vcf <out.vcf>
```
