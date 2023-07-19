# VEME 2023 NGS *De novo* Assembly Tutorial

## Taylor K. Paisie
### `https://taylorpaisie.github.io/VEME_2023_NGS_Denovo_assembly/`

### 1. Introduction to sequence assembly
#### What does "sequence assembly" mean?
#### Assembly is a “catch-all” term used to describe methods where we combine shorter individual measurements called reads into longer contiguous sequences typically called contigs
#### Because the sequencing process works by breaking the original DNA into smaller fragments, the assembly process is conceptually similar to putting together an image puzzle from its many pieces
#### The software that performs the assembly is called the assembler.
#### We will learn how to *de novo* assemble reads obtained with the Illumina sequencing platform using [SPAdes](http://cab.spbu.ru/software/spades/), an assembly toolkit containing various assembly pipelines

#### *De novo* assembly usually includes the following steps:  
1) Improving of the reads quality (remove adapters, trim reads, etc.)
2) De novo assembly of the overlapping reads into contigs
3) Joining contigs into scaffolds
4) Comparison with other known genomes
5) Filling the gaps
6) Verification of the assembled genome
7) Annotation of the assembled genome

#### 
#### 
### 2. 
<figure>
    <img src="variant_calling_steps.png" width="230" height="300">
    <figcaption>Variant Calling Workflow</figcaption>
</figure>