# ssDripPipeline
This pipeline is used to analyze ssDRIP-seq data. The following operations can be automatically performed:

**BaseAnalysis**
- Alignment of reads to reference sequences
- Duplicates removing
- Strand splitting
- Peak calling
- Normalized bam to bigwig file

**DeseqAnalysis**
- DESeq2 for peaks(merge the peaks of all samples)

**DownstreamAnalysis**
- Mfuzz cluster(peak with qvalue<=0.01)
- Correlation of samples
- Motif for peaks
- Peaks length distribution
- GCskew and ATskew
- Sense and Antisense metaplot
- Peaks content distribution(the proportion of peaks in TSS,TTS and gene body region)
## Installation
ssDripPipeline can be installed using the [conda](http://conda.pydata.org/docs/intro.html) package manager [Bioconda](https://bioconda.github.io/). To install ssDripPipeline using Bioconda, download and install Anaconda Python 3.7, following the instructions at: https://www.anaconda.com/distribution/.

To install ssDripPipeline into the current conda environment:
```bash
conda install -c bioconda ssDripPipeline
```
Alternately, to create a new environment named `ssDripPipeline_env` with ssDripPipeline:
```bash
conda create -n ssDripPipeline_env -c bioconda ssDripPipeline python=3.7
```
Activate your conda environment:
```bash
conda activate ssDripPipeline_env
```
Verify that ssDripPipeline is installed using the command:
```
ssDRIPSeqAnalysis.py

Usage:

ssDRIPSeqAnalysis.py <DripConfig.json> <BaseAnalysis|DeseqAnalysis|DownstreamAnalysis|AllPip>

```
## Commands
**ssDRIPSeqAnalysis.py** contains four subcommands:

**BaseAnalysis**
```bash
ssDRIPSeqAnalysis.py DripConfig.json BaseAnalysis
```
**DeseqAnalysis**

The results of **BaseAnalysis** is required
```bash
ssDRIPSeqAnalysis.py DripConfig.json DeseqAnalysis
```
**DownstreamAnalysis**

The results of **DeseqAnalysis** and **BaseAnalysis** are required
```bash
ssDRIPSeqAnalysis.py DripConfig.json DownstreamAnalysis
```
**AllPip**

Execute BaseAnalysis, DeseqAnalysis, DownstreamAnalysis in turn
```bash
ssDRIPSeqAnalysis.py DripConfig.json AllPip
```
## Inputs
ssDRIPSeqAnalysis.py requires a configuration file in json format as input.\
If you copy the following json configuration file, please delete the comments, because json does not support comments.
```
{
        //The name of your project.
        "ProjectName":"test",
        //A file that contains the sample name and its corresponding data.
        "Target":"target",
        //Path of sequencing data, it must be an absolute path.
        "SeqDataPath":"/data/yeat/cleandata/",
        //bowtie2-build outputs, it must be an absolute path. GenomeIndex indicates the prefix name of bowtie2-build results.
        "GeomeIndex":"/data/yeast/genome/GenomeIndex",
        //Effective genome size, it used for peak calling, normalization, and peaks proportion calculation.
        "GenomeSize":"12345",
        //Threads of some third-party software.
        "Thread":"10",
        //Filter out chromosomes that do not need to be analyzed, such as chloroplasts, mitochondria and spike-in chromosomes.
        "FilterChromFile":"/data/yeast/Analysis/filter.txt",
        //Size of the bins, in base, for the bamCoverage and plotCorrelation in the deeptools.
        "BinSize":"5",
        //The number of times a set random regions are selected as a background for motif analysis and normalization. The number and length distribution of a set of random regions are consistent with the number and length distribution of peak regions.
        "RepeatNum":"10",
        //The file contains two columns. The first column is the chromosome name, and the second column is the chromosome length. The chromosome name must be consistent with the name in the genome fasta file. The filtered chromosomes should not be included in this file. In addition, the absolute path cannot be missing.
        "ChromSize":"/data/yeast/Analysis/chrom.size",
        //Fasta file of genome. The absolute path cannot be missing.
        "GenomeFastaFile":"/data/yeast/genome/Genome.fasta",
        //If there is no control sample for DESeq2 analysis, you can write a random name, such as hehe.
        "ControlSample":"hehe",
        //Bed file of gene. The absolute path cannot be missing.
        "GeneBeb":"/data/yeast/genome/gene.bed",
        //
        "MetaplotExtend":"1000",
        //
        "ContentExtend":"150",
        //The name of the chromosome in the fasta file of the Escherichia coli genome. If you do not do spike-in normalization, please ignore this.
        "SpikeNameList":["Ecoli"],
        //Debug
        "Debug":"False",
        //Divide each input interval to fixed-sized windows in chromosome.
        "Win":"100",
        //How many base pairs to step before creating a new window. This parameter is used together with the win parameter to calculate ATskew and GCskew.
        "Step":"50"
}
```
## Outputs
```

├── Callus1
│   ├── Callus1_align.info
│   ├── Callus1_fwd.bam
│   ├── Callus1_fwd.bam.bai
│   ├── Callus1_fwd.bw
│   ├── Callus1_fwd.gz
│   ├── Callus1_fwd_nucl
```
## Third-party software
The following software is used by this pipeline. When installing ssDripPipeline, these software will be installed automatically, and you do not need to install them yourself.
- bowtie2
- picard
- samtools
- bedtools
- deeptools
- MACS2
- Homer
- DESeq2
- mfuzz
- bedGraphToBigWig
## Something
1. Only supports paired-end sequencing data.
2. ssDripPipeline does not include quality control, adapter cutting and tail trimming.
3. [Effective genome size](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html)
4. Step by step protocols [wiki](https://github.com/PEHGP/ssDripPipeline/wiki).
5. [Sam/Bam format](https://samtools.github.io/hts-specs/SAMv1.pdf)
6. [Bigwig format](https://genomebrowser.wustl.edu/goldenPath/help/bigWig.html)
7. [Bed format](https://genome-asia.ucsc.edu/FAQ/FAQformat.html#format1)
8. [Metaplot](https://deeptools.readthedocs.io/en/latest/content/tools/plotProfile.html)
9. A detailed explanation of the ssDripPipeline results can be found in these two articles([paper1](https://www.nature.com/articles/s41477-017-0004-x),[paper2](https://academic.oup.com/plcell/article/32/4/888/6115756)).
## Citing this work
