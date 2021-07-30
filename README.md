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
- Motif for peaks(used HOMER)
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
ssDRIPSeqAnalysis.py requires a configuration file in json format as input.

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
        //Effective genome size, it used for peak calling, normalization and peaks proportion calculation.
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
        //Distance upstream of the start site of the regions and distance downstream of the end site of the regions. Used for skew metaplot and sense/antisense metaplot.
        "MetaplotExtend":"1000",
        //TSS/TTS left and right extension distance. TSS extends a certain distance from left and right as a promoter region. TTS extends a certain distance from left and right as a terminator region.
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
The target file format is as follows:

If you don’t use spike-in, the target file will contain four columns. Each column is separated by Tab. The first column is the group name. The second column is the sample name. Repeat samples have the same group name. The third and fourth columns are paired-end sequencing data in fastq format.
```
SampleA	SampleA_rep1	SampleA_rep1_R1.fastq.gz	SampleA_rep1_R2.fastq.gz
SampleA	SampleA_rep2	SampleA_rep2_R1.fastq.gz	SampleA_rep2_R2.fastq.gz
SampleB	SampleB_rep1	SampleB_rep1_R1.fastq.gz	SampleB_rep1_R2.fastq.gz
SampleB	SampleB_rep2	SampleB_rep2_R1.fastq.gz	SampleB_rep2_R2.fastq.gz
```
If you use spike-in, the targe file will contain seven columns. The first four columns are the same as above. The fifth column is the sample name of the input. The sixth and seventh columns are the paired-end sequencing data of input.
```
SampleA	SampleA_rep1	SampleA_rep1_R1.fastq.gz	SampleA_rep1_R2.fastq.gz	SampleA_rep1_input	SampleA_rep1_input_R1.fastq.gz	SampleA_rep1_input_R2.fastq.gz
SampleA	SampleA_rep2	SampleA_rep2_R1.fastq.gz	SampleA_rep2_R2.fastq.gz	SampleA_rep2_input	SampleA_rep2_input_R1.fastq.gz	SampleA_rep2_input_R2.fastq.gz
SampleB	SampleB_rep1	SampleB_rep1_R1.fastq.gz	SampleB_rep1_R2.fastq.gz	SampleB_rep1_input	SampleB_rep1_input_R1.fastq.gz	SampleB_rep1_input_R2.fastq.gz
SampleB	SampleB_rep2	SampleB_rep2_R1.fastq.gz	SampleB_rep2_R2.fastq.gz	SampleB_rep2_input	SampleB_rep2_input_R1.fastq.gz	SampleB_rep2_input_R2.fastq.gz
```
## Outputs
```
├── SampleA_rep1
│   ├── SampleA_rep1_fwd_nucleus_norm.bw
│   ├── SampleA_rep1_fwd_peaks.bed
│   ├── SampleA_rep1_nucleus_norm.bw
│   ├── SampleA_rep1_peaks.bed
│   ├── SampleA_rep1_pip.sh #Scripts used in BaseAnalysis
│   ├── SampleA_rep1_rev_nucleus_norm.bw
│   ├── SampleA_rep1_rev_peaks.bed
│   ├── SampleA_rep1_scale.xls
├── SampleA_rep2
├── SampleB_rep1
├── SampleB_rep2
├── test_deseq #Analysis results of DESeq2. *_anno.xls is the result of difference analysis. *_norm_anno.xls is the normalized R-loop abundance. *_counts_final_anno.xls is the original read counts.
│   ├── all #No split strand
│   │   ├── merge_test.bed
│   │   ├── test_Callus_Flagleaf_diffexpr_results_anno.xls
│   │   ├── test_Callus_Flagleaf_diffexpr_results.xls
│   │   ├── test_Callus_Spike_diffexpr_results_anno.xls
│   │   ├── test_Callus_Spike_diffexpr_results.xls
│   │   ├── test_counts_final_anno.xls
│   │   ├── test_counts_final.xls
│   │   ├── test_counts.npz
│   │   ├── test_counts.xls
│   │   ├── test_deseq.r #Scripts used in DESeq2
│   │   ├── test_deseq_scale.txt
│   │   ├── test_norm_anno.xls
│   │   ├── test_norm.xls
│   │   ├── test_pca_deseq2.pdf
│   │   ├── test_Spike_Flagleaf_diffexpr_results_anno.xls
│   │   └── test_Spike_Flagleaf_diffexpr_results.xls
│   ├── fwd #Fwd strand
│   └── rev #Rev strand
├── test_analysis # The results of DownstreamAnalysis
│   ├── cluster
│   ├── correlation
│   ├── motif
│   ├── peaks_content_distribution
│   │   ├── test_peaks_content_distribution.xls
│   │   ├── tss_150.bed
│   │   └── tts_150.bed
│   ├── peaks_length_distribution
│   │   └── test_peaks_length_distribution.npy
│   ├── sense_antisense
│   │   ├── SampleA_rep1_antisense.gz
│   │   ├── SampleA_rep1_sense.gz
│   │   ├── SampleA_rep2_antisense.gz
│   │   ├── SampleA_rep2_sense.gz
│   │   ├── SampleB_rep1_antisense.gz
│   │   ├── SampleB_rep1_sense.gz
│   │   ├── SampleB_rep2_antisense.gz
│   │   └── SampleB_rep2_sense.gz
│   └── skew
│       ├── SampleA_rep1_fwd_GCATSkew.gz
│       ├── SampleA_rep1_GCATSkew.gz
│       ├── SampleA_rep1_rev_GCATSkew.gz
│       ├── SampleA_rep2_fwd_GCATSkew.gz
│       ├── SampleA_rep2_GCATSkew.gz
│       ├── SampleA_rep2_rev_GCATSkew.gz
│       ├── SampleB_rep1_fwd_GCATSkew.gz
│       ├── SampleB_rep1_GCATSkew.gz
│       ├── SampleB_rep1_rev_GCATSkew.gz
│       ├── SampleB_rep2_fwd_GCATSkew.gz
│       ├── SampleB_rep2_GCATSkew.gz
│       ├── SampleB_rep2_rev_GCATSkew.gz
│       ├── gene_GCATSkew.gz
│       ├── test_ATSkew.bw
│       └── test_GCSkew.bw
├── test_stat.xls #Sample statistics
├── DripConfig.json
├── filter.txt #DripConfig[FilterChromFile]
├── test_ana.sh #Scripts used in DownstreamAnalysis
└── test_bwscale.xls #Scale factor for normalization
```
## Third-party software
The following software is used by this pipeline. When installing ssDripPipeline, these software will be installed automatically, and you do not need to install them yourself.
- [Bowtie2](https://github.com/BenLangmead/bowtie2)
- [Picard](https://broadinstitute.github.io/picard/)
- [Samtools](http://www.htslib.org/)
- [BEDTools](https://bedtools.readthedocs.io/en/latest/index.html)
- [deepTools](https://deeptools.readthedocs.io/en/develop/)
- [MACS2](http://github.com/taoliu/MACS/)
- [HOMER](http://homer.ucsd.edu/homer/)
- [DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [Mfuzz](https://www.bioconductor.org/packages/release/bioc/html/Mfuzz.html)
- [bedGraphToBigWig](http://rohsdb.cmb.usc.edu/goldenPath/help/bigWig.html)
## Others
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
