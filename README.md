# BStools
README.txt for BStools

**Citation:**

Wang X, et al. A trimming-and-retrieving alignment scheme for bisulfite sequencing data. Bioinformatics (2015) 31(12):2040-2.



1. Installation: All needed software are compliled in resource folder. The toolkit is ready to use once unzipped.
                 If users prefer to use the software (FastQc, brat, fastx, cutadapt) by themself, the path for these software can be specified in main function: RRBS.pl

2. BStools.pl: the main function from data quality control to methylation level calling

Example command line:

perl BStools.pl  -i example.fastq -t F -d ~/BStools/ -p SAMPLE -R ~/hg19/  -r chr.list.name.txt -H 2 -Q ~/bin/FastQC/./fastqc -M ~/bin/brat-1.2.4/./trim -K ~/bin/brat-1.2.4/./brat-large -J ~/bin/brat-1.2.4/./acgt-count -T retrieve -N 5 -n 0 

Note:
1) please refer to methyQA for detailed description on input preparations and paramter settings
2) the retrieving step is only appled when fixed trimming (-T) is chosen by users
3) If users only want to trim off based at the 5' end (for NextSeq platform, the sequencing reads usually have low quality at the 5' end, but not 3' end), please specify -n as 0, such that no bases from 3' end will be trimmed off.

Main Input:
1) example.fastq: the fastq sequencing reads 
2) chr.list.name.txt: The chromosome list file shows the location of the reference genome for the chromosomes that users would like to analyze. See chr.list.name.txt for example. Please note that the location of chromosomes references should be corresponding to the reference directory which is specified in options "-d".

Main output:	
1) _fastqc.zip: FastQc results
2) .CG.base: methylation level for all CpG sites
3) BS.ps: distribution plot of the bisulfite conversion rate
4) CG.coverage.table.txt: coverage distribution of CpG sites
5) mC.distribution.ps:  distribution plot of methylation levels across genome


3. DMR identification: we provide a example script using Bioconductor package BiSeq to detect DMRs between tumor and normal samples. Please refer to 'biseq.example.script.R' for more detail

4. DMR visualization: we provide a R fucntions and script for DMR visualization. An example plot is shown in DMR.example.jpeg

5. We provide an example dataset (Examples/example.fastq)and unix shell scripts (example.code.txt) for applying BStools to the example dataset.




**Citation:**

Wang X, et al. A trimming-and-retrieving alignment scheme for bisulfite sequencing data. Bioinformatics (2015) 31(12):2040-2.
