# PrimerID
David Seifert (david.seifert@bsse.ethz.ch)

## Introduction
The programs and scripts provided here can be used to reproduce the results of the PrimerID paper **A Comprehensive Analysis of Primer IDs to Study Heterogeneous HIV-1 Populations** by Seifert et al. (2015).

## Prerequisites
If you wish to reproduce all results, you will need the following programs for compilation

1.  **Autoconf**; latest 2.69 release (http://www.gnu.org/software/autoconf/)

    GNU Autoconf produces the ./configure script from configure.ac.

2.  **Automake**; latest 1.15 release (http://www.gnu.org/software/automake/)

    GNU Automake produces the Makefile.in precursor, that is processed with ./configure to yield the final Makefile.

3.  **Autoconf Archive**; latest 2015.02.24 release (http://www.gnu.org/software/autoconf-archive/)

    Our configure.ac requires a number of m4 macros from the Autoconf archive.

4.  **GSL**; latest 1.16 release (http://www.gnu.org/software/gsl/)

    The GNU Scientific Library is required for random number generating and probability density functions.

5.  **Boost**; _at least 1.50_ (http://www.boost.org/)

    Boost provides the necessary abstractions for many different types.

6.  **SeqAn**; _at least 2.0_ (http://github.com/seqan/seqan/releases)

    SeqAn is the basis of the pIDalign semi-global aligner.

7.  **standard Unix utilities**; such as zip/unzip, wget, etc...

    If you cannot execute a command, chances are that you are missing one of the more common utilities we require in addition to the tools listed above.

A number of external programs are also required:

1.  **FastQC**; latest 0.11.3 release (http://www.bioinformatics.babraham.ac.uk/projects/download.html)

    FastQC produces the tile quality plots.

2.  **PRINSEQ**; latest 0.20.4 release (http://prinseq.sourceforge.net)

    PRINSEQ is a versatile perl script for trimming and quality filtering of NGS reads.

3.  **FASTX-Toolkit**; latest 0.0.14 release (http://hannonlab.cshl.edu/fastx_toolkit/download.html)

    FASTX-Toolkit is required for the _fastq_masker_ which masks internal bases of low quality.

4.  **AmpliconClipper**; latest release (http://github.com/SoapZA/AmpliconClipper)

    AmpliconClipper trims an alignment to a specific region, clipping reads in the process. It can remove the primer and ID regions from an alignment.

5.  **Picard tools**; 1.130 release (http://broadinstitute.github.io/picard/)

    The Picard toolset is used to transform the clipped SAM alignment back into raw fastq files.

6.  **bwa**; 0.7.12-r1044 release (http://bio-bwa.sourceforge.net)

    To emulate a standard viral haplotype reconstruction pipeline, we employ bwa as a standard and well-accepted aligner.

7.  **samtools**; 1.2 release (http://github.com/samtools/samtools/releases)

    Samtools is required for converting the SAM alignments produced by bwa to the BAM format.

8.  **HaploClique**; latest release (http://github.com/cbg-ethz/haploclique/releases)

    HaploClique is a cutting-edge haplotype reconstruction tool that we use to infer frequencies when not provided with tagged DNA fragments.

9.  **R**; latest release (http://www.r-project.org)

    R is used to produce most plots. In addition, you will require the following R packages:
    - _plotrix_
    - _RColorBrewer_
    - _sfsmisc_ (use the latest 1.0-28 version for proper LaTeX output)
    - _tikzDevice_ (optional, used instead of PDF output)
    - _VennDiagram_

10. **MATLAB**; R2015a release (http://www.mathworks.com/products/matlab)

    MATLAB is used to perform the likelihood ratio tests for comparing biases and variances between estimators.

Furthermore, you will require a compiler that can handle **C++0x** (which includes all C++11 compilers). Our pipeline has been developed on OS X 10.10 and has employed the LLVM/Clang C++ toolchain provided by XCode.

## Folder structure
In order to ease reproduction, you should retain the folder structure as given by the git repository
```
scripts/
├── Alignments
├── Analysis
├── Comparing estimators
├── Figure 1 - qPCR Efficiency
├── Figure 2,S10 - Histograms
├── Figure 3,S18 - Nucleotide distributions
├── Figure 4,5 - pID bias
├── Figure S17 - HMM LogLik
├── Figure S21,S22 - Bias and Variance
├── PreprocessedData
├── RawData
├── References
└── Table S7,S8 - Minimum heterozygous Hamming distances
```

## Preparing the programs
Compile the required programs by running
```
git clone https://github.com/cbg-ethz/PrimerID.git
cd PrimerID/
./autogen.sh
./configure SEQAN_INCLUDEDIR=<PATH TO SEQAN> CXX=clang++
make -j2
```
We have used Clang as CXX compiler here; GCC should work equally well.

## Retrieving the data
Download the data by performing
```
cd scripts/RawData/
wget https://n.ethz.ch/~dseifert/download/PrimerID.zip
unzip PrimerID.zip # you will be prompted for a password, which is given in the manuscript and is at this stage not public
cd ..
```

## Preprocessing
You will need to remove low-quality reads in order not to contaminate the analysis. Preprocess the data by running
```
cd PreprocessedData/
./preprocess.sh
cd ..
```
Take care to adjust the environmental flags _PRINSEQ_ and _FASTQ_MASKER_ in the bash script `preprocess.sh` to match the location of your executables.

## Alignment
Perform the alignment by doing
```
cd Alignments/
./align.sh
cd ..
```

## pIDalyse
The core of the analysis is based on pIDalyse. Run it by
```
mkdir Analysis/
cd Analysis/
../../pidalyse --r3223 ../References/5VM_3223_N_withPrimer.fasta --r3236 ../References/5VM_3236_N_withPrimer.fasta ../Alignments/32???/32???_nucMask_2.sam > OUTPUT.txt
cd ..
```
This will produce a plethora of information and statistics in **OUTPUT.txt**, like the clone frequencies and mutant frequency estimates.

## Reproducing statistics
### Figures
- Figure 1:
  ```
  cd Figure\ 1\ -\ qPCR\ Efficiency
  R CMD BATCH Figure1.R
  cd ..
  ```

- Figure 2 & S10 (**requires output of pIDalyse**):
  ```
  cd Figure\ 2,S10\ -\ Histograms
  R CMD BATCH Figure2_S10.R
  cd ..
  ```

- Figure 3 & S18 (**requires output of pIDalyse**, will require up to **1 h**):
  ```
  cd Figure\ 3,S18\ -\ Nucleotide\ distributions
  R CMD BATCH Figure3_S18.R
  cd ..
  ```

- Figure 4 & 5 (**requires output of pIDalyse**, will require up to **30 min**):
  Fitting the PCR bottleneck model requires a precomputed grid of parameter values in order to perform ordinary least squares fitting. While generating these precomputed values is possible, it requires enormous efforts that can only be replicated on a high performance computing cluster. We recommend using our precomputed values
  ```
  cd Figure\ 4,5\ -\ pID\ bias
  wget https://n.ethz.ch/~dseifert/download/RT-PCR_Simulation.zip
  unzip RT-PCR_Simulation.zip
  R CMD BATCH Figure4_5.R
  cd ..
  ```

- Figure S3-S8:
  ```
  cd RawData/
  ./analyse.sh
  cd ..
  ```
  Take care to adjust the environmental flag _FASTQC_ in the bash script `analyse.sh` to match the location of your executables.

- Figure S17 (**requires output of pIDalyse**):
  ```
  cd Figure\ S17\ -\ HMM\ LogLik
  R CMD BATCH FigureS17.R
  cd ..
  ```

- Figure S21 & S22 (**requires output of pIDalyse**):
  1. open MATLAB
  2. open the file `pIDalyse_DATA.m` (the frequencies should be identical to the frequencies from pIDalyse)
  3. evaluate all three matrices (i.e., load them into memory)
  4. change to the folder of the script by right-clicking on the script and selecting `Change Current Folder to [...]`
  5. open the file `PerformTests.m` and execute everything

### Tables
- Table 1:
  produced with **pIDalign** and **pIDalyse**

- Table 2:
  produced by **Figure4_5.R** and contained in `Figure4_5.Rout`

- Table 3:
  Raw and pID frequencies are produced by **pIDalyse**. HaploClique frequencies are produced by first running HaploClique over the trimmed alignments and then calculating frequencies:
  1. Run `prepare.sh` in **Comparing estimators** (Adjust the environmental variables _AMPLICONCLIPPER_, _PICARDTOOLS_ and _BWA_ to match the locations of your executables)
  2. Run `HC.sh` (this our local script, you will need to adjust this). HaploClique will likely require a larger server, as the algorithm is very computationally intensive.
  3. Compile the file `DetermineFreqs/determineHCfreq.cpp` using either GCC or Clang and give it the resulting `quasispecies.fasta` file to retrieve the frequencies of the 5 clones.