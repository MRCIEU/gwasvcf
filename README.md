# Reading, querying and writing GWAS summary data in VCF format

<!-- badges: start -->
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![codecov](https://codecov.io/github/mrcieu/gwasvcf/branch/master/graphs/badge.svg)](https://codecov.io/github/mrcieu/gwasvcf)
[![R build status](https://github.com/MRCIEU/gwasvcf/workflows/R-CMD-check/badge.svg)](https://github.com/MRCIEU/gwasvcf/actions)
<!-- badges: end -->

Complete GWAS summary datasets are now abundant. A large repository of curated, harmonised and QC'd datasets is available in the [IEU GWAS database](https://gwas.mrcieu.ac.uk/). They can be queried via the [API](http://gwasapi.mrcieu.ac.uk/docs/) directly, or through the [ieugwasr](https://github.com/mrcieu/ieugwasr) R package, or the [ieugwaspy](https://github.com/mrcieu/ieugwaspy) python package. However, for faster querying that can be used in a HPC environment, accessing the data directly and not through cloud systems is advantageous. 

We developed a format for storing and harmonising GWAS summary data known as [GWAS VCF format](https://github.com/MRCIEU/gwas-vcf-specification/releases/tag/1.0.0) which can be created using [gwas2vcf](https://github.com/mrcieu/gwas2vcf). All the data in the [IEU GWAS database](https://gwas.mrcieu.ac.uk/) is available for download in this format. This R package provides fast and convenient functions for querying and creating GWAS summary data in GWAS VCF format (v1.0). See also [pygwasvcf](https://github.com/mrcieu/pygwasvcf) a Python3 parser for querying GWAS VCF files.

This package includes:

- a wrapper around the [bioconductor/VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html) package, providing functions tailored to GWAS VCF for reading, querying, creating and writing GWAS VCF format files
- some LD related functions such as using a reference panel to extract proxies, create LD matrices and perform LD clumping
- functions for harmonising a dataset against the reference genome and creating GWAS VCF files.

See also the [gwasglue](https://github.com/MRCIEU/gwasglue) R package for methods to connect the VCF data to Mendelian randomization, colocalisation, fine mapping etc.

## Installation

```r
remotes::install_github("mrcieu/gwasvcf")
```

## Usage

See vignettes here: [https://mrcieu.github.io/gwasvcf](https://mrcieu.github.io/gwasvcf).

## Citation

If using GWAS-VCF files please reference the studies that you use and the following paper:

**The variant call format provides efficient and robust storage of GWAS summary statistics.** Matthew Lyon, Shea J Andrews, Ben Elsworth, Tom R Gaunt, Gibran Hemani, Edoardo Marcora. bioRxiv 2020.05.29.115824; doi: https://doi.org/10.1101/2020.05.29.115824 


## Reference datasets

Example GWAS VCF (GIANT 2010 BMI):

- [http://fileserve.mrcieu.ac.uk/vcf/IEU-a-2.vcf.gz](http://fileserve.mrcieu.ac.uk/vcf/IEU-a-2.vcf.gz)
- [http://fileserve.mrcieu.ac.uk/vcf/IEU-a-2.vcf.gz.tbi](http://fileserve.mrcieu.ac.uk/vcf/IEU-a-2.vcf.gz.tbi)

1000 genomes reference panels for LD for each superpopulation - used by default in OpenGWAS:

- [http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz](http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz)

1000 genomes European reference panel for LD (legacy):

- [http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz](http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz)


1000 genomes vcf harmonised against human genome reference:

- [http://fileserve.mrcieu.ac.uk/vcf/1kg_v3_nomult.vcf.gz](http://fileserve.mrcieu.ac.uk/vcf/1kg_v3_nomult.vcf.gz)
- [http://fileserve.mrcieu.ac.uk/vcf/1kg_v3_nomult.vcf.gz.tbi](http://fileserve.mrcieu.ac.uk/vcf/1kg_v3_nomult.vcf.gz.tbi)

---

### Notes

#### Example data

data.vcf.gz and data.vcf.gz.tbi are the first few rows of the Speliotes 2010 BMI GWAS

The eur.bed/bim/fam files are the same range as data.vcf.gz, from here http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz
