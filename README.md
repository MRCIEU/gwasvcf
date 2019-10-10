# Reading, querying and writing GWAS summary data in VCF format

Some utilities to create and query GWAS VCF files. These files are used to store GWAS summary data, and are not used to store genotype data. They are used for the utility of being able to quickly extract metadata from the fixed fields.


## Reference datasets

Example vcf (GIANT 2010 BMI):

```
wget -O IEU-a-2.vcf.gz https://www.dropbox.com/s/9bl20x537315f7u/IEU-a-2.vcf.gz?dl=0
wget -O IEU-a-2.vcf.gz.tbi https://www.dropbox.com/s/dkrqqbeg3pswl1b/IEU-a-2.vcf.gz.tbi?dl=0
```

1kg European reference panel for LD:

```
wget -O ld_files.tgz https://www.dropbox.com/s/yuo7htp80hizigy/ld_files.tgz?dl=0
tar xzvf ld_files.tgz -C app/
rm ld_files.tgz
```

1kg vcf harmonised against human genome reference: (PENDING)

## Getting started

GWAS summary data can be stored in VCF format files for consistency and querying speed. See the [GWAS VCF specification](https://github.com/MRCIEU/gwas_vcf_spec) for more information.

This package includes:

- a wrapper around the [bioconductor/VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html) package, providing functions tailored to GWAS VCF for reading, querying, creating and writing GWAS VCF format files
- some LD related functions such as using a reference panel to extract proxies, create LD matrices and perform LD clumping
- functions for harmonising a dataset against the reference genome and creating GWAS VCF files.

Examples of what can be done are as follows.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
devtools::install_github("mrcieu/gwasvcftools")
```


### Reading in everything

```r
vcf <- readVcf("IEU-a-2.vcf.gz")
```

### Reading in with filters

The `query_vcf` function takes either a filename to a vcf file, or vcf object as the main argument. You can then query on `rsid`, `pval` or `chrompos`.

```r
vcf <- query_vcf("IEU-a-2.vcf.gz", rsid=c("rs3128126", "rs3121561", "rs3813193"))
vcf <- query_vcf("IEU-a-2.vcf.gz", pval=5e-8)
vcf <- query_vcf("IEU-a-2.vcf.gz", chrompos=c("1:1097291-1099437"))
```

Chaining together e.g.

```r
vcf <- query_vcf("IEU-a-2.vcf.gz", rsid=c("rs3128126", "rs3121561", "rs3813193")) %>%
    query_vcf(pval=5e-8)
```

It's possible to have multiple GWAS studies per vcf. You can specify specific GWAS studies to read in using e.g.

```r
vcf <- query_vcf("IEU-a-2.vcf.gz", rsid=c("rs3128126", "rs3121561", "rs3813193"), id="IEU-a-2")
```

### LD proxies

If a set of rsids are requested from a vcf but some are absent, a reference panel can be used to search for LD proxies, extract them, and align the effects and alleles against the original variants that were requested.

```r
vcf <- proxy_match(
    "IEU-a-2.vcf.gz", 
    rsid, 
    bfile, 
    proxies="yes", 
    tag_kb=5000, tag_nsnp=5000, tag_r2=0.6, threads=1)
```

You may also extract only the best available proxies even if the requested rsids are present, by using `proxies="only"`.

### A note about chrompos

The fastest way to query VCFs is by specifying chromosome and position. Can specify specific positions, or ranges. e.g.

```r
cp <- c("1:10000", "2:10000-20000")
```

or as a data frame

```r
cp <- data_frame(chrom=c(1,2), start=c(10000,10000), end=c(10000, 20000))
```

You can check what will be parsed out with:

```r
parse_chrompos(cp)
```

Querying by p-value or rsid is also possible but is slower as only chrompos is indexed. On Mac and Linux, rsid and p-value queries are performed by calls to bcftools. On Windows it uses VariantAnnotation directly, because bcftools binaries are not available. This is unfortunately somewhat slower. If many operations are being performed it might be faster to read in the whole dataset and perform queries that way.

### Converting from the VCF object

The VCF object is somewhat complex and you can read more about it in the 'VariantAnnotation package documentation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html). You can create various other formats that might be easier to use from it. For example, create a `GRanges` object which is great for fast chromosome-position operations

```r
vcf_to_granges(vcf)
```

Create a data frame (tibble, note that it will drop rsids)

```
vcf_to_granges(vcf) %>% as_tibble
```

Create the `exposure_dat` or `outcome_dat` objects for the TwoSampleMR package

```
format_from_vcf(vcf, "exposure")
```

or 

```
format_from_vcf(vcf, "outcome")
```


### Creating the VCF object from a data frame

If you have GWAS summary data in a text file or data frame, this can be converted to a VCF object.

#### From a data frame

You may want to first harmonise the data so that all the non-effect alleles are aligned to the human genome reference. To do this, obtain the 1000 genomes vcf manifest (link) to use as an allele reference, and run

```
ref <- read_reference("filename.vcf.gz")
gwas        # This needs to be in exposure_dat format
dat <- harmonise_against_ref(gwas, ref)
```

Once done, you can create the vcf

```
dat %$% create_vcf(...)
```



