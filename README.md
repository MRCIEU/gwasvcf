# gwasvcftools

Some utilities to create and query GWAS BCF files. These files are used to store GWAS summary data, and are not used to store genotype data. They are used for the utility of being able to quickly extract metadata from the fixed fields.

## Example

In the file `/tests/harmonise_against_ref.r` shows an example script for creating a BCF file that is harmonised against a VCF reference file. The objective here is to use a common human genome build such that the reference allele is the same across all studies.

To run `/tests/harmonise_against_ref.r`:

```
Rscript harmonise_against_ref.r \
--ref_file /mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/1000g/1kg_v3_nomult.bcf \
--ref_build b37 \
--ref_info "1000 genomes phase 3 variants generated using https://github.com/MRCIEU/vcf-reference-datasets" \
--mrbase_id 2 \
--gwas_file /mnt/storage/private/mrcieu/research/mr-eve/gwas-files/2/elastic.gz \
--delimiter $'\t' \
--gzipped 1 \
--skip 0 \
--dbsnp_field 1 \
--ea_field 2 \
--nea_field 3 \
--ea_af_field 4 \
--effect_field 5 \
--se_field 6 \
--pval_field 7 \
--n_field 8 \
--out /mnt/storage/private/mrcieu/research/mr-eve/gwas-files/2/harmonised \
--out_type bcf
```


---

# GWAS BCF specification

Storing GWAS summary data in BCF format has the advantage that

- it uses a pre-existing, well known and well defined format
- many tools exist that can be used for manipulation
- binary format is fast and relatively small
- indexing makes looking up by chromosome and position extremely fast
- indexing time is very fast
- we can treat each GWAS as a distinct unit rather than storing everything in a database which is less nimble

## Metadata

- Contig information provides build and genome reference
- Provide MR-Base GWAS ID in the form of `##gwas.id=<mrbid>`
- Provide information about QC and harmonising against reference procedure in the form of `##counts.<metric>=<value>`


## Fields

To do

