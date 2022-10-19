FROM bioconductor/bioconductor_docker:devel

RUN R -e 'remotes::install_github("mrcieu/gwasvcf")'

# Get bcftools
ENV BCFTOOLS_VERSION 1.16
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -xf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    make && \
    mv bcftools /bin/ && \
    cd ../ && \
    rm -r bcftools-${BCFTOOLS_VERSION} bcftools-${BCFTOOLS_VERSION}.tar.bz2

# Get plink2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20220814.zip && \
    unzip plink2_linux_avx2_20220814.zip && \
    mv plink2 /bin/ && \
    rm plink2_linux_avx2_20220814.zip
