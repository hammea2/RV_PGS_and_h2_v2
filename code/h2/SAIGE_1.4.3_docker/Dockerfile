# Use the pre-built SAIGE image as the base
FROM wzhou88/saige:1.4.3

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install additional dependencies required for bcftools, vcftools, and PLINK
RUN apt-get update && \
    apt-get install -y \
    wget unzip gzip tabix build-essential \
    libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev zlib1g-dev \
    autoconf automake libtool pkg-config git cmake \
    libboost-all-dev libcurl4-openssl-dev liblzma-dev \
    libgsl-dev libgmp-dev libprotobuf-dev protobuf-compiler \
    r-base r-base-dev && \
    apt-get clean

# Install bcftools from Github (v1.21, released 2024-9-12)
RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
    tar -xjf bcftools-1.21.tar.bz2 && \
    cd bcftools-1.21 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools-1.21 bcftools-1.21.tar.bz2

# Install vcftools from GitHub (v0.1.16, released 2018-8)
RUN git clone https://github.com/vcftools/vcftools.git && \
    cd vcftools && \
    ./autogen.sh && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf vcftools

# Install PLINK 1.9 (beta 7.7, released on 2024-10-22)
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip && \
    unzip plink_linux_x86_64_20241022.zip -d /usr/local/bin/ && \
    rm plink_linux_x86_64_20241022.zip && \
    chmod +x /usr/local/bin/plink

# Install PLINK 2.0 (alpha 6.5, released on 2024-12-22)
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241222.zip && \
    unzip plink2_linux_x86_64_20241222.zip -d /usr/local/bin/ && \
    rm plink2_linux_x86_64_20241222.zip && \
    chmod +x /usr/local/bin/plink2


# Default command to keep the container running or provide an interactive shell
CMD ["bash"]












