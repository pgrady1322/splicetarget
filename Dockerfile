FROM python:3.11-slim

LABEL maintainer="Patrick Grady"
LABEL description="splicetarget: Long-read splice isoform analysis and ASO target nomination"

# System deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    curl \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libncurses5-dev \
    && rm -rf /var/lib/apt/lists/*

# Install minimap2
RUN wget -q https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
    && tar -xjf minimap2-2.28_x64-linux.tar.bz2 \
    && cp minimap2-2.28_x64-linux/minimap2 /usr/local/bin/ \
    && rm -rf minimap2-2.28*

# Install samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2 \
    && tar -xjf samtools-1.19.tar.bz2 \
    && cd samtools-1.19 && ./configure --prefix=/usr/local && make -j$(nproc) && make install \
    && cd .. && rm -rf samtools-1.19*

WORKDIR /app
COPY . /app

RUN pip install --no-cache-dir -e ".[dev]"

ENTRYPOINT ["splicetarget"]
