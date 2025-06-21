FROM ubuntu:20.04
LABEL maintainer="https://github.com/ming-en-ho/longphase-s"
LABEL version="1.0.0"

RUN apt-get update && \
    apt-get install -y git g++ gcc autoconf make zlib1g-dev libbz2-dev liblzma-dev && \
    rm -rf /var/lib/apt/lists/* 

WORKDIR /opt/longphase-s
RUN git clone https://github.com/ming-en-ho/longphase-s.git /opt/longphase-s && \
    autoreconf -i && \
    ./configure && \
    make -j 4 && \
    rm -rf /opt/longphase-s/.git

ENV PATH="${PATH}":${HOME}/bin:/opt/longphase-s

CMD ["longphase-s", "somatic_haplotag", "--help"]
# docker run longphase-s:latest longphase-s somatic_haplotag --help