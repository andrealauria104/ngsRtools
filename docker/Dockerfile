# Docker inheritance - bioconductor ----
FROM bioconductor/bioconductor_docker:RELEASE_3_10

# Labels ----
LABEL name="anlauria/ngsrtools" \
      version="0.0.1" \
      url="https://github.com/andrealauria104/ngsRtools" \
      mantainer="andrea.lauria104@gmail.com" \
      description="ngsRtools docker image"

# Install R dependencies ----
WORKDIR /opt

RUN git clone https://github.com/andrealauria104/ngsRtools.git \
	&& cd ngsRtools \
	&& ./scripts/install-dependencies.r

# Install ngsRtools packages/scripts ----
RUN cd ngsRtools \
	&& ./install.sh

WORKDIR /bin
RUN ln -s /opt/ngsRtools/scripts/compute-barcode-distance.r \
	&& ln -s /opt/ngsRtools/scripts/get-illumina-ssheet.r \
	&& ln -s /opt/ngsRtools/scripts/scrnaseq-qc.r \
	&& ln -s /opt/ngsRtools/scripts/xlsx-to-tab.r

# Install additional packages ----
WORKDIR /tmp

RUN Rscript -e 'install.packages("rmarkdown", repos = "http://cran.us.r-project.org")'

CMD ["/init"]

# To run interactive bash in current directory
# docker run --rm -it -u `id -u`:`id -g` -v $PWD:/tmp/ anlauria/ngsrtools:0.0.1 bash
# To run interactive R console in current directory
# docker run --rm -it -u `id -u`:`id -g` -v $PWD:/tmp/ anlauria/ngsrtools:0.0.1 R
# To run RStudio server in current directory - localhost:8787
# docker run --rm -v $PWD:/tmp/ -e PASSWORD="ngsrtools" -p 8787:8787 anlauria/ngsrtools:0.0.1
