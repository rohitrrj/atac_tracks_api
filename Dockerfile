# Use the official R image
# https://hub.docker.com/_/r-base
FROM rstudio/plumber

# install dependency required by samtools
RUN apt-get install --yes wget libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libxml2-dev libcurl4-gnutls-dev libssl-dev libjpeg-dev libxt6

# Install R packages that are required
# TODO: add further package if you need!
RUN R -e "install.packages(c('jpeg','Hmisc','BiocManager','githubinstall','httr','jsonlite'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('biomaRt','GenomicFeatures','Gviz','rtracklayer','trackViewer','org.Hs.eg.db','TxDb.Hsapiens.UCSC.hg19.knownGene'))"
RUN R -e "library(githubinstall); githubinstall(c('plumber'),ask=F)"

# Create and change to the app directory.
WORKDIR /usr/src/app

# Copy local code to the container image.
COPY . .

# Set port on build time
ARG PORT=8080
ENV PORT=${PORT}
EXPOSE ${PORT}

# Run the web service on container startup.
CMD [ "Rscript", "server.R"]
