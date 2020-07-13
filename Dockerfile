FROM ubuntu:20.04

MAINTAINER Kristian Peters (kpeters@ipb-halle.de)

# Add cran R backport
ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get -y update && apt-get -y dist-upgrade && apt-get -y install ca-certificates apt-transport-https perl locales gnupg2
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

# Generate locales
ENV LC_ALL="en_US.UTF-8"
ENV LC_CTYPE="en_US.UTF-8"
RUN locale-gen $LC_ALL
RUN dpkg-reconfigure locales

# Install packages
RUN apt-get -y update && apt-get -y dist-upgrade && apt-get -y --allow-unauthenticated install apt-transport-https make gcc gfortran g++ libblas-dev liblapack-dev libxml++2.6-dev libexpat1-dev libxml2-dev libnetcdf-dev libssl-dev pkg-config wget curl git unzip zip python3 python3-pip r-base r-base-dev pandoc pandoc-data

# Install R packages
RUN R -e 'install.packages(c("irlba","igraph","ggplot2","digest","lattice","XML","Rcpp","reshape2","plyr","stringi","stringr","intervals","devtools","RColorBrewer","plyr","RANN","knitr","ncdf4","microbenchmark","RUnit","foreach","doMC","curl","jsonlite","treemap","colourpicker","htmltools","shiny","shinyBS","shinyjs","DT","FactoMineR","slam","cba","squash","plotrix","plotly","circlize","matrixStats","Matrix","tools","ape","data.tree","plyr","data.table","kableExtra","htmlwidgets","sunburstR","SparseM","ROCR","PRROC","klaR","e1071","kohonen","nnet","rda","caret","caretEnsemble","RCurl","ontologyIndex"), repos="https://cloud.r-project.org/")'

# Install Bioconductor
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("multtest","MSnbase","mzR","MassSpecWavelet","S4Vectors","BiocStyle","faahKO","msdata","xcms","CAMERA","mixOmics","pcaMethods"), ask=FALSE)'

# Install development versions of R packages
RUN R -e 'devtools::install_github("https://github.com/MassBank/RMassBank")' && \
    R -e 'devtools::install_github(repo="CDK-R/rinchi@master")'

# Install outdated packages
WORKDIR /tmp
RUN wget "https://cran.r-project.org/src/contrib/Archive/rda/rda_1.0.2-2.1.tar.gz" && \
    R CMD INSTALL rda_1.0.2-2.1.tar.gz && \
    rm -f rda_1.0.2-2.1.tar.gz

# Cleanup
RUN apt-get -y --purge --auto-remove remove make gcc gfortran g++
RUN apt-get -y clean && apt-get -y autoremove && rm -rf /var/lib/{cache,log}/ /tmp/* /var/tmp/*

# Add data
ADD . /usr/local/share/mFam
WORKDIR /usr/local/share/mFam/data

