FROM continuumio/miniconda3:23.5.2-0

MAINTAINER Kiran V Garimella

# install conda packages
ENV CONDA_ENV="my-env"
RUN conda create --name $CONDA_ENV -c bioconda samtools minimap2 && conda clean -a

# install some basic packages
RUN apt-get -y update \
    && apt-get -y install make gcc g++ tree wget curl \
    && apt-get clean

# activate conda environment
RUN echo "source activate $CONDA_ENV" > ~/.bashrc
