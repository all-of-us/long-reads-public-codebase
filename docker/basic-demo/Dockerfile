FROM continuumio/miniconda3

MAINTAINER Kiran V Garimella

# install conda packages
RUN conda create --name my-env samtools minimap2 && conda clean -a

# install some basic packages
RUN apt-get -y update \
    && apt-get -y install make gcc g++ tree wget curl \
    && apt-get clean

# activate conda environment
RUN echo "source activate my-env" > ~/.bashrc
