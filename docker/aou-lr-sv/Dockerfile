FROM continuumio/miniconda3:23.5.2-0

MAINTAINER Kiran V Garimella

# copy other resources
COPY ./environment.yml /

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/aou-lr-sv/bin/:${PATH}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/envs/aou-lr-sv/lib/

RUN apt-get -y update \
	&& apt-get -y install git make cmake protobuf-compiler gcc g++ zlib1g-dev libcurl4-openssl-dev libbz2-dev tree python3-pip liblzma-dev wget curl \
	&& apt-get clean

RUN echo "source activate aou-lr-sv" > ~/.bashrc
