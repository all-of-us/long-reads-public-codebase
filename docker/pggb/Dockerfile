FROM continuumio/miniconda3

MAINTAINER William T Harvey

COPY ./environment.yml /

RUN conda env create -f /environment.yml && conda clean -a

ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/conda/envs/pggb/lib/

ENV PATH=/opt/conda/envs/pggb/bin/:${PATH}

RUN echo "source activate pggb" > ~/.bashrc