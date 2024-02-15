FROM continuumio/miniconda3

MAINTAINER Julie Wertz

COPY ./environment.yml /

RUN conda update -n base conda && \
        conda install -n base conda-libmamba-solver && \
        conda config --set solver libmamba && \
        conda env create -f /environment.yml && \
        conda clean -a

ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/conda/envs/augustus/lib/

ENV PATH=/opt/conda/envs/augustus/bin/:${PATH}

RUN echo "source activate augustus" > ~/.bashrc
