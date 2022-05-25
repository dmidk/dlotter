FROM python:3.9.7

ARG UID
ARG GID
ARG USER

ENV UID=$UID
ENV GID=$GID
ENV USER=$USER

RUN apt-get update \
    && apt-get install -yq libeccodes-tools libeccodes-dev

RUN apt-get install -y proj-bin libproj-dev libgeos++-dev libgeos-c1v5 libgeos-dev libgeos-doc
RUN apt-get -yq install unzip

RUN groupadd -g ${UID} ${USER}
RUN useradd -u ${UID} -g ${GID} -d /home/${USER}/ -m -s /bin/bash $USER && echo "${USER}:${USER}" | chpasswd

USER ${USER}

ENV PATH=$PATH:/home/${USER}/.local/bin/

#RUN pip install --upgrade pip
RUN mkdir -p ~/miniconda3 \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh \
    && bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3

ENV CONDA_DIR $HOME/miniconda3
ENV PATH=$CONDA_DIR/bin:$PATH
ENV PATH=/home/${USER}/miniconda3/bin:$PATH

WORKDIR /home/${USER}/
COPY dlotter.yml /home/${USER}/
COPY entrypoint.sh /home/${USER}/entrypoint.sh
COPY dlotter /home/${USER}/dlotter/

# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile

RUN conda env create -f dlotter.yml && conda clean -afy

ENV CARTOPY_DIR=/home/${USER}/.local/share/cartopy/

ENV NE_PHYSICAL=${CARTOPY_DIR}/shapefiles/natural_earth/physical
RUN mkdir -p ${NE_PHYSICAL} \ 
    && wget https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip -P ${CARTOPY_DIR} \
    && unzip ${CARTOPY_DIR}/ne_10m_coastline.zip -d  ${NE_PHYSICAL} \
    && rm ${CARTOPY_DIR}/*.zip

ENV NE_CULTURAL=${CARTOPY_DIR}/shapefiles/natural_earth/cultural
RUN mkdir -p ${NE_CULTURAL} \
    && wget https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_0_boundary_lines_land.zip -P ${CARTOPY_DIR} \
    && unzip ${CARTOPY_DIR}/ne_10m_admin_0_boundary_lines_land.zip -d  ${NE_CULTURAL} \
    && rm ${CARTOPY_DIR}/*.zip

ENV NE_CULTURAL_LOWRES=${CARTOPY_DIR}/shapefiles/natural_earth/cultural
RUN mkdir -p ${NE_CULTURAL_LOWRES} \
    && wget https://naturalearth.s3.amazonaws.com/50m_cultural/ne_50m_admin_0_boundary_lines_land.zip -P ${CARTOPY_DIR} \
    && unzip ${CARTOPY_DIR}/ne_50m_admin_0_boundary_lines_land.zip -d  ${NE_CULTURAL_LOWRES} \
    && rm ${CARTOPY_DIR}/*.zip

#DownloadWarning: Downloading: https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_0_boundary_lines_land.zip
# COPY entrypoint.sh /home/${USER}/entrypoint.sh

ENV ECCODES_DEFINITION_PATH /home/${USER}/dlotter/ec_definitions/:/home/${USER}/miniconda3/envs/dlotter/share/eccodes/definitions/

# ENTRYPOINT [ "/bin/bash", "entrypoint.sh" ]
ENV PATH /home/${USER}/miniconda3/envs/dlotter/bin:$PATH
CMD ["python", "-m", "dlotter"]
