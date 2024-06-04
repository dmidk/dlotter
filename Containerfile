FROM python:3.9.7

RUN apt-get update \
    && apt-get install -yq libeccodes-tools libeccodes-dev

RUN apt-get install -y proj-bin libproj-dev libgeos++-dev libgeos-c1v5 libgeos-dev libgeos-doc
RUN apt-get -yq install unzip

RUN curl -sSL https://pdm-project.org/install-pdm.py | python3 -

WORKDIR $HOME/dlotter

ADD dlotter $HOME/dlotter/dlotter/
COPY pyproject.toml $HOME/dlotter/pyproject.toml
COPY README.md $HOME/dlotter/README.md
COPY setup.py $HOME/dlotter/setup.py

ENV PATH="${PATH}:/root/.local/bin/"

RUN pdm install

ENV CARTOPY_DIR=${HOME}/.local/share/cartopy/

ENV NE_PHYSICAL=${CARTOPY_DIR}/shapefiles/natural_earth/physical
RUN mkdir -p ${NE_PHYSICAL} \
    && wget https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_coastline.zip -P ${CARTOPY_DIR} \
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

    https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_0_boundary_lines_land.zip
    https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_coastline.zip


    #DownloadWarning: Downloading: https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_0_boundary_lines_land.zip
# COPY entrypoint.sh /home/${USER}/entrypoint.sh

#ENV ECCODES_DEFINITION_PATH /home/${USER}/dlotter/ec_definitions/:/home/${USER}/miniconda3/envs/dlotter/share/eccodes/definitions/

# ENTRYPOINT [ "/bin/bash", "entrypoint.sh" ]
RUN ls -l
ENTRYPOINT ["pdm", "run", "python", "-m", "dlotter"]
CMD ["--help"]