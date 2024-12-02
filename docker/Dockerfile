FROM rocker/verse:4.3.3

RUN apt-get update -qq
RUN apt-get -y install cmake
RUN apt-get -y --no-install-recommends install libglu1-mesa-dev
RUN apt-get -y install libxml2-dev
RUN apt-get -y install ghostscript
RUN apt-get -y install libxt-dev
RUN apt-get -y install libfontconfig1-dev
RUN apt-get -y install libpng-dev
RUN apt-get -y install libtiff-dev
RUN apt-get -y install libharfbuzz-dev
RUN apt-get -y install libfribidi-dev
RUN apt-get -y install -y python3
RUN apt-get -y install python3-pip
RUN apt-get -y install python3-venv

RUN pip3 install numpy --target /usr/lib/python3.10
RUN pip3 install scipy --target /usr/lib/python3.10
RUN pip3 install matplotlib --target /usr/lib/python3.10
RUN pip3 install pandas --target /usr/lib/python3.10
RUN pip3 install mne --target /usr/lib/python3.10
RUN pip3 install psignifit --target /usr/lib/python3.10
RUN pip3 install pingouin --target /usr/lib/python3.10
RUN pip3 install statsmodels --target /usr/lib/python3.10

RUN install2.r –error \
--deps TRUE \
quarto \
reticulate \
osfr

RUN R 'reticulate::use_python("/usr/bin/python3.10")'
