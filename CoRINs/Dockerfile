FROM python:3.7.7

RUN apt-get update
RUN apt-get install -y r-base
RUN Rscript -e 'install.packages("packrat")'

ENV PYTHONUNBUFFERED 1
ENV PIP_NO_CACHE_DIR 1
WORKDIR /var/www/html/server_corins/
COPY . /var/www/html/
RUN pip install -r requirements.txt
