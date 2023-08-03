# FROM opencpu/ubuntu-18.04:v2.1.2
FROM opencpu/base:v2.2.8-5

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 4242737E1C226BE8
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y libmariadbclient-dev
RUN R -e "install.packages('cli', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('devtools', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('plotly', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages(c('fastcluster', 'reshape', 'stringr'), repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('RMySQL', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('gplots', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install()"
RUN R -e "BiocManager::install(c('Biobase', 'KEGGREST', 'graphite', 'marray', 'limma'))"
RUN R -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db', 'org.Rn.eg.db'))"
RUN R -e "BiocManager::install(c('gtools', 'gdata', 'qvalue', 'fgsea', 'preprocessCore'))"
RUN R -e "BiocManager::install(c('SPIA', 'IRanges', 'genefilter'))"
RUN R -e "BiocManager::install(c('ctc'))"
RUN R -e "devtools::install_github('uc-bd2k/CLEAN.Hs')"
RUN R -e "devtools::install_github('uc-bd2k/CLEAN.Mm')"
RUN R -e "devtools::install_github('uc-bd2k/CLEAN.Rn')"
RUN R -e "devtools::install_github('uc-bd2k/CLEAN')"
RUN apt-get install -y opencpu-cache
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
RUN dpkg-reconfigure --frontend noninteractive tzdata

EXPOSE 80
EXPOSE 443
EXPOSE 8006
EXPOSE 8007
CMD service cron start && service nginx start && apachectl -DFOREGROUND
