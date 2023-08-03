
FROM shiny.ilincs.org:5000/ilincsr_base:latest
# FROM shiny.ilincs.org:5000/ilincsr-opencpu:base

RUN echo "CXX20=clang++-10" > ~/.R/Makevars
RUN echo "CXX20STD=-std=c++20" >> ~/.R/Makevars
RUN echo "CXX20FLAGS=-fPIC -O3" >> ~/.R/Makevars

RUN R -e "install.packages('sessioninfo', repos = 'http://cran.us.r-project.org', lib='/usr/lib/R/library')"

COPY . /tmp/ilincsR
RUN [ -f /tmp/ilincsR/iLincsCor* ] && R CMD INSTALL -l /usr/lib/R/library $(ls /tmp/ilincsR/iLincsCor*) && rm /tmp/ilincsR/iLincsCor*
RUN cd /tmp && R CMD build ilincsR && rm -r ilincsR && R CMD INSTALL -l /usr/lib/R/library $(ls ilincsR*)
#RUN R -e "install.packages('/tmp/ilincsR', repos=NULL, type='source')"

RUN sed -i 's/300,/3000,/g' /etc/opencpu/server.conf
RUN sed -i 's/4e9,/2e10,/g' /etc/opencpu/server.conf
RUN sed -i 's/1e9,/1e10,/g' /etc/opencpu/server.conf
RUN sed -i 's/100,/1000,/g' /etc/opencpu/server.conf
RUN sed -i 's/60,/1200,/g' /etc/opencpu/server.conf
RUN sed -i 's/90,/1200,/g' /etc/opencpu/server.conf
RUN sed -i 's/"]/","ggplot2","plotly","Biobase","ilincsR"]/g' /etc/opencpu/server.conf
RUN sed -i 's/\$scheme\:\/\/\$host\/\;/\$scheme\:\/\/\$http_host\/\;/g' /etc/nginx/opencpu.d/ocpu.conf

# Start non-daemonized webserver
CMD service cron start && service nginx start && apachectl -DFOREGROUND
