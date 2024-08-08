FROM ubuntu:22.04 AS base
RUN apt update && apt install -y wget build-essential
RUN wget https://sourceforge.net/projects/weinberg-r2r/files/R2R-1.0.7.tgz
RUN tar -xvzf R2R-1.0.7.tgz
RUN cd R2R-1.0.7 && ./configure && make && make install


FROM r-base
# RUN R -e "install.packages(c('ggplot2', 'shiny', 'shinyBS'))"
RUN R -e "install.packages(c('remotes', 'ggplot2', 'shiny', 'shinyBS', 'stringi'))"
RUN R -e "remotes::install_github('JPSieg/R2easyR')"

COPY --from=base /usr/local/bin/r2r /usr/local/bin/r2r
COPY --from=base /usr/local/bin/R2R_Stockholm.pm /usr/local/bin/R2R_Stockholm.pm
COPY --from=base /usr/local/bin/SelectSubFamilyFromStockholm.pl /usr/local/bin/SelectSubFamilyFromStockholm.pl

RUN mkdir -p /srv/shiny-server/www
COPY shinyapp/app.R /srv/shiny-server/app.R
# RUN chown -R shiny:shiny /srv/shiny-server/www
EXPOSE 8080
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app.R')"]
