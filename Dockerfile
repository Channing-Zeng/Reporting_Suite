FROM stackbrew/ubuntu:14.04
MAINTAINER Justin Johnson "https://github.com/bioinfo

RUN apt-get update && apt-get install -y

RUN cd /usr/local/src && \
	git clone https://github.com/AstraZeneca-NGS/Reporting_Suite.git && \
	ln -s Reporting_Suite az.reporting
