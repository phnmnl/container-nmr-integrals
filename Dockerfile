FROM artemklevtsov/r-alpine:3.3.1

MAINTAINER PhenoMeNal-H2020 Project ( phenomenal-h2020-users@googlegroups.com )

LABEL Description="R-based script to compute NMR spectrum deviation for sample QC"
LABEL software.version="0.1.1"
LABEL version="0.1"
LABEL software="nmr-integrals"

LABEL base.image="artemklevtsov/r-alpine:3.3.1"
LABEL website="https://github.com/phnmnl/container-nmr-integrals"
LABEL documentation="https://github.com/phnmnl/container-nmr-integrals"
LABEL license="https://github.com/phnmnl/container-nmr-integrals"
LABEL tags="Metabolomics"

RUN \
     apk update \
  && apk add ca-certificates unzip wget \
  && rm -rf /var/cache/apk/*

RUN Rscript -e 'install.packages(c("caTools", "optparse"))'

WORKDIR /nmr-integrals

COPY \
  integrals.R \
  integrals_galaxy_wrapper.sh \
  runTest1.sh \
/usr/local/bin/

RUN \
  chmod a+rx \
     /usr/local/bin/integrals.R \
     /usr/local/bin/integrals_galaxy_wrapper.sh \
     /usr/local/bin/runTest1.sh

ENTRYPOINT ["integrals_galaxy_wrapper.sh"]
