# R with the GSbench package preinstalled.
#
# Run:  docker run --rm -it ghcr.io/mqfarooqi1/GSbench
FROM rocker/r-ver:4.5.2

LABEL org.opencontainers.image.title="GSbench"
LABEL org.opencontainers.image.description="R with GSbench preinstalled: genomic selection benchmarking (GBLUP, machine learning and a stacked ensemble), including the optional ML back-ends."
LABEL org.opencontainers.image.source="https://github.com/mqfarooqi1/GSbench"
LABEL org.opencontainers.image.licenses="MIT"

RUN Rscript -e 'install.packages(c("withr","glmnet","ranger","xgboost","rrBLUP"), repos = "https://cloud.r-project.org")'

COPY . /build/GSbench
RUN R CMD INSTALL --clean /build/GSbench && rm -rf /build
RUN Rscript -e 'library(GSbench); cat("GSbench installed\n")'

WORKDIR /work
CMD ["R"]
