FROM dynverse/dyngen

ARG GITHUB_PAT

RUN cd /home/rstudio && \
  rm getting_started.R && \
  git clone --depth=1 --branch=master https://github.com/dynverse/dyngen_manuscript.git && \
  rm -rf dyngen_manuscript/.git && \
  Rscript -e 'devtools::install("dyngen_manuscript/package", dependencies = TRUE)'

RUN Rscript -e 'remotes::install_cran(c("dtw", "ggpubr", "pheatmap", "RANN"))' && \
  Rscript -e 'remotes::install_github("rcannood/GENIE3bis")' && \
  Rscript -e 'reticulate::install_miniconda()' && \
  Rscript -e 'reticulate::py_install("Cython", pip = TRUE)' && \
  Rscript -e 'reticulate::py_install("velocyto", pip = TRUE)' && \
  Rscript -e 'remotes::install_github("dynverse/scvelo@devel")' && \
  Rscript -e 'scvelo:::install_scvelo()' && \
  Rscript -e 'reticulate::py_install("git+https://github.com/aertslab/pySCENIC.git", pip = TRUE)' && \
  Rscript -e 'reticulate::py_install("git+https://github.com/tmoerman/arboreto.git", pip = TRUE)'
  
