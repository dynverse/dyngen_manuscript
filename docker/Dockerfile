FROM dynverse/dyngen:devel

ARG GITHUB_PAT

# install r deps
RUN apt-get update && \
   apt-get install -y libglpk-dev

# install miniconda for python dependencies
RUN Rscript -e 'reticulate::install_miniconda()'

# install more dependenies related to running the different usecases
RUN Rscript -e 'remotes::install_cran(c("ggpubr", "pheatmap", "RANN", "dyndimred", "ggbeeswarm", "grImport"))' && \
  Rscript -e 'remotes::install_github("dynverse/dynfeature@devel")' && \
  Rscript -e 'remotes::install_github("dynverse/dynwrap@devel")' && \
  Rscript -e 'remotes::install_github("dynverse/dynplot@devel")'


# install dtw
RUN Rscript -e 'remotes::install_cran("dtw")'

# install anndata
RUN Rscript -e 'remotes::install_github("rcannood/anndata")' && \
  Rscript -e 'anndata:::install_anndata()'
  
# install velocyto & scvelo
RUN Rscript -e 'reticulate::py_install("Cython", pip = TRUE)' && \
  Rscript -e 'reticulate::py_install("velocyto", pip = TRUE)' && \
  Rscript -e 'remotes::install_github("dynverse/scvelo@devel")' && \
  Rscript -e 'scvelo:::install_scvelo()'
  
# install genie3bis (contains metric), pyscenic & arboreto
RUN Rscript -e 'remotes::install_github("rcannood/GENIE3bis")' && \
  Rscript -e 'reticulate::py_install("git+https://github.com/aertslab/pySCENIC.git", pip = TRUE)' && \
  Rscript -e 'reticulate::py_install("git+https://github.com/tmoerman/arboreto.git", pip = TRUE)'

# clone dyngen_manuscript and install dependencies
ADD . /home/rstudio/dyngen_manuscript
RUN cd /home/rstudio && \
  rm getting_started.R && \
  Rscript -e 'devtools::install("dyngen_manuscript/package", dependencies = TRUE)' && \
  rm -r dyngen_manuscript
