# dyngen manuscript

This repo contains all code necessary to reproduce the dyngen manuscript.

# Install dependencies
Clone the repository as follows.
```
# clone the repo
git clone https://github.com/dynverse/dyngen_manuscript
cd dyngen_manuscript
```

Open the R project file inside Rstudio and run the following code.
```r
# install the R package and its dependencies
devtools::install("package/", dep = TRUE)
```

# Reproducing an analysis
For each separate experiment, the code to reproduce the analysis can be found in the [analysis](analysis) folder. Any output generated by those scripts will be stored in the [results](results) folder. Note that each of the scripts assume that the working directory is the dyngen_manuscript folder.

# Build the manuscript
The complete manuscript (including supplementary files) can be built from the `render.Rmd` file, which can by build by running the following code inside R.

```r
library(rmarkdown)
render("render.Rmd")
```
