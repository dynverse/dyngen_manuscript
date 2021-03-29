library(tidyverse)
library(googledrive)
library(rmarkdown)
library(textreadr)

if (Sys.info()[["user"]] == "rcannood") {
  # httr::set_config(httr::config(http_version = 0)) # avoid http2 framing layer bug
  drive <- drive_download(as_id("1gU6--vq988dQS3qYQRlIqC1QO96A4XM81s54K6YJKIk"), type = "text", overwrite = TRUE, path = tempfile())

  # read docx
  c(readr::read_lines(drive$local_path), readr::read_lines("manuscript/vignettes.Rmd")) %>%
    str_replace_all("([^ ])(\\[@[^\\]]*\\])", "\\1 \\2") %>% # add spaces before citations
    str_replace_all("^#", "\n#") %>%
    write_lines("manuscript/render.Rmd")
}

# render the manuscript
render("manuscript/render.Rmd")

system("pdftk manuscript/render.pdf cat 1-16 output manuscript/manuscript.pdf")
system("pdftk manuscript/render.pdf cat 17-45 output manuscript/supplementary_files.pdf")

# system("latexdiff manuscript/v3/render.tex manuscript/render.tex > manuscript/track_changes.tex")
# file.remove("manuscript/render.pdf")
