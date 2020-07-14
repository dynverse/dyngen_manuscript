library(tidyverse)
library(googledrive)
library(rmarkdown)
library(textreadr)

if (Sys.info()[["user"]] == "rcannood") {
  httr::set_config(httr::config(http_version = 0)) # avoid http2 framing layer bug
  drive <- drive_download(as_id("15ZkzzB-XSYFfZutO4YPkUIw6iUvxJS92fSPfzqR1cLA"), type = "text", overwrite = TRUE, path = tempfile())

  # read docx
  readr::read_lines(drive$local_path) %>%
    str_replace_all("([^ ])(\\[@[^\\]]*\\])", "\\1 \\2") %>% # add spaces before citations
    str_replace_all("^#", "\n#") %>%
    write_lines("manuscript/render.Rmd")
}

# render the manuscript
render("manuscript/render.Rmd", output_dir = "manuscript/")

system("pdftk manuscript/render.pdf cat 1-17 output manuscript/manuscript.pdf")
system("pdftk manuscript/render.pdf cat 18-25 output manuscript/supplementary_files.pdf")

file.remove("manuscript/render.pdf")
