library(tidyverse)
library(googledrive)
library(rmarkdown)
library(textreadr)

# download as docx instead of txt, because otherwise comments get pushed into the document
httr::set_config(httr::config(http_version = 0)) # avoid http2 framing layer bug
drive <- drive_download(as_id("15ZkzzB-XSYFfZutO4YPkUIw6iUvxJS92fSPfzqR1cLA"), type = "docx", overwrite = TRUE, path = tempfile())

# read docx
textreadr::read_docx(drive$local_path, remove.empty = FALSE, trim = FALSE) %>%
  str_replace_all("([^ ])(\\[@[^\\]]*\\])", "\\1 \\2") %>% # add spaces before citations
  str_replace_all("^#", "\n#") %>%
  write_lines("render.Rmd")

# render the manuscript
render("render.Rmd")

system("pdftk render.pdf cat 1-17 output manuscript.pdf")
system("pdftk render.pdf cat 18-25 output supplementary_files.pdf")

file.remove("render.pdf")
