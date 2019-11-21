library(tidyverse)
library(googledrive)
library(rmarkdown)

# copy zotero library if on rcannood's computer
if (Sys.info()[["user"]] == "rcannood") {
  file.copy("~/Workspace/library.bib", "library.bib", overwrite = TRUE)
}

# download Rmd from google drive
httr::set_config(httr::config(http_version = 0)) # avoid http2 framing layer bug
drive <- drive_download(as_id("1BOltsUJBQ8NhF5s9kU_Q470U0U00QblbG8HBgQwEzRE"), type="text/plain", overwrite=TRUE, path = tempfile())
read_lines(drive$local_path) %>%
  str_replace_all(" *(\\[@[^\\]]*\\])", "\\1") %>% # remove spaces before citations
  write_lines("manuscript.Rmd")

# render the manuscript
render("manuscript.Rmd")

# copy content to content.tex
lines <- read_lines("manuscript.tex")
content <- lines[seq(which(lines == "\\hypertarget{sec:dyngen-introduction}{%"), which(lines == "\\printbibliography")-1)]
write_lines(content, "content.tex")
