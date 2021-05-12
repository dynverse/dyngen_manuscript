library(tidyverse)
library(googledrive)
library(rmarkdown)
library(textreadr)

if (Sys.info()[["user"]] == "rcannood") {
  # httr::set_config(httr::config(http_version = 0)) # avoid http2 framing layer bug
  drive <- drive_download(as_id("1aJTm23fUmUH4ad7Hl3eHoNFfK4Ltqj4LwLfP6pwmfi4"), type = "text", overwrite = TRUE, path = tempfile())

  # read docx
  readr::read_lines(drive$local_path) %>%
    str_replace_all("([^ ])(\\[@[^\\]]*\\])", "\\1 \\2") %>% # add spaces before citations
    str_replace_all("^#", "\n#") %>%
    write_lines("manuscript/manuscript.Rmd")
}

# take the relevant citatiosn and write it to a second file
bib_keys <-
  readr::read_lines("manuscript/manuscript.Rmd") %>%
  paste(collapse = "\n") %>%
  str_extract_all("\\[@[^\\]]*\\]") %>%
  first() %>%
  gsub("^\\[@", "", .) %>%
  gsub("\\]$", "", .) %>%
  str_split("; *") %>%
  unlist() %>%
  gsub("@", "", .) %>%
  unique

knitcitations::cite_options(check.entries = FALSE)
bib <- knitcitations::read.bibtex("manuscript/library.bib")
knitcitations::write.bibtex(bib[bib_keys], file = "manuscript/knitcitations.bib")

# render the manuscript
render("manuscript/manuscript.Rmd")



render("manuscript/supplementary_files.Rmd")

# system("latexdiff manuscript/v3/render.tex manuscript/render.tex > manuscript/track_changes.tex")
# file.remove("manuscript/render.pdf")
