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
system(str_glue("sed -i '1s/^.//' {drive$local_path}")) # remove first character, because this is some strange unicode character added by google
system(str_glue("sed -i 's/ *\\(\\[@[^\\]]*\\]\\)/\\1/' {drive$local_path}")) # remove spaces before citations
system(str_glue("sed -i 's/^\\(\\#\\#* \\)/\\n\\1/' {drive$local_path}")) # add enters before headings
system(str_glue("cat {drive$local_path} > manuscript.Rmd"))

# render the manuscript
render("manuscript.Rmd")

# copy content to content.tex
lines <- read_lines("manuscript.tex")
content <- lines[seq(which(lines == "\\maketitle")+1, which(lines == "\\printbibliography")-1)]
write_lines(content, "content.tex")
