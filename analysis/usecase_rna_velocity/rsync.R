library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_rna_velocity")

local_dir <- exp$temporary()
remote_dir <- paste0("/scratch/irc/shared/dyngen_manuscript/", local_dir)

# UPLOAD
# qsub::rsync_remote(
#   remote_src = FALSE,
#   path_src = local_dir,
#   remote_dest = "prism",
#   path_dest = dirname(remote_dir)
# )

# DOWNLOAD
# qsub::rsync_remote(
#   remote_src = "prism",
#   path_src = remote_dir,
#   remote_dest = FALSE,
#   path_dest = dirname(local_dir)
# )
