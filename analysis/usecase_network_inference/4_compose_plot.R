library(tidyverse)
library(dyngen.manuscript)

exp <- start_analysis("usecase_network_inference")

explots <- read_rds(exp$result("casewise_grn.rds"))
summplots <- read_rds(exp$result("score_summary.rds"))

grnh <- 5
bch <- 2.5
g <- patchwork::wrap_plots(
  explots + labs(tag = "A"),
  summplots$casewise_casewise +
    scale_x_continuous(limits = c(.4, .7)) +
    scale_y_continuous(limits = c(0, .06)) +
    labs(tag = "B", title = NULL) +
    # theme_classic() +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold")),
  heights = c(grnh, bch),
  ncol = 1
) & theme(plot.tag.position = c(0,1))

ggsave(exp$result("cni.pdf"), g, width = 8, height = grnh + 1 * bch)

