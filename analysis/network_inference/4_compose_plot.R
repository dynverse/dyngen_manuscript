library(tidyverse)

explots <- read_rds("fig/network_inference/casewise_grn.rds")
summplots <- read_rds("fig/network_inference/score_summary.rds")

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

ggsave("fig/network_inference/cni.pdf", g, width = 8, height = grnh + 1 * bch)

