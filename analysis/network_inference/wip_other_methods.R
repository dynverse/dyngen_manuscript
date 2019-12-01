
library(pandaR)
# exprm <- as.matrix(expr)
regulators <- unique(model$feature_network$from)
targets <- colnames(exprm)
samples <- rownames(exprm)

# motif <- crossing(tf = regulators, gene = targets) %>% mutate(score = 1) %>% as.data.frame()
#
# out <- lioness(motif, t(exprm), progress = TRUE, hamming = .1)
#
#
# data(pandaToyData)
# linonessRes <- lioness(pandaToyData$motif,
#                        pandaToyData$expression[,1:20],pandaToyData$ppi,hamming=.1,progress=FALSE)
#

base <- dynutils::calculate_similarity(
  expr[,regulators],
  expr[,targets],
  method = "spearman",
  margin = 2
)

sample_i <- 1
sample_sp <- dynutils::calculate_similarity(
  expr[-sample_i,regulators],
  expr[-sample_i,targets],
  method = "spearman",
  margin = 2
)

pred <-
  reshape2::melt(base - sample_sp, varnames = c("from", "to"), value.name = "importance") %>%
  arrange(desc(importance)) %>%
  head(10000)

ground_truth_sample <-
  with(
    dataset$feature_network_sc[sample_i,,drop=FALSE] %>% Matrix::summary(),
    tibble(
      from = factor(dataset$feature_network$from[j], levels = regulators),
      to = factor(dataset$feature_network$to[j], levels = targets),
      gold_effect = dataset$feature_network$effect[j],
      gold_value = x,
      gold = 1
    )
  )


pred2 <- pred %>%
  left_join(ground_truth_sample, by = c("from", "to")) %>%
  mutate(gold = gold %|% 0, gold_value = gold_value %|% 0, gold_effect = gold_effect %|% 0)


eval_pyscenic_sample <- GENIE4::evaluate_ranking_direct(
  values = pred2$importance,
  are_true = pred2$gold,
  num_positive_interactions = nrow(ground_truth_sample),
  num_possible_interactions = length(regulators) * length(targets)
)

