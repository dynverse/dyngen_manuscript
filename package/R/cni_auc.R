#' @importFrom rlang %|%
#'
#' @export
cni_auc <- function(dataset, model) {
  regulators <- dataset$regulators
  targets <- dataset$targets
  cell_ids <- dataset$cell_ids

  # evaluate static NI
  regulatory_network <-
    dataset$regulatory_network %>%
    mutate(gold = 1L) %>%
    rename(gold_strength = strength, gold_effect = effect)

  eval_static_static <- with(
    model$regulatory_network %>%
      left_join(regulatory_network, by = c("regulator", "target")) %>%
      mutate(gold = gold %|% 0L, gold_strength = gold_strength %|% 0, gold_effect = gold_effect %|% 0L) %>%
      as_tibble(),
    {
      evaluate_ranking_direct(
        values = strength,
        are_true = gold,
        num_positive_interactions = nrow(regulatory_network),
        num_possible_interactions = length(regulators) * length(targets)
      )
    }
  )
  static_static_auc <- eval_static_static$area_under[c("auroc", "aupr")]

  # evaluate casewise NI
  regulatory_network_sc <-
    dataset$regulatory_network_sc %>%
    mutate(gold = 1L) %>%
    rename(gold_strength = strength)

  casewise_casewise_auc <- map_df(
    cell_ids,
    function(cell_id) {
      gold_sc <- regulatory_network_sc %>%
        filter(cell_id == !!cell_id)
      reg_sc <-
        model$regulatory_network_sc %>%
        filter(cell_id == !!cell_id) %>%
        mutate(strength = strength + runif(n(), 0, 1e-8))

      eval_sc <- with(
        reg_sc %>%
          left_join(gold_sc, by = c("cell_id", "regulator", "target")) %>%
          mutate(gold = gold %|% 0L, gold_strength = gold_strength %|% 0) %>%
          as_tibble(),
        {
          evaluate_ranking_direct(
            values = strength,
            are_true = gold,
            num_positive_interactions = nrow(gold_sc),
            num_possible_interactions = length(regulators) * length(targets)
          )
        }
      )
      eval_sc$area_under %>% mutate(cell_id) %>% select(cell_id, everything())
    }
  )
  static_casewise_auc <- map_df(
    cell_ids,
    function(cell_id) {
      gold_sc <- regulatory_network_sc %>%
        filter(cell_id == !!cell_id)
      reg_sc <-
        model$regulatory_network

      eval_sc <- with(
        reg_sc %>%
          left_join(gold_sc, by = c("regulator", "target")) %>%
          mutate(gold = gold %|% 0L, gold_strength = gold_strength %|% 0) %>%
          as_tibble(),
        {
          evaluate_ranking_direct(
            values = strength,
            are_true = gold,
            num_positive_interactions = nrow(gold_sc),
            num_possible_interactions = length(regulators) * length(targets)
          )
        }
      )

      eval_sc$area_under %>% mutate(cell_id) %>% select(cell_id, everything())
    }
  )

  eval <- bind_rows(
    static_casewise_auc %>% mutate(method = "static_casewise"),
    static_static_auc %>% mutate(method = "static_static"),
    casewise_casewise_auc %>% mutate(method = "casewise_casewise")
  )

  # ggplot(eval) + geom_point(aes(auroc, aupr, colour = method))

  list(
    cc_auroc = mean(casewise_casewise_auc$auroc),
    cc_aupr = mean(casewise_casewise_auc$aupr),
    sc_auroc = mean(static_casewise_auc$auroc),
    sc_aupr = mean(static_casewise_auc$aupr),
    evals = list(eval)
  )
}


#' Evaluate a ranking
#'
#' @param values A vector of importance values of predicted interactions.
#' @param are_true A vector denoting whether the corresponding predicted interactions are true.
#' @param num_positive_interactions The total number of positives.
#' @param num_possible_interactions The total number ranked values.
#' @param extend_by The number of steps with which to fill the ranking as if random, if only a part of the ranking is given.
#'
#' @return A list containing two items, the ranked evaluation and the area under the curve scores
#'
#' @import dplyr
#' @importFrom pracma trapz
#'
#' @export
evaluate_ranking_direct <- function(
  values,
  are_true,
  num_positive_interactions,
  num_possible_interactions,
  extend_by = 10000
) {

  ord <- order(rank(values, ties.method = "random"), decreasing = TRUE)
  values <- values[ord]
  are_true <- are_true[ord]

  # calculate base statistics
  num_selected <- seq_along(are_true)
  tp <- cumsum(are_true)
  fp <- num_selected - tp
  length_ranking <- length(tp)
  num_negative_interactions <- num_possible_interactions - num_positive_interactions

  # extend base statistics, if necessary
  if (extend_by > 0 && length_ranking != num_possible_interactions) {
    diff.predictions <- num_possible_interactions - length_ranking
    diff.trues <- num_positive_interactions - tail(tp, 1)
    diff.negs <- num_negative_interactions - tail(fp, 1)

    multiplier <- seq_len(extend_by) / extend_by

    extra_num_selected <- multiplier * diff.predictions + tail(num_selected, 1)
    extra_tp <- multiplier * diff.trues + tail(tp, 1)
    extra_fp <- multiplier * diff.negs + tail(fp, 1)

    num_selected <- c(num_selected, extra_num_selected)
    are_true <- c(are_true, rep(NA, extend_by))
    tp <- c(tp, extra_tp)
    fp <- c(fp, extra_fp)
  }

  # calculate extended statistics
  metrics <- tibble(
    num_selected = c(0, num_selected),
    are_true = c(NA, are_true),
    tp = c(0, tp),
    fp = c(0, fp),
    fn = num_positive_interactions - tp,
    tn = num_negative_interactions - fp,
    acc = (tp + tn) / (num_positive_interactions + num_negative_interactions),
    tpr = tp / num_positive_interactions,
    spec = tn / num_negative_interactions,
    prec = ifelse(num_selected == 0, 1, tp / (tp + fp)),
    npv = tn / (tn + fn),
    f1 = 2 * tp / (2 * tp + fp + fn),
    mcc = ifelse(num_selected == 0, 0, (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))),
    informedness = tpr + spec - 1,
    markedness = prec + npv - 1
  )

  # calculate area under the curves
  area_under <- tibble(
    auroc = pracma::trapz(1 - metrics$spec, metrics$tpr),
    aupr = abs(pracma::trapz(metrics$tpr, metrics$prec)),
    F1 = ifelse(auroc + aupr != 0, 2 * auroc * aupr / (auroc + aupr), 0)
  )

  # generate output
  list(
    metrics = metrics,
    area_under = area_under
  )
}
