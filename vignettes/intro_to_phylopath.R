## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(dev = "png")

## ------------------------------------------------------------------------
library(phylopath)

models <- list(
  one   = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ DD),
  two   = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ LS + DD),
  three = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ NL),
  four  = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ BM + NL),
  five  = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ BM + NL + DD),
  six   = DAG(LS ~ BM, NL ~ BM + RS, DD ~ NL, RS ~ BM),
  seven = DAG(LS ~ BM, NL ~ BM + RS, DD ~ NL, RS ~ LS + BM),
  eight = DAG(LS ~ BM, NL ~ BM + RS, DD ~ NL),
  nine  = DAG(LS ~ BM, NL ~ BM + RS, DD ~ NL, RS ~ LS)
)

## ------------------------------------------------------------------------
models$one

## ------------------------------------------------------------------------
plot(models$one)

## ------------------------------------------------------------------------
result <- phylo_path(models, data = rhino, tree = rhino_tree, 
                     order = c('BM', 'NL', 'DD', 'LS', 'RS'))

## ------------------------------------------------------------------------
result

## ------------------------------------------------------------------------
summary(result)

## ------------------------------------------------------------------------
(best_model <- best(result))

## ------------------------------------------------------------------------
plot(best_model)

## ------------------------------------------------------------------------
average_model <- average(result)
plot(average_model)

## ------------------------------------------------------------------------
average_model_full <- average(result, method = "full")
plot(average_model_full)

## ---- fig.height=6-------------------------------------------------------
coef_plot(best_model)

## ---- fig.width=6--------------------------------------------------------
coef_plot(average_model_full, reverse_order = TRUE) + 
  ggplot2::coord_flip() + 
  ggplot2::theme_bw()

## ------------------------------------------------------------------------
result$d_sep$one

