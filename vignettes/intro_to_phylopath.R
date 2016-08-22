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

## ---- echo=FALSE---------------------------------------------------------
plot(result$average_model)

## ------------------------------------------------------------------------
plot(result$best_model)

## ------------------------------------------------------------------------
result$d_sep$one

