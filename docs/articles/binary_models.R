## ---- fig.align='center', fig.width=10, fig.height=8, out.height="600px", fig.dpi = 600----
library(phylopath)

models <- list(
  A = DAG(C ~ M + P + G + D),
  B = DAG(C ~ P + G + D,      M ~ M),
  C = DAG(C ~ P + G + D,                P ~ M),
  D = DAG(C ~ P + G + D,      M ~ P,               G ~ P),
  E = DAG(C ~ P + G + D,                P ~ M,     G ~ P),
  F = DAG(C ~ P + G + D,                P ~ M + G),
  G = DAG(C ~ P + G + D,      M ~ P,    P ~ G),
  H = DAG(C ~ P + G + D,      M ~ P),
  I = DAG(C ~ P + G + D,      M ~ M,               G ~ P),
  J = DAG(C ~ P + G,          M ~ P,               G ~ D),
  K = DAG(C ~ P + G,                    P ~ M,     G ~ D),
  L = DAG(C ~ M + P + G + D,            P ~ M + G)
)

plot_model_set(models)

## ---- eval=FALSE---------------------------------------------------------
#  cichlids_results <- phylo_path_binary(models, cichlids, cichlids_tree, parallel = "SOCK")

## ---- echo=FALSE---------------------------------------------------------
message("15 rows were dropped because they contained NA values.")
message("Pruned tree to drop species not included in dat.")
cichlids_results <- phylopath:::cichlids_results

## ------------------------------------------------------------------------
summary(cichlids_results)

## ---- eval=FALSE---------------------------------------------------------
#  best_cichlids <- best(cichlid_results)

## ---- echo=FALSE---------------------------------------------------------
best_cichlids <- phylopath:::best_cichlids

## ------------------------------------------------------------------------
best_cichlids

## ---- fig.align='center', fig.width=8, fig.height=4, out.width="600px", fig.dpi = 300----
plot(best_cichlids, algorithm = 'kk')

