## ----define models, fig.align='center', fig.width=10, fig.height=8, out.height="600px", fig.dpi = 600----
library(phylopath)

models <- define_model_set(
  A = c(C~M+D),
  B = c(C~D),
  C = c(C~D, P~M),
  D = c(C~D, M~P, G~P),
  E = c(C~D, P~M, G~P),
  F = c(C~D, P~M+G),
  G = c(C~D, M~P, P~G),
  H = c(C~D, M~P),
  I = c(C~D, M~M, G~P),
  J = c(M~P, G~D),
  K = c(P~M, G~D),
  L = c(C~M+D, P~M+G),
  .common = c(C~P+G)
)

plot_model_set(models, algorithm = 'kk')

## ----fit models----------------------------------------------------------
(cichlids_results <- phylo_path(models, cichlids, cichlids_tree))

