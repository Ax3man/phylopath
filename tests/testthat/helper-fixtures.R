# Test fixtures.
#
# Everything here is fully deterministic: the tree is a fixed Newick string and
# the simulated traits are hard-coded. Nothing depends on the state of the RNG,
# so the numeric regression tests cannot drift if R's random number generation
# changes between versions.

# A 40 tip tree, used for the fast end-to-end tests. 40 species is enough for
# the method to reliably distinguish the causal models used in the tests (see
# test-phylo_path.R) while keeping every fit well under a second.
test_tree <- function() {
  ape::read.tree(text = paste0(
    "((((((t1:0.0231,t2:0.0231):0.0197,((t3:0.0305,(t4:0.01,t5:0.01):0.0205):0.0026,t6:0.033)",
    ":0.0098):0.1438,((t7:0.05,((t8:0.0044,t9:0.0044):0.0031,t10:0.0075):0.0424):0.1331,(((t1",
    "1:0.0053,t12:0.0053):0.02,(t13:0.0118,t14:0.0118):0.0135):0.0884,(t15:0.0348,t16:0.0348)",
    ":0.0789):0.0694):0.0036):0.2028,(((t17:0.0067,t18:0.0067):0.023,(t19:0.0102,t20:0.0102):",
    "0.0194):0.0388,t21:0.0684):0.321):0.4231,((((t22:0.0754,t23:0.0754):0.0215,(t24:0.0704,t",
    "25:0.0704):0.0266):0.0985,t26:0.1954):0.6058,((t27:0.0138,(t28:0.0004,t29:0.0004):0.0134",
    "):0.1488,((t30:0.0275,t31:0.0275):0.0844,(t32:0.1089,(t33:0.0377,(t34:0.0205,t35:0.0205)",
    ":0.0173):0.0711):0.0031):0.0506):0.6387):0.0113):0.6857,(((t36:0.0147,t37:0.0147):0.0126",
    ",t38:0.0274):0.0713,(t39:0.0258,t40:0.0258):0.0729):1.3995);"
  ))
}

# Traits simulated on test_tree() under a known causal model:
#   A -> B -> C   and   D as an isolate
# (Brownian motion, each child regressed on its parent with slope 0.9.)
# Any analysis of these data should prefer A -> B -> C over models that reroute
# the chain. Note that it cannot be expected to prefer it over C -> B -> A,
# which is Markov equivalent and so implies exactly the same independencies.
test_data <- function() {
  data.frame(
    A = c(-2.1632, -2.1944, -1.9804, -2.1172, -2.2491, -2.3869, -2.8867, -2.9909,
          -2.8937, -2.8402, -1.1486, -0.9810, -1.1674, -1.0710, -1.2592, -1.2429,
          -0.0547, -0.0522, -0.3782, -0.3166, -1.0537, -0.8547, -1.1727, -1.2058,
          -0.9253, -1.2423, -0.7987, -0.8433, -0.8264, -0.7792, -0.8515, -1.1272,
          -0.3744, -0.6173, -0.6120, 0.0470, -0.2118, -0.2143, -0.4367, -0.3689),
    B = c(-0.6159, -0.8486, -0.6949, -0.6256, -0.7745, -1.1004, -1.7052, -1.8911,
          -1.6873, -1.6934, -0.4794, -0.3508, -0.4983, -0.2842, -0.1015, 0.5763,
          1.0623, 0.9897, 0.7551, 0.8577, 0.6407, -0.4899, -0.2094, -0.8544,
          -0.1373, -0.3948, 0.5106, 0.3240, 0.3264, -0.0751, -0.3589, -0.5340,
          -0.9355, -0.7154, -0.7236, -0.1096, -0.4451, -0.9183, -0.8285, -0.8139),
    C = c(-1.7190, -1.5711, -2.0110, -1.5231, -1.9053, -2.2338, -3.4434, -2.9787,
          -2.8154, -3.0071, -1.4081, -1.2720, -1.3978, -1.2505, -1.0540, -0.8192,
          -0.9223, -0.8335, -1.6490, -1.1731, -1.5554, -1.7267, -1.1379, -1.7934,
          -1.2615, -1.6137, -0.5081, -0.8124, -0.8898, -0.8536, -0.9320, -1.3383,
          -1.8951, -1.5553, -1.7347, -1.9126, -2.3293, -2.4272, -2.0751, -1.9843),
    D = c(0.7025, 0.4487, 0.4846, 0.4542, 0.6304, 0.2842, 0.4295, 0.2718, 0.2803,
          0.3459, 0.3233, 0.0972, 0.1119, 0.1119, 0.8455, 0.5473, -0.1832, -0.1270,
          0.3126, 0.3482, 0.0431, -0.9322, -0.8491, -0.8712, -0.7164, -0.4835,
          -0.5353, -0.7493, -0.7838, 0.2687, 0.2895, 0.5555, 0.0161, -0.0004,
          0.0787, 2.0939, 1.9046, 1.9155, 1.4432, 1.6740),
    row.names = paste0("t", 1:40)
  )
}

# A trivially small tree/data pair, for input-validation tests where the actual
# numbers are irrelevant.
tiny_tree <- function() {
  ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
}

tiny_data <- function(...) {
  data.frame(..., row.names = c("a", "b", "c", "d"))
}

# Build a random DAG on `v` nodes, as an upper triangular adjacency matrix with
# the node order permuted, which guarantees acyclicity. Used by the property
# based tests. Caller is responsible for seeding.
random_DAG <- function(v, p = 0.4) {
  nm <- LETTERS[seq_len(v)]
  perm <- sample(nm)
  m <- matrix(0L, v, v, dimnames = list(perm, perm))
  for (i in seq_len(v - 1)) {
    for (j in seq(i + 1, v)) {
      if (stats::runif(1) < p) m[i, j] <- 1L
    }
  }
  m <- m[nm, nm]
  class(m) <- c("matrix", "array", "DAG")
  m
}

# Construct a fitted_DAG by hand, so that averaging can be checked against
# values computed on paper. `coefs`/`ses` are named like "AB" for the path A -> B.
make_fitted_DAG <- function(vars, coefs, ses) {
  m <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  cm <- sm <- m
  for (path in names(coefs)) {
    from <- substr(path, 1, 1)
    to   <- substr(path, 2, 2)
    cm[from, to] <- coefs[[path]]
    sm[from, to] <- ses[[path]]
  }
  # These are synthetic continuous coefficients, so record that: without `binary`
  # every function that reports them warns that their scale is unknown.
  structure(list(coef = cm, se = sm, binary = stats::setNames(rep(FALSE, length(vars)), vars)),
            class = "fitted_DAG")
}

# A deliberately under-powered subset: few enough species that CICc is undefined
# for the denser causal models (it needs n > q + 1).
small_tree <- function() ape::keep.tip(test_tree(), paste0("t", 1:8))

small_data <- function() test_data()[paste0("t", 1:8), ]
