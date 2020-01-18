#### Install ####
requiredPackages <-
  c("igraph",
    "rstudioapi")

for (pac in requiredPackages) {
  if (!require(pac,  character.only = TRUE)) {
    install.packages(pac, repos = "http://cran.rstudio.com")
    library(pac,  character.only = TRUE)
  }
}

rm(pac)
rm(requiredPackages)

# set pwd to current directory, must load rstudioapi before.
if(rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

graph.fb <- read.graph("data/facebook_combined.txt", directed = FALSE)
graph.tw <- read.graph("data/twitter_combined.txt", format = "ncol", directed = TRUE)
graph.ep <- read.graph("data/soc-Epinions1.txt", directed = TRUE)

n <- 10000
ba <- barabasi.game(n, 1, directed=FALSE)
ws <- watts.strogatz.game(1, n, 3, 0.1)
er <- erdos.renyi.game(n, 0.05)

#degree.org <- degree(graph)
#betweeness.org <- betweenness(graph)
#pagerank.org <- page.rank(graph)
#eigen_centrality(graph)

plot_easy <- function(graph) {
  par(mar = c(0, 0, 0, 0))
  plot(graph, vertex.label=NA, vertex.size=3, edge.arrow.mode="-")
}

plot_with_color <- function(graph, infected, colors = c("white", "red")) {
  set.seed(123)
  par(mar = c(0, 0, 0, 0))
  colors <- colors[as.numeric(as_ids(V(er)) %in% infected)+1]
  plot(graph, vertex.label=NA, vertex.size=10, vertex.color=colors, edge.arrow.mode="-")
}

plot_history <- function(graph, history, steps = 10) {
  len <- length(history)
  indexes <- seq(1, len, by = floor((len - 1) / steps))
  for(infected in history) {
    plot_with_color(graph, infected)
  }
}

simulate <- function(graph, vac=integer(0), iter=20, inf=0.001, p_i=0.005, p_h=0.3) {
  n <- gorder(graph)
  infected <- sample(setdiff(seq(1, n), vac), round(n * inf), replace = FALSE)
  res <- c(length(infected))
  history <- list(infected)
  for (r in seq(iter)) {
    healed <- infected[sample(c(TRUE, FALSE), length(infected), prob = c(p_h, 1-p_h), replace = TRUE)]
    neighbours <- c()
    for (v in infected) {
      nei <- as.numeric(neighbors(graph, v, mode = "OUT"))
      nei <- setdiff(nei, append(infected, vac))
      neighbours <- append(neighbours, nei)
    }
    new_infected <- unique(as.numeric(neighbours[sample(c(TRUE, FALSE), length(neighbours), prob = c(p_i, 1-p_i), replace = TRUE)]))
    infected <- append(infected[-healed], new_infected)
    #history[[length(history)+1]] <- infected
    res <- append(res, c(length(infected)))
  }
  return(res)
}

sample_random_walk_edges <- function(graph, num, metropolis_hastings = FALSE, mode = "OUT") {
  sampled <- make_empty_graph(directed = FALSE)
  v <- sample(as.numeric(V(graph)), 1)
  v.id <- 1
  v_na <- NA
  sampled <- add_vertices(sampled, 1, name=v)
  for (iter in seq(num)){
    neighbours <- as.numeric(neighbors(graph, v, mode = mode))
    if (metropolis_hastings) {
      repeat {
        v_n <- sample(neighbours, 1)
        if (runif(1) < min(1, degree(graph, v)/degree(graph, v_n))) {
          break
        }
      }
    } else {
      v_n <- sample(neighbours, 1)
    }
    v_n.id <- match(v_n, as.numeric(V(sampled)))
    if (is.na(v_n.id)) {
      sampled <- add_vertices(sampled, 1, name=v_n)
      v_n.id <- gorder(sampled)
    }
    if(!are.connected(sampled, v.id, v_n.id)) {
      sampled <- add_edges(sampled, c(v.id, v_n.id))
    }
    v <- v_n
    v.id <- v_n.id
  }
  return(sampled)
}

sample_random_walk_nodes <- function(graph, num, metropolis_hastings = FALSE, mode = "OUT") {
  v <- sample(as.numeric(V(graph)), 1)
  nodes <- v
  for (iter in seq(num)){
    neighbours <- as.numeric(neighbors(graph, v, mode = mode))
    if (metropolis_hastings) {
      repeat {
        v_n <- sample(neighbours, 1)
        if (runif(1) < min(1, degree(graph, v)/degree(graph, v_n))) {
          v <- v_n
          break
        }
      }
    } else {
      v <- sample(neighbours, 1)
    }
    nodes <- append(nodes, v)
  }
  return(induced.subgraph(graph, nodes))
}

sample_expansion_snowball <- function(graph, num) {
  nodes <- sample(as.numeric(V(graph)), 1)
  neighbours <- as.numeric(neighbors(graph, nodes, mode = "OUT"))
  for (iter in seq(num)){
    max_n <- NA
    max_nei <- c()
    for (n in sample(neighbours)) {
      new_neighbourhood <- setdiff(append(as.numeric(neighbors(graph, n, mode = "OUT")), neighbours), append(nodes, n))
      if (length(new_neighbourhood) > length(max_nei)) {
        max_n <- n
        max_nei <- new_neighbourhood
      }
    }
    nodes <- append(nodes, max_n)
    neighbours <- max_nei
  }
  return(induced.subgraph(graph, nodes))
}

sample_expansion_mcmc <- function(graph, num) {
  n <- gorder(graph)
  p <- 10 * gsize(graph) / n * log10(n)
  nodes <- sample(as.numeric(V(graph)), 1)
  neighbours <- as.numeric(neighbors(graph, nodes, mode = "OUT"))
  max_nodes <- NA
  max_exp <- 0
  exp <- length(neighbours) / (n - length(nodes))
  for (iter in seq(num)){
    repeat {
      v <- sample(setdiff(1:n, nodes), 1)
      new_nodes <- append(nodes, v)
      neighbours <- unique(setdiff(append(neighbours, as.numeric(neighbors(graph, v))), nodes))
      exp_new <- length(neighbours) / (n - length(nodes))
      if (runif(1) < min(1, exp_new / exp) ** p) {
        nodes <- new_nodes
        exp <- exp_new
        if (exp > max_exp) {
          max_nodes <- nodes
          max_exp <- exp
        }
        break
      }
    }
  }
  return(induced.subgraph(graph, max_nodes))
}


graph <- graph.fb
simulate(graph, iter=15, p_i=0.01, p_h=0.9)
vac <- rev(order(page.rank(graph)$vector))[1:100]
simulate(graph, vac=vac, iter=15, p_i=0.01, p_h=0.9)


# Twitter
simulate(graph.tw, iter=10, inf=0.001, p_i=0.005, p_h=0.3)

# 5 Vaccination Selection

vac_base1 <- rev(order(page.rank(graph.tw)$vector))[1:500]
sim_vac_base1 <- simulate(graph.tw, vac=vac_base1, iter=10, inf=0.001, p_i=0.005, p_h=0.3)


vac_base2 <- sample(1:gorder(graph.tw), 500)
sim_vac_base2 <- simulate(graph.tw, vac=vac_base2, iter=10, inf=0.001, p_i=0.005, p_h=0.3)

vac_base3 <- rev(order(degree(graph.tw)))[1:500]
sim_vac_base3 <- simulate(graph.tw, vac=vac_base3, iter=10, inf=0.001, p_i=0.005, p_h=0.3)

vac_base4 <- rev(order(betweenness(graph.tw)$vector))[1:500]
sim_vac_base4 <- simulate(graph.tw, vac=vac_base4, iter=10, inf=0.001, p_i=0.005, p_h=0.3)

vac_base5 <- rev(order(eigen_centrality(graph.tw)$vector))[1:500]
sim_vac_base5 <- simulate(graph.tw, vac=vac_base5, iter=10, inf=0.001, p_i=0.005, p_h=0.3)


# 4 Samples

sample_mcmc <- sample_expansion_mcmc(graph.tw, 2000)
vac1 <- rev(order(page.rank(sample_mcmc)$vector))[1:500]
simulate(graph.tw, vac=vac1, iter=10, inf=0.001, p_i=0.005, p_h=0.3)

sample_snowball <- sample_expansion_mcmc(graph.tw, 2000)
vac2 <- rev(order(page.rank(sample_snowball)$vector))[1:500]
simulate(graph.tw, vac=vac2, iter=10, inf=0.001, p_i=0.005, p_h=0.3)

sample_random_node <- sample_random_walk_nodes(graph.tw, 2000, mode="ALL")
vac3 <- rev(order(page.rank(sample_random_node)$vector))[1:500]
simulate(graph.tw, vac=vac3, iter=10, inf=0.001, p_i=0.005, p_h=0.3)

sample_random_edge <- sample_random_walk_edges(graph.tw, 5000, mode="ALL")
vac4 <- rev(order(page.rank(sample_random_edge)$vector))[1:500]
simulate(graph.tw, vac=vac4, iter=10, inf=0.001, p_i=0.005, p_h=0.3)

for (i in 0:15) {  
  inf <- append(inf, mean(vaccinated[seq(1,160,by=16)+i]))
  v <- append(v, mean(infected[seq(1,160,by=16)+i]))
}

plot(inf, main = "Immunization ")
points(v, pch=15)
