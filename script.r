setwd("C:\\Users\\Lukas\\OneDrive\\Uni\\Master\\3_HWS19\\Complex and Social Networks\\Project")

library("igraph")

#graph <- read.graph("Simple_pairwise-London_tube_map.txt", format = "ncol", directed = FALSE)

n <- 100
ba <- barabasi.game(n, 1, directed=FALSE)
ws <- watts.strogatz.game(1, n, 3, 0.1)
er <- erdos.renyi.game(n, 0.3)
tree <- make_tree(n, children=2, mode="undirected")

#degree.org <- degree(graph)
#closeness.org <- closeness(graph)
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

plot_history <- function(graph, history) {
  for(infected in history) {
    plot_with_color(graph, infected)
  }
}

simulate <- function(graph, vac=integer(0), iter=10, inf=0.01, p_i=0.05, p_h=0.5) {
  n <- gorder(graph)
  infected <- sample(setdiff(seq(1, n), vac), round(n * inf), replace = FALSE)
  res <- c(length(infected))
  history <- list(infected)
  for (r in seq(iter)) {
    healed <- infected[sample(c(TRUE, FALSE), length(infected), prob = c(p_h, 1-p_h), replace = TRUE)]
    neighbours <- c()
    for (v in infected) {
      nei <- as_ids(neighbors(graph, v))
      nei <- setdiff(nei, append(infected, vac))
      neighbours <- append(neighbours, nei)
    }
    new_infected <- unique(neighbours[sample(c(TRUE, FALSE), length(neighbours), prob = c(p_i, 1-p_i), replace = TRUE)])
    infected <- append(infected[-healed], new_infected)
    history[[length(history)+1]] <- infected
    res <- append(res, c(length(infected)))
  }
  return(history)
}

sample_random_walk <- function(graph, num) {
  sampled <- make_empty_graph(directed = FALSE)
  v <- sample(as_ids(V(graph)), 1)
  v.id <- 1
  sampled <- add_vertices(sampled, 1, name=v)
  for (iter in seq(num)){
    v_n <- sample(as_ids(neighbors(graph, v)), 1)
    v_n.id <- match(v_n, as_ids(V(sampled)))
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

sample_metropolis_hastings <- function(graph, num) {
  sampled <- make_empty_graph(directed = FALSE)
  v <- sample(as_ids(V(graph)), 1)
  v.id <- 1
  v_n <- NA
  sampled <- add_vertices(sampled, 1, name=v)
  for (iter in seq(num)){
    neighbours <- as_ids(neighbors(graph, v))
    repeat {
      v_n <- sample(neighbours, 1)
      if (runif(1) < min(1, degree(graph, v)/degree(graph, v_n))) {
        break
      }
    }
    v_n.id <- match(v_n, as_ids(V(sampled)))
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

sample_expansion <- function(graph, num) {
  
}
  
vac <- rev(order(page.rank(er)$vector))[1:10]
