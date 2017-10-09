library(igraph)
library(readr)
library(DTMCPack)
NUMBER_OF_NODES <- 5000
NUMBER_OF_EDGES <- 10000
NUMBER_OF_GRAPHS <- 20
DIRECTORY_NAME = '~/edges_datasets/'
MAXIMUM_NUMBER_OF_EDGES <- (NUMBER_OF_NODES * (NUMBER_OF_NODES - 1)) / 2

random_walk_entropy_rate <- function(graph){
  degrees <- degree(graph)
  stationary_distribution <- degrees / (2 * ecount(graph))
  return (log2(2*ecount(graph)) - entropy(stationary_distribution))
}
random_walk_stationary_distribution <- function(graph){
  return (degree(graph) / (2 * ecount(graph)))
}
random_walk_entropy <- function(graph){
  return (entropy(random_walk_stationary_distribution(graph)))
}
entropy_of_page_rank <- function(graph){
  save <- page.rank(graph)$vector
  return (entropy(save))
}

distances_page_rank <- function(graphs){
  distances_page <- c()
  distances <- c()
  stepi <- 0
  stepj <- 0
  for (i in graphs){
    print('Step i')
    print(stepi)
    step <- stepi + 1
    stepj <- 0
    for (j in graphs){
      stepj <- stepj +1
      print('stepj')
      print(stepj)
      #stat_i <- random_walk_stationary_distribution(i)
      #stat_j <- random_walk_stationary_distribution(j)
      #middle <- (stat_i + stat_j) / 2
      stat_i_page <- page.rank(i)$vector
      stat_j_page <- page.rank(j)$vector
      middle_page <- (stat_i_page + stat_j_page) / 2
      distances_page <- c(distances_page, jensen_shannon_divergence(stat_i_page, stat_j_page, middle_page))
      #distances <- c(distances, jensen_shannon_divergence(stat_i, stat_j, middle))
    }
  }
  #save_distances <<- matrix(distances, nrow = length(graphs), ncol = length(graphs))
  return (matrix(distances_page, nrow = length(graphs), ncol = length(graphs)))
}
graphs_from_multiplex_dataset <- function(dataset){
  graphs <- list()
  min_node <- min(dataset$from, dataset$to)
  max_node <- max(dataset$from, dataset$to)
  for (i in unique(dataset$layer)){
    new_set <- dataset[dataset$layer == i, ]
    new_set$layer <- NULL
    vertex_id <- min_node:max_node
    vertex_frame <- as.data.frame(vertex_id)
    graphs[[i]] <- simplify(graph_from_data_frame(d = new_set, directed = FALSE, vertices = vertex_frame))
  }
  return (graphs)
}

read_dataset_as_multiplex_network <- function(path){
  dataset <<- read_delim(path, " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  attr_names <- c('layer', 'from', 'to', 'weight')
  colnames(dataset) <- attr_names
  graphs <- graphs_from_multiplex_dataset(dataset)
  return (graphs)
}
layer_set <- data.frame(Name = character(), num_layer = integer())
analyze_datasets <- function(){
  names <- list.files(DIRECTORY_NAME)
  for (i in names){
    print(i)
    name = i
    path = paste0(DIRECTORY_NAME, name)
    graphs <- read_dataset_as_multiplex_network(path)
    g <- data.frame(name, length(graphs))
    layer_set <<- rbind(layer_set, g)
    #test <<- graphs
    #dists <- (distances_page_rank(graphs))
    #heatmap(dists, Rowv = NA, Colv = NA, hclustfun = function(x) hclust(x, method = 'ward.D2'))
    #dists <- as.dist(dists)
    #cluster_graphs_page_rank(dists, graphs, name)
  }
}
disjoint_random_graphs <- function(number_of_disjoint_graphs){
  snapshotter <- list()
  number_of_nodes <- 10000
  number_of_edges <- as.integer(0.05 * number_of_nodes)
  g <- create_random_graph(number_of_nodes = number_of_nodes, num_of_edges = number_of_edges)
  result <- list()
  for (i in 1:number_of_disjoint_graphs){
    result[[i]] <- g
    g <- rewire(g, each_edge(0.6))
  }
  return (result)
} 
analyze_random_graphs <-function(){
  graphs <<- create_multiple_random_graphs(NUMBER_OF_NODES, NUMBER_OF_EDGES)
  dists <- distances_page_rank(graphs)
  save_dists <<- dists
  dists <- as.dist(dists)
  cluster_graphs_page_rank(dists, graphs, 'name')
}
create_regular_graphs <- function(num_of_layers){
  result <- list()
  for (i in 1:num_of_layers){
    result[[i]] <- sample_k_regular(NUMBER_OF_NODES, i)
  }
  return (result)
}
analyze_graph <- function(graphs){
  save_graphs <<- graphs
  dists <- (distances_page_rank(graphs))
  save_disjoin_dists <<- dists
  dists <- as.dist(dists)
  cluster_graphs_page_rank(dists, graphs, 'name')
}
#results <- disjoint_random_graphs(20)
#dists <- distances_page_rank(results)
#heatmap(dists, Rowv = NA, hclustfun = function(x) hclust(x, method = 'ward.D2'))
#graphs <- create_multiple_random_graphs(NUMBER_OF_NODES, NUMBER_OF_EDGES)
#dists <- distances_page_rank(graphs)
#new_dists <- as.dist(dists)
#disjoint_graphs <- disjoint_random_graphs(number_of_disjoint_graphs = NUMBER_OF_GRAPHS)
#analyze_graph(disjoint_graphs)
analyze_datasets()