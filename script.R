library(igraph)
library(ggplot2)

final_data <- data.frame(Name = character(), 
                         Violation_SubAddivity = integer(),
                         Average_Information_loss = double(), 
                         Information_loss_layer = character())
create_random_graph <- function(number_of_nodes, num_of_edges){
  graph <- erdos.renyi.game(number_of_nodes, num_of_edges, type = 'gnm', directed = TRUE)
  return (graph)
}
transition_row <- function(row){
  if (sum(row) == 0)
    return (rep(1/length(row), length(row)))
  return(row/sum(row))
}
markov_transition_matrix <- function(graph){
  adj <- as.matrix(as_adjacency_matrix(graph))
  adj <- apply(adj, 1, transition_row)
  alpha <- 0.85
  len <- nrow(adj) * ncol(adj)
  I_N <- matrix(1/ncol(adj), nrow= nrow(adj), ncol=ncol(adj))
  return (adj * alpha + (1 - alpha)*I_N)
}
jensen_shannon_divergence <- function(first, second, middle){
  entropy_first <- entropy(first)
  entropy_second <- entropy(second)
  entropy_middle <- entropy(middle)
  return (entropy_middle - (entropy_first + entropy_second) / 2)
}
entropy_graph <- function(graph){
  g <- laplacian_scaled(graph)
  distr <- eigen(g)$values
  return (entropy(distr))
}
entropy <- function(distribution){
  return (-sum(distribution * log2(distribution)))
}
laplacian_scaled <- function(graph){
  laplacian <- as.matrix(laplacian_matrix(graph))
  result <- laplacian / sum(diag(laplacian))
  return (result)
}
distances_eigen <- function(graphs){
  distances <- c()
  for (i in graphs){
    for (j in graphs){
      first <- eigen(laplacian_scaled(i))$values
      second <- eigen(laplacian_scaled(j))$values
      middle <- eigen((laplacian_scaled(i) + laplacian_scaled(j))/2)$values
      distances <- c(distances, jensen_shannon_divergence(first, second, middle))
    }
  }
  return (matrix(distances, nrow = length(graphs), ncol = length(graphs)))
}
create_multiple_random_graphs <- function(num_of_nodes, num_of_edges){
  result <- list()
  
  for (i in 1:NUMBER_OF_GRAPHS){
    current_edge <- num_of_edges
    result[i] <- list(create_random_graph(num_of_nodes, current_edge))
  }
  return (result)
}
eigen_values <- function(matrix){
  return (eigen(matrix)$values)
}
entropy_of_network_per_layer <- function(graphs){
  laplacians <- lapply(graphs, laplacian_scaled)
  eigens <- lapply(laplacians, eigen_values)
  entropies <- lapply(eigens, entropy)
  return (mean(unlist(entropies)))
}
aggregated_entropy_rate <- function(graphs){
  agg <- do.call(union, graphs)
  return (entropy_of_page_rank(agg))
}
aggregated_eigen_entropy <- function(graphs){
  agg <- do.call(union, graphs)
  return (entropy_graph(agg))
}



cluster_graphs_page_rank <- function(distances, graphs, name){
  
  distances <- as.dist(distances)
  clusters <- hclust(distances, "ward.D2")
  plot(clusters)
  mergers <- clusters$merge
  agg_entropy <- aggregated_entropy_rate(graphs)
  remember_merges <- list()
  edge_difference <- c()
  entropy_rates_per_layer <- c()
  temp_graphs <- graphs
  total_edges <- 0
  for (i in 1:(length(graphs))){
    total_edges <- total_edges + ecount(graphs[[i]])
  }
  for (i in 1:(length(graphs) - 1)){
    merge_indices <- mergers[i,]
    print('clustering first')
    print(merge_indices)
    if (merge_indices[1] < 0){
      index <- merge_indices[1] * (-1)
      first_merge <- graphs[[index]]
      
    }else {
      index <- merge_indices[1]
      first_merge <- remember_merges[[index]]
    }
    if (merge_indices[2] < 0){
      index <- merge_indices[2] * -1
      second_merge <- graphs[[index]]
    }else {
      index <- merge_indices[2]
      second_merge <- remember_merges[[index]]
    }
    first_adj <- as_adjacency_matrix(first_merge)
    second_adj <- as_adjacency_matrix(second_merge)
    sum_adj <- first_adj + second_adj
    union_graphs <- graph_from_adjacency_matrix(sum_adj, mode='undirected')
    union_graphs <- simplify(union_graphs)
    edge_difference <- c(edge_difference,
                         (ecount(union_graphs) - max(ecount(first_merge), ecount(second_merge))) /
                           min(ecount(first_merge), ecount(second_merge)))
    remember_merges[[i]] <- union_graphs
  }
  information_loss_per_layer <- edge_difference
  entropy_rates <- lapply(graphs, entropy_of_page_rank)
  entropy_rates <- c(unlist(entropy_rates))
  entropy_rates_per_layer <- c(entropy_rates_per_layer, mean(entropy_rates))
  print(entropy_rates)
  for (i in seq(1, length(graphs) - 1, 1)){
    acquire_indices <- c(mergers[1:i,])
    negative <- acquire_indices[acquire_indices < 0]
    positives <- acquire_indices[acquire_indices > 0]
    f_positives <- 1:i
    if (length(positives) > 0){
      f_positives <- f_positives[-positives]
    }
    negative_graphs <- graphs[negative]
    positive_graphs <- remember_merges[f_positives]
    graphs_final <- c(negative_graphs, positive_graphs)
    count_edges <- 0
    for (i in 1:length(graphs_final)){
      count_edges <- count_edges + ecount(graphs_final[[i]])
    }
    print(i)
    total_edges <- c(total_edges, count_edges)
    entropy_rates <- lapply(graphs_final, entropy_of_page_rank)
    entropy_rates <- c(unlist(entropy_rates))
    entropy_rates_per_layer <- c(entropy_rates_per_layer, mean(entropy_rates))
  }
  ilrp <- ""
  for (i in information_loss_per_layer){
    ilrp <- paste0(ilrp, ', ', round(i, 2)) 
  }
  save_total <<- total_edges
  #path = paste0("~/save_plots/",name,".png")
  #path_quality = paste0("~/save_plots/", name,"_quality_function.png")
  #png(filename = path)
  #rint('entropy_rates_per_layer')
  #print(entropy_rates_per_layer)
  #print('Quality function is')
  number_of_nodes <- length(V(graphs[[1]]))
  cutting_point <- 0.01 * log2(number_of_nodes)
  entropy_plot <- entropy_rates_per_layer
  save_entropy <<- entropy_plot
  save_edge <<- edge_difference
  entropy_rates_per_layer <- diff(entropy_rates_per_layer, 1)
  quality_function <- (entropy_rates_per_layer / log2(number_of_nodes))
  final_decision <- cumsum(quality_function)
  violation <- diff(entropy_rates_per_layer, 1)
  numViol <- length(violation < 0)
  #print('Cumsum final')
  #print(final_decision)
  decision <- head(which(final_decision > cutting_point), 1)
  #print(decision)
  #print('Agg entropy')
  #print(agg_entropy)
  add_row <- data.frame(name, numViol, round(mean(information_loss_per_layer), 2), ilrp)
  final_data <<- rbind(final_data, add_row)
  plot(entropy_plot, xlab = 'Layers aggregated', ylab='Entropy', type='l')
  plot(edge_difference, xlab = 'Layers aggregated', ylab='Information loss', type='l')
  plot(total_edges, xlab = 'Layers aggregated', ylab='Total number of edges', type='l')
}


