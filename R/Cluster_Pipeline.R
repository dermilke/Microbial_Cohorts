#### Pipeline ####
# This pipeline requires the following packages:
# tidyverse, cluster, mclust, vegan, igraph
#
# This pipeline clusters a network in a robust and reproducible way. It builds on 
# the Leiden-Clustering algorithm and conducts a consensus-clustering approach within a user-defined interval 
# of clustering resolutions. The algorithm has two modes: 
# find_parameters = TRUE -> User needs to define optimal clustering interval and methodology
# find_parameters = FALSE -> requires res_param, cluster_method and cluster_num != NULL
# This will automatically run the pipeline based on the input parameters.
#
# The cluster_method must be either "Consensus_hierarchical" or "Consensus_network"
# 
# The algorithm conducts the following steps (find_parameters = TRUE):
# 1. build a network based on user-provided edgelist-table
# 2. cluster network with all clustering-resolutions within the interval of 0.1 and 5 in 0.1 steps
# 3. evaluate clustering-performance within this interval by computing four different metrics
# (penalized modularity, shilouette-width, cluster-similarity between adjacent resolution-values, 
# cluster-size distribution)
# 4. The user then selects the optimal interval of clustering-resolution values based on these metrics
# 5. Recluster the network for every clustering resolution within the selected interval (in steps of 0.05)
# 6. Create a co-occurrence consensus-matrix (count the number of shared clusters for all ASV-pairs) and 
# use it as adjacency-matrix for consensus-network
# 7. Cluster the consensus network in 2 ways: 1) network clustering 2) hierarchical clustering
# 8. Ask the user to select best clustering approach
# 9. Return the final network, clustering-table, clustering resolution interval, number of clusters.

complete_cluster_pipeline <- function(network_edgelist_table, 
                                      find_parameters = TRUE, vertex_num_drop = 5,
                                      res_param = NULL, cluster_method = NULL, cluster_num = NULL) {
  
  set.seed(123)
  
  # Create a weighted, undirected network using the provided edgelist-table. 
  # Normalize edge-weights to the range of 0-1.
  # Drop small, unconnected networks (subnetwork < vertex_num_drop).
  network <- network_edgelist_table %>%
    filter(weight > 0 & weight < 1) %>% 
    mutate(weight = vegan::decostand(weight, method = "range")) %>%
    select(From, To, weight) %>%
    igraph::graph_from_data_frame(directed = F) %>%
    igraph::simplify(edge.attr.comb = "max") %>%
    igraph::induced_subgraph(., which(igraph::components(.)$membership %in% 
                                which(table(igraph::components(.)$membership) > vertex_num_drop)))
  
  adjacency_mat <- as.matrix(igraph::as_adjacency_matrix(network, attr = "weight"))
  
  if (find_parameters == TRUE) {
    
    # To find the optimal clustering-resolution interval, we first need to evaluate all clustering-resolution results
    res_test_vec <- seq(0.1, 5, 0.1)
    
    cluster_list_res <- purrr::map(res_test_vec, function(res) {
      igraph::cluster_leiden(network, objective_function = "modularity", resolution_parameter = res, 
                             weights = igraph::E(network)$weight)
    })
    
    ## Compute four different metrics to evaluate clustering results:
    
    # Modularity penalized by cluster-size standard-deviation (otherwise it will always yield smallest resolution)
    sd_res <- purrr::map_dbl(seq_along(res_test_vec), function(ind) {
      clusters <- cluster_list_res[[ind]]
      sd <- table(igraph::membership(clusters)) %>%
        sd()
      return(sd)
    })
    scores <- purrr::map2_dbl(seq_along(res_test_vec), sd_res/max(sd_res, na.rm = T), function(ind, sd_val) {
      clusters <- cluster_list_res[[ind]]
      mod <- igraph::modularity(network, clusters$membership)
      mod - 1 * sd_val 
    })
    cluster_num_highmod <- cluster_list_res[[which.max(scores)]]$membership %>%
      max()
    
    # Shilouette-width: The resolution with the highest silhouette width balances cluster cohesion and separation.
    sil_scores <- purrr::map_dbl(seq_along(res_test_vec), function(ind) {
      clusters <- cluster_list_res[[ind]]$membership
      if (length(unique(clusters)) > 1) {
        mean(cluster::silhouette(clusters, as.dist(1 - adjacency_mat))[, 3])
      } else {NA}
    })
    cluster_num_shilwid <- cluster_list_res[[which.max(sil_scores)]]$membership %>% 
      max()
    
    # Compare cluster results and look for "stable" solutions, where adjacent resolutions produce similar results
    stability_scores <- purrr::map_dbl(2:length(cluster_list_res), function(ind) {
      mclust::adjustedRandIndex(cluster_list_res[[ind-1]]$membership, 
                                cluster_list_res[[ind]]$membership)
    })
    cluster_num_stabscore <- cluster_list_res[[which.max(stability_scores)]]$membership %>% 
      max()
    
    # Compare cluster-size distribution between varying resolutions
    imbalances <- purrr::map_dbl(seq_along(res_test_vec), function(ind) {
      tmp <- table(cluster_list_res[[ind]]$membership)
      sd(tmp) / mean(tmp)
    })
    cluster_num_imbal <- cluster_list_res[[which.min(imbalances)]]$membership %>%
      max()
    
    # Visualize clustering metrics along complete range of clustering resolutions
    par(mfrow = c(2,2), mar = c(5,5,4,1))
    
      plot(res_test_vec, scores, type = "b", xlab = "Resolution", ylab = "Penalized modularity",
           main = paste0("Modularity\nBest resolution: ", res_test_vec[which.max(scores)],
                         " - Cluster: ", cluster_num_highmod))
      points(res_test_vec[which.max(scores)], scores[which.max(scores)], bg = "red",
             pch = 21, col = "black", cex = 1.5)
      
      plot(res_test_vec, sil_scores, type = "b", xlab = "Resolution", ylab = "Shilouette-width",
           main = paste0("Shilouette width\nBest resolution: ", res_test_vec[which.max(sil_scores)],
                         " - Cluster: ", cluster_num_shilwid))
      points(res_test_vec[which.max(sil_scores)], sil_scores[which.max(sil_scores)], bg = "red",
             pch = 21, col = "black", cex = 1.5)
      
      plot(res_test_vec[-1], stability_scores, type = "b", xlab = "Resolution", ylab = "ARI",
           main = paste0("Cluster member stability\nBest resolution: ", res_test_vec[which.max(stability_scores)+1],
                         " - Cluster: ", cluster_num_stabscore))
      points(res_test_vec[which.max(stability_scores)+1], stability_scores[which.max(stability_scores)], bg = "red",
             pch = 21, col = "black", cex = 1.5)
      
      plot(res_test_vec, imbalances, type = "b", xlab = "Resolution", ylab = "Coefficient of variation\nof cluster-size",
           main = paste0("Cluster-size distribution\nBest resolution: ", res_test_vec[which.min(imbalances)],
                         " - Cluster: ", cluster_num_imbal))
      points(res_test_vec[which.min(imbalances)], imbalances[which.min(imbalances)], bg = "red",
             pch = 21, col = "black", cex = 1.5)
    
    # Ask user to select the optimal range of clustering-resolutions based on the above plot.
    # The algorithm proposes an optimal range based on maxima and minima of the computed curves, 
    # but this rarely yields optimal results. Hence, user needs to judge plots.
    res_range <- range(c(res_test_vec[which.max(scores)], res_test_vec[which.max(sil_scores)],
                         res_test_vec[which.max(stability_scores)+1], res_test_vec[which.min(imbalances)]))
    
    message(paste0("Based on these analyses the proposed resolution of Leiden-clustering is in the range of ",
                   res_range[1], " and ", res_range[2]))
    
    res_range[1] <- readline(prompt = "Enter new range based on figures (lower boundary): " ) %>%
      as.numeric()
    
    res_range[2] <- readline(prompt = "(upper boundary): " ) %>%
      as.numeric()
    
    # This is the new interval of clustering-resolution values for the consensus-analysis
    resolution_vec <- seq(res_range[1], res_range[2], by = 0.05)
    
    # The user should input a guess for the expected number of clusters. This is only a first guess
    # and can later be changed iteratively.
    optimal_k_sil <- readline(prompt = "Enter expected number of clusters: " ) %>%
      as.integer()
    
    # Recluster the network only within the interval of clustering-resolutions
    memberships <- purrr::map(resolution_vec, function(x) {
      igraph::cluster_leiden(network, objective_function = "modularity", 
                             resolution_parameter = x, weights = igraph::E(network)$weight)$membership
    }) %>%
      bind_cols() %>%
      magrittr::set_colnames(paste0("res_", resolution_vec))
    
    # Calculate the co-occurrence consensus-matrix from the clustering result
    co_occurrence <- matrix(0, nrow = igraph::vcount(network), ncol = igraph::vcount(network)) %>%
      magrittr::set_rownames(names(igraph::V(network))) %>%
      magrittr::set_colnames(names(igraph::V(network)))
    
    membership_matrix <- as.matrix(memberships)
    
    # Vectorized co-occurrence calculation for turbo-speed using outer matrix-product
    for (col_ind in seq_len(ncol(membership_matrix))) {
      data_tmp <- membership_matrix[, col_ind]
      co_occurrence <- co_occurrence + outer(data_tmp, data_tmp, FUN = "==")
    }
    
    # Normalize by the number of resolutions
    co_occurrence <- co_occurrence / length(resolution_vec)
    
    # Cluster consensus matrix using hierarchical clustering for user-defined guess of cluster number
    hc <- hclust(as.dist(1 - co_occurrence), method = "ward.D2")
    
    satisfied = FALSE
    
    # This loop will be run for as long as the user changes the number of clusters 
    # (it will recluster the consensus-matrix each time)
    while (satisfied == FALSE) {
      
      # Cut the tree to expected number of clusters
      consensus_clusters <- cutree(hc, k = optimal_k_sil)
      
      # Compare consensus clustering to original resolutions
      # (this is for the user to see if his clustering result is very sensitive to the clustering resolution)
      stability_scores <- sapply(memberships, function(clusters) {
        mclust::adjustedRandIndex(consensus_clusters, clusters)
      })
      
      # This plot shows both, the comparison between consensus solution based on
      # user defined cluster-number and the original cluster-results
      #
      # And the hierarchical dendogram that yielded these clusters. 
      # This plot is only an intermediate product for the user to see how sensitive his 
      # consensus clustering result is to the clustering-resolution. (less sensitive = better)
      par(mfrow = c(1,2))
      # Visualize the consensus clustering
        plot(hc, labels = FALSE, main = "Consensus Clustering", xlab = "")
        rect.hclust(hc, k = optimal_k_sil, border = "red")
        
        # Visualize stability across resolutions
        plot(resolution_vec, stability_scores, type = "b", xlab = "Resolution", ylab = "ARI with Consensus",
             main = "Agreement with clusters\nof varying resolution")
      
      void <- readline(prompt = paste0("Now showing the derived consensus-clusters. ", 
                                       "Look whether hierarchical clusters are reasonable\nand agreement with ",
                                       "Leiden-Clusters is high. Continue with Enter."))
      
      # Create networks using the consensus-matrix as adjacency-matrix. 
      # Similar to original networks, we drop very small isolated networks that are likely spurious.
      # Cluster the network using the walktrap-algorithm. This yields the second
      # consensus-clustering approach and will be compared to the hierarchical clustering approach.
      network_consensus <- igraph::graph_from_adjacency_matrix(co_occurrence, mode = "undirected", weighted = TRUE) %>%
        igraph::simplify() %>%
        igraph::set_vertex_attr(name = "Consensus_hierarchical", value = consensus_clusters) %>%
        igraph::induced_subgraph(., which(igraph::components(.)$membership %in% 
                                    which(table(igraph::components(.)$membership) > vertex_num_drop))) %>%
        igraph::set_vertex_attr(name = "Consensus_network", value = igraph::membership(igraph::cluster_walktrap(.)))
      
      network_clusters <- igraph::vertex_attr(network_consensus, "Consensus_network")
      
      cluster_tab_consensus <- tibble(ASV = names(consensus_clusters), 
                                      Cluster = as.integer(consensus_clusters))
      cluster_tab_network <- tibble(ASV = names(igraph::V(network_consensus)), 
                                    Cluster = as.integer(network_clusters))
      
      # Additionally, the original co-occurrence network will be visualized and its nodes
      # coloured according to the two different consensus-clustering approaches.
      network_raw <- network_edgelist_table %>%
        select(From, To, weight) %>%
        igraph::graph_from_data_frame(directed = F) %>%
        igraph::simplify(edge.attr.comb = "max") %>%
        igraph::set_vertex_attr(., name = "Consensus_hierarchical", value = left_join(tibble(ASV = names(igraph::V(.))), cluster_tab_consensus, by = "ASV")$Cluster) %>%
        igraph::set_vertex_attr(., name = "Consensus_network", value = left_join(tibble(ASV = names(igraph::V(.))), cluster_tab_network, by = "ASV")$Cluster) %>%
        igraph::induced_subgraph(., which(igraph::components(.)$membership %in% 
                                    which(table(igraph::components(.)$membership) > vertex_num_drop))) %>%
        igraph::subgraph(vids = cluster_tab_network$ASV)
      
      # Use Fruchterman-Reingold layout
      layout_raw <- igraph::layout_with_fr(network_raw, grid = "nogrid")
      layout_consensus <- igraph::layout_with_fr(network_consensus, grid = "nogrid")
      
      # Here, we visualize all four networks:
      # top 2 networks = consensus networks 
      # (left = clustering based on hierarchical clusters, right = clustering based on consensus-network clusters)
      #
      # bottom 2 networks = observed networks
      # Colors in the observed networks are according to the clustering results from the two consensus-clustering approaches
      par(mfrow = c(2,2))
        igraph::plot.igraph(network_consensus, layout = layout_consensus, 
                            vertex.color = RColorBrewer::brewer.pal(n = 12, name = 'Paired')[igraph::vertex_attr(network_consensus, "Consensus_hierarchical")],
                            vertex.label = NA, vertex.size = 5, main = "Consensus-network - Hierarchical")
        igraph::plot.igraph(network_consensus, layout = layout_consensus, 
                            vertex.color = RColorBrewer::brewer.pal(n = 12, name = 'Paired')[igraph::vertex_attr(network_consensus, "Consensus_network")],
                            vertex.label = NA, vertex.size = 5, main = paste0("Consensus-network - Network clustering\nFound cluster: ", max(cluster_tab_network$Cluster)))
        
        igraph::plot.igraph(network_raw, layout = layout_raw, 
                            vertex.color = RColorBrewer::brewer.pal(n = 12, name = 'Paired')[igraph::vertex_attr(network_raw, "Consensus_hierarchical")],
                            vertex.label = NA, vertex.size = 5, main = "Raw network with hierarchical clusters")
        igraph::plot.igraph(network_raw, layout = layout_raw, 
                            vertex.color = RColorBrewer::brewer.pal(n = 12, name = 'Paired')[igraph::vertex_attr(network_raw, "Consensus_network")],
                            vertex.label = NA, vertex.size = 5, main = "Raw network with network clusters")
        
      # Output in the console the distribution of cluster-sizes for both consensus-approaches
      message(paste0("Number of Cluster-members for Hierarchical-Cluster ", cluster_tab_consensus$Cluster %>% table() %>% 
                       paste0(names(.), " = ", .), "\n"))
      message("")
      message(paste0("Number of Cluster-members for Network-Cluster ", cluster_tab_network$Cluster %>% table() %>% 
                       paste0(names(.), " = ", .), "\n"))
      
      cluster_method <- NA
      
      # Ask user to either enter new number of clusters to re-do the hierarchical clustering, or continue
      optimal_k_sil <- readline(prompt = "Enter different number of clusters if desired (if satisfied leave empty): ")
      
      if (!grepl(pattern = "[0-9]*", x = optimal_k_sil) | optimal_k_sil == "") {satisfied = T} else {optimal_k_sil = as.integer(optimal_k_sil)}
      
    # If continue: you will leave the while-loop that iterates over the consensus-clustering
    # and which uses each iteration a new user-defined number of clusters.
    }
    
    # Ask user for the preferred consensus-approach
    while (!cluster_method %in% c("1", "2")) {
      
      cluster_method <- readline(prompt = "Choose clustering-method (1 - Hierarchical clustering ; 2 - Network clustering): ")
      
      if (cluster_method == 1) {
        cluster_tab <- cluster_tab_consensus
      } else if (cluster_method == 2) {
        cluster_tab <- cluster_tab_network
      } 
    }
    
  # Done! The next else-condition is the case for
  # find_parameters = FALSE 
  # and basically runs the predefined consensus-clustering approach only for the predefined parameters
  } else {
    
    # User-defined values
    res_range = res_param
    optimal_k_sil = cluster_num
    
    resolution_vec <- seq(res_param[1], res_param[2], by = 0.05)
    
    # Run clustering within this interval
    memberships <- purrr::map(resolution_vec, function(x) {
      igraph::cluster_leiden(network, objective_function = "modularity", 
                             resolution_parameter = x, weights = igraph::E(network)$weight)$membership
    }) %>%
      bind_cols() %>%
      magrittr::set_colnames(paste0("res_", resolution_vec))
    
    # Build consensus co-occurrence-matrix
    co_occurrence <- matrix(0, nrow = igraph::vcount(network), ncol = igraph::vcount(network)) %>%
      magrittr::set_rownames(names(igraph::V(network))) %>%
      magrittr::set_colnames(names(igraph::V(network)))
    
    # Vectorized co-occurrence calculation
    membership_matrix <- as.matrix(memberships)
    for (col_ind in seq_len(ncol(membership_matrix))) {
      data_tmp <- membership_matrix[, col_ind]
      co_occurrence <- co_occurrence + outer(data_tmp, data_tmp, FUN = "==")
    }
    
    # Normalize by the number of resolutions
    co_occurrence <- co_occurrence / length(resolution_vec)
    
    # Choose consensus-cluster approach 
    if (cluster_method == "Consensus_hierarchical") {
      
      consensus_clusters <- hclust(as.dist(1 - co_occurrence), method = "ward.D2") %>%
        cutree(k = optimal_k_sil)
      
      cluster_tab <- tibble(ASV = names(consensus_clusters), 
                            Cluster = as.integer(consensus_clusters))
      
    } else if (cluster_method == "Cluster_network") {
      
      network_consensus <- igraph::graph_from_adjacency_matrix(co_occurrence, mode = "undirected", weighted = TRUE) %>%
        igraph::simplify() %>%
        igraph::induced_subgraph(., which(igraph::components(.)$membership %in% 
                                    which(table(igraph::components(.)$membership) > vertex_num_drop))) %>%
        igraph::set_vertex_attr(name = "Consensus_network", value = igraph::membership(igraph::cluster_walktrap(.)))
      
      network_clusters <- igraph::vertex_attr(network_consensus, "Consensus_network")
      
      cluster_tab <- tibble(ASV = names(igraph::V(network_consensus)), 
                            Cluster = as.integer(network_clusters))
      
      optimal_k_sil <- max(cluster_tab$Cluster)
      
    }
    
    network_raw <- network_edgelist_table %>%
      select(From, To, weight) %>%
      igraph::graph_from_data_frame(directed = F) %>%
      igraph::simplify(edge.attr.comb = "max") %>%
      igraph::set_vertex_attr(., name = "Cluster", value = left_join(tibble(ASV = names(igraph::V(.))), cluster_tab, by = "ASV")$Cluster) %>%
      igraph::induced_subgraph(., which(igraph::components(.)$membership %in% 
                                  which(table(igraph::components(.)$membership) > vertex_num_drop))) %>%
      igraph::subgraph(vids = cluster_tab$ASV)
    
  }
  
  # Return the results
  return(list(Cluster_Tab = cluster_tab, Network = network_raw, Resolutionrange = res_range,
              Clusternumber = optimal_k_sil))
  
}
