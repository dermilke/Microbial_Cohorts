### Verification of network clusters via Modularity significance and stability analysis ###
#
# This function requires the following packages:
# tidyverse, igraph, mclust
#
# This function has two functionalities: 
# 1) it tests if the modularity of the inferred cluster-membership yields significantly higher modularity
# scores than by random chance
# 2) it tests if the network clusters are stable in the light of node-removal (network attack) and
# goes one step further: checks whether this stability is higher than random expectations
#
# Random networks are created by keeping the node-degree distribution and rewiring the existing network based on
# niter_rewire iterations (default = 2500).
#
# We infer modularity of random networks using the same total number of clusters as in the observed
# dataset, as we believe this is a fair comparison here. To do that, we use the cluster_fluid_communities()
# function from the igraph package since this approach allows us to specify the total number of clusters.
#
# Network attacks are implemented via node removal (remove_nodes = TRUE) but can also be realized
# by randomly deleting edges (remove_nodes = FALSE). However, we highly recommend to use node-removal
# which can be biologically interpreted as extinction events.
#
# The function requires a network object as from the igraph package. And it requires a cluster-table
# that is the output of our clustering pipeline.
# 
# Author: Felix Milke (felix.milke@uni-oldenburg.de)
# Date: 26.02.2025

verify_network_cluster <- function(network, cluster_tab, res_param, bootstrap_num, remove_nodes = T,
                                   random_expectation = F, niter_rewire = 2500, attack_bootstrap = 50,
                                   dataset_name) {
  
  ## 1) Significance of cluster modularity
  
  # Total number of cluster
  k <- max(cluster_tab$Cluster)
  
  # Infer significance of cluster modularity by comparing against random networks.
  modularity_random <- replicate(bootstrap_num, {
    rewired_network <- rewire(network, with = keeping_degseq(niter = niter_rewire)) %>%
      induced_subgraph(., which(components(.)$membership == 1))
    #modularity(rewired_network, cluster_fluid_communities(rewired_network, no.of.communities = k)$membership)
    modularity(rewired_network, cluster_fast_greedy(rewired_network)$membership)
  })
  
  # Compute observed modularity based on the clustering results from our clustering pipeline
  #observed_modularity <- modularity(network, left_join(tibble::enframe(V(network), name = "ASV"), cluster_tab, by = "ASV")$Cluster)
  observed_modularity <- modularity(network, cluster_fast_greedy(network)$membership)
  
  # Output the significance of our observed modularity
  print(paste0(dataset_name, ": P-value: ", mean(modularity_random >= observed_modularity), 
               " (total of ", sum(modularity_random < observed_modularity), " out of ", 
               bootstrap_num, " iterations below observation) - ", 
               "Observed: ", round(observed_modularity, digits = 4), " - Avg. random: ", round(mean(modularity_random), digits = 4), 
               " - SD random: ", round(sd(modularity_random), digits = 4)))
  
  # Visualize the distribution of modularity values from random networks and highlight the observed modularity
  p1 <- tibble(Modularity = modularity_random, bootstrap_ID = seq(1, bootstrap_num)) %>%
    ggplot(aes(x = modularity_random)) +
    geom_density(fill = "steelblue") +
    geom_vline(xintercept = observed_modularity, col = "forestgreen", lwd = 1.5, lty = 1) +
    theme_bw() +
    labs(x = "Modularity", y = "Density") +
    lims(x = c(0, 1))
  
  ## 2) Cluster stability against network perturbation
  
  # Use inferred cluster-membership from clustering-pipeline as observed values
  cluster_obs_network <- cluster_tab$Cluster
  
  # If attack-mode = node-removal
  if (remove_nodes == T) {
    
    # Bootstrap network attack with varying attack-strength (keep between 20 and 100% of all nodes)
    attack_test <- purrr::map(seq(0.2, 1, by = 0.05), function(attack_strength) {
      # This bootstrap-iteration calculates observed cluster-stability
      result_tmp_obs <- replicate(bootstrap_num, {
        network_subsample <- subgraph(network, sample(V(network), size = attack_strength * vcount(network), replace = F))
        mclust::adjustedRandIndex(cluster_obs_network[match(names(V(network_subsample)), names(V(network)))], 
                                  cluster_leiden(network_subsample, resolution_parameter = res_param, objective_function = "modularity")$membership)
      })
      
      # Comparison with random-network stability
      if (random_expectation) {
        
        # This bootstrap-iteration calculates cluster-stability of randomized networks
        result_tmp_rewired <- replicate(attack_bootstrap, {
          
          # First create a random network
          network_rewired <- rewire(network, with = keeping_degseq(niter = niter_rewire)) %>%
            induced_subgraph(., which(components(.)$membership == 1))
          
          # Cluster the random network
          cluster_rewired <- cluster_leiden(network_rewired, resolution_parameter = res_param, objective_function = "modularity")$membership
          
          # Assess cluster-stability of randomized network under current attack_strength as the average
          # ARI score of 10 random attacks
          replicate(10, {
            network_subsample <- subgraph(network_rewired, sample(V(network_rewired), size = attack_strength * vcount(network_rewired), replace = F))
            mclust::adjustedRandIndex(cluster_rewired[match(names(V(network_subsample)), names(V(network_rewired)))], 
                                      cluster_leiden(network_subsample, resolution_parameter = res_param, objective_function = "modularity")$membership)
          }) %>%
            mean()
        })
        
        tibble(Attack_Strength = attack_strength, ARI = c(result_tmp_obs, result_tmp_rewired), 
               Type = c(rep("Observed", length(result_tmp_obs)), rep("Rewired", length(result_tmp_rewired))))
        
      # Without comparison with random networks
      } else {
        tibble(Attack_Strength = attack_strength, ARI = result_tmp_obs, Type = "Observed")
      }
      
    }) %>%
      dplyr::bind_rows()
    
  # Attack edges instead of nodes
  } else if (remove_nodes == F) {
    
    # As above: Calculate stability of clusters when attacked (edge-removal)
    attack_test <- purrr::map(seq(0.2, 1, by = 0.05), function(attack_strength) {
      result_tmp_obs <- replicate(bootstrap_num, {
        network_subsample <- subgraph.edges(network, sample(E(network), size = attack_strength * ecount(network), replace = F), 
                                            delete.vertices = F)
        mclust::adjustedRandIndex(cluster_obs_network, 
                                  cluster_leiden(network_subsample, resolution_parameter = res_param, objective_function = "modularity")$membership)
      })
      
      # Comparison with random-network stability
      if (random_expectation) {
        
        # This bootstrap-iteration calculates cluster-stability of randomized networks
        result_tmp_rewired <- replicate(attack_bootstrap, {
          
          # First create a random network
          network_rewired <- rewire(network, with = keeping_degseq(niter = niter_rewire)) %>%
            induced_subgraph(., which(components(.)$membership == 1))
          
          # Cluster the random network
          cluster_rewired <- cluster_leiden(network_rewired, resolution_parameter = res_param, objective_function = "modularity")$membership
          
          # Assess cluster-stability of randomized network under current attack_strength as the average
          # ARI score of 10 random attacks
          replicate(10, {
            network_subsample <- subgraph.edges(network_rewired, sample(E(network_rewired), size = attack_strength * ecount(network_rewired), replace = F), 
                                                delete.vertices = F)  # Resample edges
            mclust::adjustedRandIndex(cluster_rewired, 
                                      cluster_leiden(network_subsample, resolution_parameter = res_param, objective_function = "modularity")$membership)
          }) %>%
            mean()
        })
        tibble(Attack_Strength = attack_strength, ARI = c(result_tmp_obs, result_tmp_rewired), 
               Type = c(rep("Observed", length(result_tmp_obs)), rep("Rewired", length(result_tmp_rewired))))
        
      # Without comparison with random networks
      } else {
        tibble(Attack_Strength = attack_strength, ARI = result_tmp_obs, Type = "Observed")
      }
      
    }) %>%
      dplyr::bind_rows()
    
  }
  
  # Visualize cluster stability as boxplots for both, observed and randomized networks
  p2 <- ggplot(attack_test, aes(x = as.factor(round(1-Attack_Strength, digits = 3)), y = ARI, fill = Type)) +
    geom_boxplot() +
    labs(x = paste0("Attack Strength (Proportion of ", ifelse(remove_nodes, "nodes", "edges"), " removed)"), y = "ARI", 
         title = paste0("Network-attack simulation by removing ", ifelse(remove_nodes, "nodes", "edges"))) +
    theme_bw()
  
  # Return both plots, cluster-significance-table, pvalue of clusters (as proportion of random iterations that
  # yielded equal or higher modularity scores as the observed one), and cluster-attack-table.
  return(list(p1 = p1, p2 = p2, Cluster_Significance_Table = tibble(Modularity_randomized = modularity_random), 
              pval = mean(modularity_random >= observed_modularity), Cluster_Attack_Table = attack_test))
  
}