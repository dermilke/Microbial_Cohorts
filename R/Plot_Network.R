plot_network <- function(network, parameter, layout_network = NULL, vertex_size_multiplier = 1, color_pal = NULL, reverse = F, n = 9) {
  
  get_colors_cont <- function(vec, palette, reverse = F, n = 9) {
    
    vec <- ifelse(is.na(vec), mean(vec, na.rm = T), vec)
    
    if (reverse) { ramp <- rev(RColorBrewer::brewer.pal(n, palette)) %>% colorRamp(.) } else 
    { ramp <- RColorBrewer::brewer.pal(n, palette) %>%  colorRamp(.) } 
    
    ramp(vegan::decostand(vec, method = "range")) %>%
      rgb(., maxColorValue = 255)
  }
  
  if (is.null(layout)) {
    layout_network <- igraph::layout_nicely(network)
  }
  
  vertex_size <- igraph::get.vertex.attribute(network, "Plot_Size")*vertex_size_multiplier
  vertex_cluster <- tibble(Cluster = igraph::get.vertex.attribute(network, "Cluster") %>% as.character())
  
  if (is.null(color_pal) & parameter == "Cluster") {
    vertex_color <- igraph::get.vertex.attribute(network, "Color")
  } else if (length(color_pal) == length(igraph::V(network))) {
    vertex_color <- color_pal
  } else {
    vertex_color <- get_colors_cont(igraph::get.vertex.attribute(network, parameter), color_pal, reverse = reverse, n = n)
  }
  
  plot(network, layout = layout_network, vertex.color = vertex_color, vertex.size = vertex_size,
       main = parameter)
  
}