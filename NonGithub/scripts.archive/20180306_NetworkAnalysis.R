## R code written by Elizabeth A. Bowman Mar. 6, 2018
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to do a network analysis of EM community betwen ranges and sites with 
## differing burn history

## Created using Jesse Sadler's blog post on network analyses
## https://www.jessesadler.com/post/network-analysis-with-r/

#========================================================================================#
# Load data and libraries----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#
install.packages('igraph')
install.packages('tidyverse')
library(ggplot2);library (vegan); library(indicspecies)
library(igraph); library(tidyverse); library(network)

#----------------------------------------------------------------------------------------#
# set up paths to directories----
#----------------------------------------------------------------------------------------#
#--path to directory where the climate data downloaded from BIOCLIM site is stored
dat.dir <- "~/Documents/PhD/2_EM_Fire_effect/data/"
fig.dir <- '~/Documents/PhD/2_EM_Fire_effect/figures_output/'
res.dir <- "~/Documents/PhD/2_EM_Fire_effect/results_output/"

#----------------------------------------------------------------------------------------#
# Load data and clean up: Tree level----
#----------------------------------------------------------------------------------------#
#--Load data
all.data <- read.csv(paste0(dat.dir,'20170806_OTU_data.csv'), as.is = T)

#--Remove sequences not assigned to an OTU
otu.data <- all.data[!is.na(all.data$otu.97), ]

#--Remove sequences not assigned to Ponderosa host
otu.data <- otu.data[otu.data$Host == 'Ponderosa',]

#--Add an OTU count column (1 sequence = 1 frequency count)
otu.data['otu.count'] <- 1

#--Creates matrix, grouping OTUs by tree number and getting a sum of EM tip abundance 
# for each OTU
otu.data %>%
  select(Sample_name,Burn_status,Range,Site,Tree,otu.97,otu.count) %>%
  spread(otu.97,otu.count) %>%
  group_by(Tree,Site,Range,Burn_status) %>%
  select(matches ("otu*")) %>%
  summarize_each(funs(sum(., na.rm = TRUE))) %>%
  as.data.frame() -> stsp.matrix

#--Creates matrix, grouping OTUs by tree number and getting a sum of EM tip abundance 
# for each OTU
otu.data %>%
  select(Sample_name,Burn_status,Range,Site,Tree,otu.97,Tip_count) %>%
  spread(otu.97,Tip_count) %>%
  group_by(Tree,Site,Range,Burn_status) %>%
  select(matches ("otu*")) %>%
  summarize_each(funs(sum(., na.rm = TRUE))) %>%
  as.data.frame() -> abund.matrix

#========================================================================================#
# Network analysis----
#========================================================================================#

#----------------------------------------------------------------------------------------#
# OTU count data----
#----------------------------------------------------------------------------------------#

#--Create node list
species <- otu.data %>%
  distinct(otu.97) %>%
  rename(label = otu.97)

sites <- otu.data %>%
  distinct(Site) %>%
  rename(label = Site)

## Join both sets of nodes
nodes <- full_join(species, sites, by = "label") %>% 
  rowid_to_column("id")
nodes

#--create edge list with node lists and add weight (count occurrence)
per_site <- otu.data %>%
  group_by(otu.97, Site) %>%
  summarise(weight = n()) %>%
  ungroup()
per_site

#--Link the species and site data to the unique ids in the nodes object
edges.orig <- per_site %>%
  left_join(nodes, by = c('otu.97' = 'label')) %>%
  rename(from = id)

edges.orig <- edges %>%
  left_join(nodes, by = c('Site' = 'label')) %>%
  rename(to = id)

edges <- select(edges.orig, from, to, weight)
edges

#<< Network it >> --------
# em.network <- network(edges, vertex.attr = nodes,
#                       matrix.type = 'edgelist',
#                       ignore.eval = F)
# em.network
# 
# library(network)
# plot(em.network, vertex.cex = 3)
# plot(em.network, vertex.cex = 3, mode = 'circle')
# 
# detach(package:network)
# rm(routes_network)
# library(igraph)
# 
# routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
# plot(routes_igraph, edge.arrow.size = 0.2)
# plot(routes_igraph, layout = layout_with_graphopt, edge.arrow.size = 0.2)

library(tidygraph)
library(ggraph)

routes_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
routes_igraph_tidy <- as_tbl_graph(routes_igraph)

routes_tidy %>% 
  activate(edges) %>% 
  arrange(desc(weight))
ggraph(routes_tidy) + geom_edge_link() + geom_node_point() + theme_graph()
ggraph(routes_tidy, layout = "graphopt") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Letters") +
  theme_graph()

install.packages('visNetwork')
install.packages('networkD3')
library(visNetwork)
library(networkD3)
visNetwork(nodes, edges)

edges.1 <- mutate(edges, width = weight/5 + 1)

visNetwork(nodes, edges.1) %>% 
  visIgraphLayout(layout = "layout_with_fr") %>% 
  visEdges(arrows = "middle")


nodes_d3 <- mutate(nodes, id = id - 1)
edges_d3 <- mutate(edges, from = from - 1, to = to - 1)
forceNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
             NodeID = "label", Group = "id", Value = "weight", 
             opacity = 1, fontSize = 16, zoom = TRUE)

sankeyNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
              NodeID = "label", Value = "weight", fontSize = 16, unit = "Letter(s)")


#----------------------------------------------------------------------------------------#
# Root tip data----
#----------------------------------------------------------------------------------------#

#--Create node list
species <- otu.data %>%
  distinct(otu.97) %>%
  rename(label = otu.97)

sites <- otu.data %>%
  distinct(Site) %>%
  rename(label = Site)

## Join both sets of nodes
nodes <- full_join(species, sites, by = "label") %>% 
  rowid_to_column("id")
nodes

#--create edge list with node lists and add weight (count occurrence)
per_site <- otu.data %>%
  group_by(otu.97, Site) %>%
  summarise(weight = sum(Tip_count)) %>%
  ungroup()
per_site

#--Link the species and site data to the unique ids in the nodes object
edges.orig <- per_site %>%
  left_join(nodes, by = c('otu.97' = 'label')) %>%
  rename(from = id)

edges.orig <- edges.orig %>%
  left_join(nodes, by = c('Site' = 'label')) %>%
  rename(to = id)

edges <- select(edges.orig, from, to, weight)
edges

#<< Network it >> --------
# em.network <- network(edges, vertex.attr = nodes,
#                       matrix.type = 'edgelist',
#                       ignore.eval = F)
# em.network
# 
# library(network)
# plot(em.network, vertex.cex = 3)
# plot(em.network, vertex.cex = 3, mode = 'circle')

# detach(package:network)
# rm(routes_network)
# library(igraph)
# 
# routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
# plot(routes_igraph, edge.arrow.size = 0.2)
# plot(routes_igraph, layout = layout_with_graphopt, edge.arrow.size = 0.2)

library(tidygraph)
library(ggraph)

routes_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
routes_igraph_tidy <- as_tbl_graph(routes_igraph)

routes_tidy %>% 
  activate(edges) %>% 
  arrange(desc(weight))
ggraph(routes_tidy) + geom_edge_link() + geom_node_point() + theme_graph()
ggraph(routes_tidy, layout = "graphopt") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  labs(edge_width = "Letters") +
  theme_graph()

library(visNetwork)
library(networkD3)
visNetwork(nodes, edges)

edges.1 <- mutate(edges, width = weight/5 + 1)

visNetwork(nodes, edges.1) %>% 
  visIgraphLayout(layout = "layout_with_fr") %>% 
  visEdges(arrows = "middle")


nodes_d3 <- mutate(nodes, id = id - 1)
edges_d3 <- mutate(edges, from = from - 1, to = to - 1)
forceNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
             NodeID = "label", Group = "id", Value = "weight", 
             opacity = 1, fontSize = 16, zoom = TRUE)

sankeyNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
              NodeID = "label", Value = "weight", fontSize = 16, unit = "Letter(s)")

