# collasing_branches_for_Habia_pub

#### PGLS for-loop #### 
library(ape)
library(phylobase)
library(phytools)
library(Rcpp)
library(geiger)
library(caper)
library(MCMCglmm)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(MuMIn)

# remotes::install_github("GuangchuangYu/tidytree")

library(tidytree)
library(ggtree)
########################################################################

setwd("C:/Research/Masters_thesis/MT_tree_analysis")

# Read in Tree and data 
mt_tree <- read.nexus("FullMCC.tree.nexus") 
plot(mt_tree)

uce_tree <- read.nexus("UCE_tree.contree")
plot(uce_tree)


######################### mt _tree ##########################

# Define the targeted branches you want to collapse
mt_targeted_taxa <- c ("Habia_fuscicauda","Habia_atrimaxillaris","Habia_gutturalis","Habia_rubica",
                      "Chlorothraupis_carmioli","Chlorothraupis_olivacea","Chlorothraupis_stolzmanni")

# Define the target branches with a representative for all other clades

outside_taxa <- c( "Cardinalis_cardinalis", "Caryothraustes_poliogaster", "Piranga_bidentata", "Pheucticus_ludovicianus",
                   "Periporphyrus_erythromelas", "Rhodothraupis_celaeno")

targeted_taxa <-    c(mt_targeted_taxa, outside_taxa)             
# collapse all to[s except for target taxa

# Identify the nodes that are not ancestral to the target tips
non_target_nodes <- setdiff(mt_tree$tip.label, targeted_taxa)

# Collapse the non-target nodes in the tree
collapsed_tree <- collapse.singles(ape::drop.tip(mt_tree, non_target_nodes))

# Plot the collapsed tree
plot(collapsed_tree)
write.tree(collapsed_tree, file = "Habia_manuscript_figures/mt_tree_collapse.tree")



####### UCE tree #####

# All taxa for habia paper
#uce_targered_taxa <- c("Chlorothraupis_carmioli_frenata", "Chlorothraupis_carmioli_lutescens", "Chlorothraupis_stolzmanni",
#                       "Chlorothraupis_olivacea", "Habia_rubica_peruviana",
#                       "Habia_fuscicauda_salvini", "Habia_atrimaxillaris","Habia_gutturalis", "Habia_cristata")
# Again with fewer taxa 

uce_targered_taxa <- c("Chlorothraupis_olivacea", "Habia_rubica_peruviana",
                       "Habia_fuscicauda_salvini", "Habia_atrimaxillaris","Habia_gutturalis", "Habia_cristata")


# Define the target branches with a representative for all other clades
outside_taxa <- c( "Cardinalis_phoeniceus", "Caryothraustes_poliogaster", "Piranga_bidentata", "Pheucticus_ludovicianus",
                   "Periporphyrus_erythromelas", "Rhodothraupis_celaeno")

targeted_taxa <-    c(uce_targered_taxa, outside_taxa)             
# collapse all to[s except for target taxa

# Identify the nodes that are not ancestral to the target tips
non_target_nodes <- setdiff(uce_tree$tip.label, targeted_taxa)

# Collapse the non-target nodes in the tree
collapsed_tree <- collapse.singles(ape::drop.tip(uce_tree, non_target_nodes))

# Plot the collapsed tree
plot(collapsed_tree)
write.tree(collapsed_tree, file = "Habia_manuscript_figures/uce_tree_collapse.tree")


########33 Attempts in ggtree ###################

# Create a ggtree object from the collapsed tree
ggtree_obj <- ggtree(tree) # collapsed_tree


p <- ggtree_obj + geom_tiplab() + ggplot2::xlim(0, 0.06)



collapse(p, 51, 'mixed', fill='steelblue', alpha=.4) %>% 
  collapse(65, 'mixed', fill='firebrick', color='blue')

######



# Get the node label of the target tip
target_node <- tree$node.label[tree$tip.label == "Cardinalis_cardinalis"]

df <- p$data

# Get the coordinates of the tip position
tip_pos <- data.frame(x = p %>% get_node_pos(target_node, "y"), y = p %>% get_node_pos(target_node, "x"))

# Define the size of the triangle to be added
triangle_size <- 0.05

# Add a triangle over the target tip
ggtree_obj <- ggtree_obj %<+% geom_polygon(aes(x = c(tip_pos$x - triangle_size,  $x + triangle_size, tip_pos$x),
                                               y = c(tip_pos$y, tip_pos$y, tip_pos$y - triangle_size)), fill = "red", color = NA)

# Plot the tree with the highlighted tip
ggtree_obj













# Create a phylogenetic tree object
ggtree_obj <- ggtree(tree)
plot(tree)

p <- ggtree(tree, aes(color=group)) + theme(legend.position='none') +
  scale_color_manual(values=c("black", "firebrick", "steelblue"))
p

p2 <- p + geom_tiplab()
node <- 21
collapse(p2, node, 'max') %>% expand(node)
collapse(p2, node, 'min') %>% expand(node)
collapse(p2, node, 'mixed') %>% expand(node)
collapse(p, 21, 'mixed', fill='steelblue', alpha=.4) %>% 
  collapse(23, 'mixed', fill='firebrick', color='blue')




###################3

ggtree_obj <- ggtree_obj %>% 
  hide_labels() %>% # hide all labels
  collapse_branch(targeted_taxa, keep_tips = TRUE) %>% # collapse branches
  geom_label(aes(label = ifelse(is.tip, label, ""), fill = node), alpha = 0.7, label.padding = unit(0.05, "lines")) # add labels back


# Display the tree
ggtree_obj



###################

# Create a data frame mapping edge labels to tip labels
edge_to_tip <- data.frame(start = tree$edge[,1],
                          end = tree$edge[,2],
                          tip = rep(tree$tip.label, each = 2),
                          stringsAsFactors = FALSE)
edge_to_tip$tip[is.na(edge_to_tip$tip)] <- tree$node.label[edge_to_tip$start[is.na(edge_to_tip$tip)]]



# Print the matrix
edge_to_tip


# Highlight the targeted branches by creating a triangle over them
ggtree_obj <- ggtree_obj + geom_tiplab(align = TRUE) +
  geom_label2(aes(subset = !is.na(label)), 
              fill = "white", 
              label.padding = unit(0.2, "lines"), 
              size = 3.5, 
              alpha = 0.8, 
              color = "black") +
  geom_nodepoint(size = 0) +
  geom_label2(aes(subset = is.na(label)), 
              fill = "white", 
              label.padding = unit(0.2, "lines"), 
              size = 2.5, 
              alpha = 0.8, 
              color = "black") +
  theme_tree2()
plot(ggtree_obj)
# Connect edge labels to tip labels
for (i in 1:nrow(edge_to_tip)) {
  start_node <- edge_to_tip[i, 1]
  end_node <- edge_to_tip[i, 2]
  start_label <- tree$node.label[start_node]
  end_label <- tree$node.label[end_node]
  ggtree_obj <- ggtree_obj + geom_segment2(aes(x = x, y = y, xend = xend, yend = yend),
                                           data = data.frame(x = start_node, y = 0, xend = end_node, yend = 0),
                                           size = 1, color = "gray") +
    geom_text(aes(x = (start_node + end_node) / 2, y = -0.1, label = end_label),
              data = data.frame(x = start_node, y = 0, xend = end_node, yend = 0),
              size = 3.5, color = "black", check_overlap = TRUE)
}





# Highlight the targeted branches by creating a triangle over them
ggtree_obj %<+% subset(tree, edge %in% targeted_branches) +
  geom_tippoint() + geom_tiplab(align = TRUE) +
  geom_label2(aes(subset = !is.na(label)), 
              fill = "white", 
              label.padding = unit(0.2, "lines"), 
              size = 3.5, 
              alpha = 0.8, 
              color = "black") +
  geom_nodepoint(size = 0) +
  geom_label2(aes(subset = is.na(label)), 
              fill = "white", 
              label.padding = unit(0.2, "lines"), 
              size = 2.5, 
              alpha = 0.8, 
              color = "black") +
  theme_tree2()

