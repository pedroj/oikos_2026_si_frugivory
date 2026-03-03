
### Analysis of strength

This is code to assess strength variation of trees in the monolayer networks and in the
multiplex network, and compare them.

```{r multilayer_strength, echo= TRUE, fig.height=12,fig.width=14}
# Get a list of igraph objects as layers igraph 
g_layers_bip <- get_igraph(multilayer = multilayer_bip_pru, bipartite = T, directed = F)
g_layers_bip$nodes

# Get strength scores of igraph layers
strength_igraph_herb<- igraph::strength(g_layers_bip$layers_igraph$'Herbivory')
strength_igraph_herb<- as.data.frame(unlist(strength_igraph_herb))
colnames(strength_igraph_herb)<- c("strength")
strength_igraph_herb$node<- rownames(strength_igraph_herb)
rownames(strength_igraph_herb)<- NULL
strength_igraph_herb$layer<- rep("Herbivory", 30)

strength_igraph_poll<- igraph::strength(g_layers_bip$layers_igraph$'Pollination')
strength_igraph_poll<- as.data.frame(unlist(strength_igraph_poll))
colnames(strength_igraph_poll)<- c("strength")
strength_igraph_poll$node<- rownames(strength_igraph_poll)
rownames(strength_igraph_poll)<- NULL
strength_igraph_poll$layer<- rep("Pollination", 62)

strength_igraph_disp<- igraph::strength(g_layers_bip$layers_igraph$'Seed dispersal')
strength_igraph_disp<- as.data.frame(unlist(strength_igraph_disp))
colnames(strength_igraph_disp)<- c("strength")
strength_igraph_disp$node<- rownames(strength_igraph_disp)
rownames(strength_igraph_disp)<- NULL
strength_igraph_disp$layer<- rep("Seed dispersal", 40)

strength <- bind_rows(strength_igraph_herb, strength_igraph_poll, strength_igraph_disp)

# Save the mono-layer strength values
strength_igraph_mono<- list(strength_igraph_herb, strength_igraph_poll,
                            strength_igraph_disp)

# Multilayer strength
# Get strength scores
# Strength centrality
strength_sam.g <- igraph::strength(graph_sam.g) # Using the SAM matrix. 264 values.
#strength_SAM<- as.numeric(strength_sam.g[cond])
strength_SAM<- as.data.frame(strength_sam.g)
colnames(strength_SAM)<- c("strength")


# Extract multilayer strength values
cond<- c(1:132) 
strength_SAM<- as.data.frame(strength_SAM$strength[cond]) # I'm extracting half of the  
# dataset, as it is duplicated  
# in the SAM array.

# Extract the nodes 
# Extract the trees
# Extract the trees
cond<- c(12:30,74:92,114:132) # Tree strengths in the monolayer networks.
strength_SAM_trees<- as.data.frame(strength_SAM[cond])

# merge to data frame eigenvector comparison
eigen_igraph<- ec_comparison[,2]
eigen_multi_trees<- eigen_igraph[1:19]

ec_comparison <- data.frame(eigen_sam.g, eigen_igraph) ## CHECK!!!
colnames(ec_comparison) <- c("multilayer","monolayer")

# Scatter plot for eigenvector centralities comparison
ggplot(data = ec_comparison, aes(monolayer, multilayer)) + 
    geom_point(color = "#00325B", size = 3) + 
    labs(title = "Eigenvector Centrality (EC) comparison", 
         x = "EC (monolayer)", y = "EC (multilayer)")+
    geom_abline()+
    scale_x_continuous(limits = c(0,1))+
    scale_y_continuous(limits = c(0,1))+
    coord_fixed() +
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=15),
          axis.text = element_text(color='black',size = 10),
          legend.text =  element_text(size=15),
          legend.title = element_text(size=20))


############
# Graph stacked bar values of mono- and multilayer eigenvector centralities and strength.
# 
levels(as.factor(eigen_centr$layer))
eigen_centr$layer<- factor(eigen_centr$layer, levels=c('Multilayer', 
                                                       'Herbivory',
                                                       'Pollination', 
                                                       'Seed dispersal'))

ggplot(eigen_centr, aes(x= reorder(node, -eigen.centrality), 
                        y= eigen.centrality, fill=layer))+
    geom_bar(stat="identity", width=1, col="black")+
    theme_bw()+
    scale_x_discrete(name="")+
    coord_flip()

```

```{r multi_stats, echo=TRUE, fig.height=9, fig.width=12}
# Read the multilayer stats.
multi_stats<- read.table(here::here("data/raw/multiplex/muxViz_multistats.csv"), sep=",", dec=".", header= TRUE)


```
