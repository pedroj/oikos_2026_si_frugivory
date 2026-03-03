my_plot_modular_matrix<-
function (x, fix_coordinates = T, axes_titles = c("Set 1", "Set 2"), 
          transpose = F, outside_module_col = "gray") 
{
    if (class(x) != "infomap_monolayer") {
        stop("x must be of class infomap_monolayer")
    }
    M_set1 <- M_set2 <- x$edge_list[1:3]
    names(M_set1) <- names(M_set2) <- names(x$edge_list)[1:3] <- 
        c("Set1", "Set2", "w")
    suppressMessages(suppressWarnings(M_set1 %<>% dplyr::left_join(x$modules, 
                                      by = c(Set1 = "node_name")) %>% 
                                      rename(module1 = module_level1)))
    suppressMessages(suppressWarnings(M_set2 %<>% dplyr::left_join(x$modules, 
                                      by = c(Set2 = "node_name")) %>% 
                                      dplyr::rename(module2 = module_level1)))
    suppressMessages(suppressWarnings(M <- dplyr::full_join(M_set1, 
                                           M_set2, by = c("Set1", "Set2", "w")) %>% 
                                           dplyr::select(Set1, Set2, w, module1, module2)))
    Set1_modules <- unique(M_set1[, c("Set1", "module1")])
    Set1_modules <- with(Set1_modules, Set1_modules[order(module1, Set1), ])
    Set2_modules <- unique(M_set2[, c("Set2", "module2")])
    Set2_modules <- with(Set2_modules, Set2_modules[order(module2, Set2), ])
    M %<>% mutate(edge_in_out = 
                      ifelse(module1 == module2, "in", 
                             "out")) %>% mutate(value_mod = 
                                         ifelse(edge_in_out == 
                                                "in", module1, 0)) %>%
                                         mutate(Set1 = factor(Set1, 
                                                       levels = Set1_modules$Set1), 
                                                       Set2 = factor(Set2, 
                                                       levels = Set2_modules$Set2))
    module_colors <- tibble(module1 = unique(M$module1), 
                            col = gg_color_hue(n = length(unique(M$module1))))
    suppressMessages(M %<>% dplyr::left_join(module_colors) %>% 
                         mutate(col = ifelse(edge_in_out == "in", col, outside_module_col)))
    if (transpose) {
        p <- ggplot() + geom_tile(data = M %>% dplyr::filter(w != 0), 
                                  aes(Set2, Set1, fill = col))
    }
    else {
        p <- ggplot() + geom_tile(data = M %>% dplyr::filter(w != 0), 
                                  aes(Set1, Set2, fill = col))
    }
    p <- p + labs(x = axes_titles[2], y = axes_titles[1]) + scale_fill_identity() + 
        theme(legend.position = "none", panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), panel.background = element_blank(), 
              axis.text.x = element_text(angle = 90), axis.ticks = element_blank())
    if (fix_coordinates) {
        p <- p + coord_fixed()
    }
    return(p)
}
