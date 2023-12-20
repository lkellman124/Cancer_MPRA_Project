
library(igraph)
library(visNetwork)
library(ggraph)
library(ggrepel)
library(RColorBrewer)
library(ggsci)
library(tidyverse)
library(cowplot)

# written by Robin Meyers, I only slightly adapted to change formatting
# plots coessential network (but does not perform coessentiality analysis)

set.seed(42)
plot_network <- function(network) {
    color_levels <- unique(V(network)$Annotation) %>% {.[!is.na(.)]}

    annotation_colors <- pal_d3(palette = "category20")(length(color_levels)) %>%
        set_names(color_levels)

    g_layout <- create_layout(network, "fr")

    ggraph(g_layout) +
        geom_edge_link(color = "grey50") +
        geom_node_point(aes(size = ifelse(is.na(Annotation), 1, 2.5)),
                        color = "grey50") +
        geom_node_point(aes(color = factor(Annotation, levels = names(annotation_colors)),
                            size = ifelse(is.na(Annotation), 1, 2.5))) +
        geom_text_repel(aes(x = x, y = y,
                            label = ifelse(is.na(Annotation), name, name),
                            size = ifelse(is.na(Annotation), 2.5, 4)),
                        max.overlaps = 1000,
                        # fontface = "bold",
                        nudge_x = 0.1,
                        nudge_y = -0.1,
                        point.padding = unit(.2, "lines"),
                        box.padding = unit(.25, "lines"),
                        min.segment.length = unit(0, "lines"),
                        show.legend = F) +
        scale_size_identity() +
        scale_color_manual(values = annotation_colors, na.translate = F) +
        scale_x_continuous(expand = c(0.1, 0)) +
        scale_y_continuous(expand = c(0.1, 0)) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme_nothing() +
        theme(legend.position = "bottom", legend.title = element_blank())
}

coessential_network <- read_rds(here("data/eGene_coessential-ppi_network.rds"))
# After some analysis changes, SMU1 dropped out of the ptGene set
# Getting rid of it manually
coess_edit <- coessential_network
coess_edit <- delete_vertices(coess_edit, "SMU1")

name_key <- read_tsv(here("data/celltype_name_colors.txt"))
name_key <- dplyr::filter(name_key, !celltype == "lx")
name_key <- name_key %>% dplyr::select(-full_name, -cancer)
name_key <- rbind(name_key, data.frame(celltype = c("multiple", "none"),
                                       display_name = c("Multiple", "None"),
                                       color = c("black", "grey")))
annot_fix <- V(coess_edit)$Annotation
annot_fix <- sapply(annot_fix,
                    function(x) ifelse(is.na(x), "None",
                                       name_key$display_name[name_key$celltype == x]))

color_vals <- name_key$color
names(color_vals) <- name_key$display_name
color_vals_2 <- color_vals
color_vals_2[["None"]] <- NA
V(coess_edit)$Annotation <- factor(annot_fix, levels = name_key$display_name,
                                   ordered = T)
  
plot_network_onlyhitnames <- function(network, color_vals) {
    g_layout <- create_layout(network, "fr")
    
    ggraph(g_layout) +
        geom_edge_link(color = "grey60", edge_width = 0.2) +
        geom_node_point(aes(size = ifelse(Annotation == "None", 0.4, 0)),
                        color = "grey50") +
        geom_node_point(aes(color = Annotation,
                            size = ifelse(Annotation == "None", 0.8, 0.8))) +
        geom_text_repel(aes(x = x, y = y,
                            label = ifelse(is.na(Annotation), "", name),
                            size = ifelse(Annotation == "None", 0, 4)),
                        max.overlaps = 1000,
                        # fontface = "bold",
                        nudge_x = 0.1,
                        nudge_y = -0.1,
                        point.padding = unit(.2, "lines"),
                        box.padding = unit(.25, "lines"),
                        min.segment.length = unit(0.4, "lines"),
                        show.legend = F) +
        scale_size_identity() +
      #  scale_color_manual(values = color_vals) +
        scale_x_continuous(expand = c(0.1, 0)) +
        scale_y_continuous(expand = c(0.1, 0)) +
        guides(color = guide_legend(override.aes = list(size = 0.1))) +
      scale_color_manual(values = color_vals_2) +
        theme_nothing() +
        theme(legend.position = "bottom", legend.title = element_blank(),
           #   legend.key.size = unit(0.2, 'cm'),
              text = element_text(size = 8))
           #   text = element_text(size = 5, family = "Helvetica"))
}

plot_network_onlyhitnames(coess_edit, color_vals)

set.seed(42)
pdf(file=here("images/figure_panels/coessentiality_network.pdf"), width = 3.7, height = 7)
plot_network_onlyhitnames(coess_edit, color_vals)
dev.off()



