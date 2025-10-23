#proteomics final figure plot
#helen pennington
#october 7, 2025

library(readxl)
library(dplyr)
library(data.table)
library(igraph)
library(ggraph)
library(viridis)

# read file (update filename/sheet as needed)
res <- read_excel("/restricted/projectnb/cteseq/projects/somascan/data/fgsea_grouped_all.xlsx", sheet = "Groups")
res_filt <- res[,1:2]
res_filt$Group[res_filt$Group == "other"] <- "Other"
#functions
#merge with groups
add_groups <- function(res_filt, fgsea){
  res_group <- merge(fgsea, res_filt, by.x = "pathway", by.y = "Pathways",all.x = TRUE)
  return(res_group)
}

#make sure all significant pathways have a group
check_groups <- function(fgsea_group){
  fgsea_sig <- fgsea_group[which(fgsea_group$padj < 0.05),]
  print(which(is.na(fgsea_sig$Group)))
  return(fgsea_sig)
}

##results files##
#first line reads in file and adds group names from res_filt
#second line (and maybe third) makes sure that all significant pathways have groups and if not adds them in the next line
#third or fourth line adds model column to designate what model the pathway result came from to help when they are all in one file

#RHI vs Low CTE fgsea results
rhivlow <- add_groups(res_filt, fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_RHIvslow.csv"))
rhivlow_sig <- check_groups(rhivlow)
rhivlow_sig[, model := "rhivlow"]

#Low CTE vs High CTE fgsea results
lowvhigh <- add_groups(res_filt,fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_lowvshigh.csv"))
lowvhigh_sig <- check_groups(lowvhigh)
lowvhigh_sig[, model := "lowvhigh"]

#AT8 Total fgsea results
AT8 <- add_groups(res_filt,fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_AT8.csv"))
AT8_sig <- check_groups(AT8)
AT8_sig[, model := "AT8_total"]

#Total Years of Play fgsea results
totyrs <- add_groups(res_filt,fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_totyrs.csv"))
totyrs_sig <- check_groups(totyrs)
totyrs_sig$Group[62] <- "GTPase"
totyrs_sig[, model := "totyrs"]

#CDS total fgsea results
cds <- add_groups(res_filt,fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_cds.csv"))
cds_sig <- check_groups(cds)
cds_sig[, model := "cds"]

#Dementia fgsea results
dementia <- add_groups(res_filt,fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_dementia.csv"))
dementia_sig <- check_groups(dementia)
dementia_sig$Group[c(100, 189, 190, 205, 208, 209)] <- "GTPase"
dementia_sig[, model := "dementia"]

#FAQ Total fgsea results
faq <- add_groups(res_filt,fread("/restricted/projectnb/cteseq/projects/somascan/results/fgsea/gseafiles/plot_format_files/fgsea_compact_faq.csv"))
faq_sig <- check_groups(faq)
faq_sig$Group[c(148, 158)] <- "GTPase"
faq_sig[, model := "faq"]

# Combine all model fgsea results files
combined <- rbindlist(list(rhivlow_sig,lowvhigh_sig,AT8_sig,totyrs_sig,cds_sig,dementia_sig,faq_sig), use.names = TRUE, fill = TRUE)
write.csv(combined, "/restricted/projectnb/cteseq/projects/somascan/final_plots/combined_fgsea_groups.csv")

# summarize counts
summary_counts <- combined %>%
  mutate(direction = ifelse(NES > 0, "Up", "Down")) %>%
  group_by(model, Group, direction) %>%
  summarise(
    n_pathways = n(),
    .groups = "drop"
  ) %>%
  group_by(model, Group) %>%
  mutate(
    total = sum(n_pathways),
    frac = n_pathways / total
  ) %>%
  ungroup()

summary_counts

edges <- summary_counts %>%
  select(from = model, to = Group, direction, weight = n_pathways)

# build graph
g <- graph_from_data_frame(edges, directed = TRUE)
V(g)$type <- V(g)$name %in% unique(edges$from)
# plot
ggraph(g, layout = "bipartite") +
  geom_edge_link(aes(color = direction, width = weight), alpha = 0.8) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_edge_color_manual(values = c("Up" = "red", "Down" = "blue")) +
  theme_void()


edges <- summary_counts %>%
  select(from = model, to = Group, direction, weight = n_pathways)

# Build graph
g <- graph_from_data_frame(edges, directed = TRUE)

# Tag node types for bipartite layout
V(g)$type <- V(g)$name %in% unique(edges$from)

scale_edge_color_gradientn(colours = okabe_ito)
# Plot
ggraph(g, layout = "bipartite") +
  geom_edge_link(aes(edge_color = weight, edge_width = weight), alpha = 0.8) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_color_gradientn(colours = rainbow(7)) +
  theme_void()



g <- g %>%
  igraph::set_edge_attr("signed_weight",
                        value = ifelse(E(g)$direction == "Up", E(g)$weight, -E(g)$weight)
  )

# plot
ggraph(g, layout = "bipartite") +
  geom_edge_link(
    aes(edge_color = signed_weight, edge_width = weight),
    alpha = 0.8
  ) +
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  # Diverging blue–gray–red scale (colorblind-safe)
  scale_edge_color_gradient2(
    low = "#4575B4",      # blue for down
    mid = "gray95",       # neutral
    high = "#D73027",     # red for up
    midpoint = 0
  ) +
  scale_edge_width(range = c(0.5, 3)) +
  theme_void() +
  theme(legend.position = "right")