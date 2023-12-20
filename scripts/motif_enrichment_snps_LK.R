library(MPRAnalyze)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggh4x)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(here)

# Originally written by Robin Meyers, I've adapted for my data


################################################################################
# format lk inputs for dynamic and differential analysis
metadata <- read_tsv(here("data/new_lib_5.20_jmprimers.txt"))
metadata <- metadata %>% 
  unite(col = "frag_id", locus, allele, sep = "_") %>%
  dplyr::select(snp = snp_id, frag_id) %>% distinct()

motif_dir <- here("output/motifs/motif_matches/")


# format lk mpranalyze quant input for diff analysis
res_mq <- read_tsv(here("output/mpranalyze_quantitative_results.tsv"))

diff_res <- dplyr::select(res_mq, element_group = oligo_id, 
                             condition = celltype, condition_fdr = fdr)


################################################################################
out_dir <- file.path(here("output/motifs/tfbs_enrichments"))
dir.create(out_dir, recursive = T)

metadata <- metadata %>%
  distinct(snp, frag_id)

tfbs_data <- tibble(
    file = list.files(motif_dir, pattern = ".tsv", full.names = T),
    motif = basename(file) %>% str_remove("\\.tsv"),
    data = map(file, read_tsv, col_types = cols(TF = col_character())))

tfbs_sites <- tfbs_data %>%
    dplyr::select(-file) %>%
    mutate(data = map(data, ~ mutate(., TF = ifelse(is.logical(TF) & as.logical(TF), "", TF)))) %>%
    unnest(data)
# editing here to make one entry per motif/fragment
# not currently using the number of sites per fragment or position of fragment
# so might as well collapse to save memory

tfbs_sites <- tfbs_sites %>% dplyr::select(motif, seqnames) %>%
  distinct()

tfbs_sites <- tfbs_sites %>%
  dplyr::rename(frag_id = seqnames) %>%
    left_join(metadata)


################################################################################
# want df of results with element_group, condition, condition_fdr

# this join throws a warning but I'm ok with it
# diff_res is unique on element_group and condition
# tfbs_sites is unique on motif and frag_id (element_group)
# there should be a motif/condition explosion
differential_snps_tfbs <- diff_res %>%
  left_join(tfbs_sites, by = c("element_group" = "frag_id"))

differential_tfbs_summary <- differential_snps_tfbs %>%
    group_by(condition) %>%
    mutate(total_snps = n_distinct(element_group),
           differential_snps = n_distinct(element_group[condition_fdr < 0.05])) %>%
    group_by(motif, condition, total_snps, differential_snps) %>%
    summarise(motif_snps = n_distinct(element_group),
              motif_differential_snps = n_distinct(element_group[condition_fdr < 0.05])) %>%
    ungroup() %>%
    mutate(fisher_test = pmap(across(c(total_snps, differential_snps, motif_snps, motif_differential_snps)),

                              ~ fisher.test(matrix(c(..1,..2,..3,..4), nrow = 2)) %>% broom::tidy())) %>%
    unnest(fisher_test) %>%
    mutate(log_odds = log(estimate),
           fdr = p.adjust(p.value, method = 'fdr'))

sum(differential_tfbs_summary$fdr < 0.05)/dim(differential_tfbs_summary)[1]

write_tsv(differential_tfbs_summary, file.path(out_dir, "mpraq_differential_tfbs_summary_snpfdr0.05.tsv"))


################################################################################
################################################################################
# Visualize
# make a heatmap for mpranalyze quant results
# one for motifs significant in differential analysis
# one for motifs not significant in differential analysis

# put in the motif names from hocomoco
hocomoco_pwms <- motifbreakR::hocomoco %>% convert_motifs("TFBSTools-PWMatrix") %>%
  set_names(map_chr(., ~ name(.))) %>% do.call(PWMatrixList, .)
h_names <- lapply(hocomoco_pwms, function(x) x@name)
h_id <- lapply(hocomoco_pwms, function(x) x@ID) 
id_conversion = data.frame(hoco_names = unlist(h_names),
                           hoco_ids = unlist(h_id))
# hoco_names <- read_tsv(here("output/hocomoco_id_name_conversion.tsv"))

activity_mq <- dplyr::select(res_mq, condition = celltype, element_group = oligo_id, zscore) %>%
  pivot_wider(names_from = condition, values_from = zscore)


motifs_rename <- left_join(differential_tfbs_summary, id_conversion,
                           by = c("motif" = "hoco_ids"))
motifs_rename <- motifs_rename %>% 
  group_by(hoco_names, motif) %>%
  summarise(sig_conditions = sum(fdr < 0.01),
            min_p = min(fdr))



snp_motif_activity <- left_join(activity_mq, tfbs_sites, by = c("element_group"="frag_id"))

sig_motifs <- unique(differential_tfbs_summary$motif[differential_tfbs_summary$fdr < 0.01])

snp_motif_activity_sig <- snp_motif_activity %>% 
  dplyr::filter(motif %in% sig_motifs)

snp_motif_activity_sig_med <- snp_motif_activity_sig %>%  
  dplyr::select(-element_group, -snp) %>%
  group_by(motif) %>%
  summarise_all(function(x) median(x, na.rm=T))

motifs_choose <- motifs_rename %>%
  dplyr::filter(motif %in% sig_motifs) %>%
  group_by(hoco_names) %>%
  summarise(id_to_choose = ifelse(length(unique(motif)) == 1, unique(motif), 
                                  ifelse(length(unique(sig_conditions)) == 1,
                                         motif[which.min(min_p)],
                                         motif[which.max(sig_conditions)])))

snp_motif_activity_sig_rename <- inner_join(snp_motif_activity_sig_med,
                                           motifs_choose,
                                           by = c("motif"="id_to_choose")) %>%
  dplyr::select(-motif) %>%
  dplyr::rename(motif = hoco_names)
  
snp_motif_activity_sig_rename <- snp_motif_activity_sig_rename %>%
  dplyr::filter(!motif == "GABPB1+GABPB2")

snp_motif_activity_sig_rename <- snp_motif_activity_sig_rename %>%
  dplyr::filter(!is.na(motif)) %>%
  column_to_rownames(var = "motif")


cell_names <- read_tsv(here("data/celltype_name_colors.txt"))

store_cols <- ifelse(cell_names$celltype == colnames(snp_motif_activity_sig_rename),
                     cell_names$display_name,
                    colnames(snp_motif_activity_sig_rename))
colnames(snp_motif_activity_sig_rename) <- store_cols

enr_annot <- differential_tfbs_summary %>%
  dplyr::filter(fdr < 0.01) %>%
  group_by(motif) %>%
  summarise(enrichment = case_when(
    sum(log_odds > 0) == length(condition) ~ "enriched",
    sum(log_odds < 0) == length(condition) ~ "depleted",
    T ~ "mixed")) %>%
  dplyr::filter(!is.na(motif))

enr_annot <- inner_join(enr_annot, motifs_choose,
                                            by = c("motif"="id_to_choose")) %>%
  dplyr::filter(!is.na(motif))

lambert <- read_tsv(here("data/HumanTFs_Lambert2018_1-s2.0-S0092867418301065-mmc2.txt"))

annot_df <- data.frame(motif = rownames(snp_motif_activity_sig_rename))
annot_df <- left_join(annot_df, dplyr::select(lambert, motif = Name,
                                              DBD))
dbd_keep <- names(table(annot_df$DBD))[table(annot_df$DBD) > 3]
annot_df <- annot_df %>% 
  mutate(DBD = case_when(DBD == "Unknown" ~ NA,
                         DBD %in% dbd_keep ~ DBD,   
                         T ~ NA)) 


enr_annot <- left_join(enr_annot, annot_df,
                       by = c("hoco_names" = "motif"))
enr_annot <- column_to_rownames(enr_annot, "hoco_names")
enr_annot <- enr_annot %>%
  dplyr::select(DBD, enrichment)


start = quantile(snp_motif_activity_sig_rename,
                 seq(0, 1, 0.01), na.rm=T)[2]
end = quantile(snp_motif_activity_sig_rename, 
               seq(0, 1, 0.01), na.rm=T)[100]
# myBreaks = c(seq(start, end, 0.01))

# try complexheatmap
col_fun = circlize::colorRamp2(c(start, 0, end), colors = c("cadetblue2", "white", "darkgoldenrod1"))

# new dimensions
enr_annot_ch <- HeatmapAnnotation(df = enr_annot[c(2,1)],
                                  simple_anno_size = unit(1.25, 'mm'),
                                  annotation_label = c("DBD", "Enrichment"),
                                  annotation_name_gp = gpar(fontsize = 5, family = "Helvetica",
                                                            fontface = "bold"),
                                  col = list(enrichment = c(enriched = "black", 
                                                            depleted = "grey", 
                                                            mixed = "grey95"),
                                             DBD = c(Homeodomain = "tomato",
                                                     bHLH = "turquoise",
                                                     bZIP = "yellow",
                                                     `C2H2 ZF` = "mediumseagreen",
                                                     `Nuclear receptor` = "purple3",
                                                     Forkhead = "goldenrod1",
                                                     `Homeodomain; POU` = "pink")),
                                  na_col = "white",
                                  annotation_legend_param = list(title_gp = gpar(fontsize = 5,
                                                                                 family = "Helvetica",
                                                                                 fontface = "bold"),
                                                                 labels_gp = gpar(fontsize = 4,
                                                                                  family = "Helvetica"),
                                                                 grid_height = unit(0.5, 'mm'),
                                                                 grid_width = unit(1, 'mm')))

p <- ComplexHeatmap::Heatmap(t(snp_motif_activity_sig_rename),
                             col = col_fun,
                             cluster_rows = T,
                             cluster_columns = T,
                             heatmap_width = unit(80, 'mm'),
                             heatmap_height = unit(30, 'mm'),
                             show_row_names = T,
                             show_row_dend = F,
                             show_column_dend = F,
                             top_annotation = enr_annot_ch,
                             heatmap_legend_param = list(title = "Activity",
                                                         title_gp = gpar(fontsize = 5,
                                                                         fontface = "bold",
                                                                         family = "Helvetica"),
                                                         labels_gp = gpar(fontsize = 5,
                                                                          family = "Helvetica"),
                                                         legend_height = unit(17, 'mm'),
                                                         grid_width = unit(2, 'mm')),
                             column_names_gp = gpar(fontsize = 1.3, family = "Helvetica"),
                             row_names_gp = gpar(fontsize = 5, family = "Helvetica"))


pdf(here('images/figure_panels/median_zscores_mpraq_snpfdr0.05_motiffdr0.01_annot_outlier01_complexheatmap.pdf'))
draw(p)
dev.off()

# write out supplementary table files

dts_format <- differential_tfbs_summary %>% rowwise() %>%
  mutate(condition = cell_names$display_name[cell_names$celltype == condition])
dts_format <- dts_format %>%
  dplyr::rename(HOCOMOCO_ID = motif)
dts_format <- left_join(dts_format, dplyr::select(id_conversion, motif= hoco_names,
                                                  HOCOMOCO_ID = hoco_ids))
dts_format <- dts_format[c(15, 1:14)]
write_tsv(dts_format,
          file.path(here("output/SupplementaryTables/differential_tfbs_activity_summary.tsv")))

