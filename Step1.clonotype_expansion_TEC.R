###############################################################################
# Import dependencies
###############################################################################
# Seurat ecosystem for scRNA-seq analysis
library(Seurat)
# Utility libraries for clustering and manipulation
library(bluster)
library(tidyverse)
library(qs)          # Fast serialization
library(ggplot2)
library(readxl)
library(purrr)
library(patchwork)
library(kableExtra)
library(data.table)
# Output directories
output_dir <- "~/TCR_project/"
#load TEC T cell object
T_cells_filt2<-qread("~/Immune_project/data_to_share_CART_TEC/T_cells_nov25_deid.qs")

remove <- unique(rownames(T_cells_filt2@meta.data[is.na(T_cells_filt2@meta.data$clonotypeID), ]))
keep<-setdiff(unique(rownames(T_cells_filt2@meta.data)), remove)
T_cells_filt_clono<-subset(T_cells_filt2, cells = keep)
T_cells_filt_clono@meta.data<-T_cells_filt_clono@meta.data%>%mutate(combined_clonotype = paste(TRA_cdr3_aa, TRB_cdr3_aa, sep="_"),
                                                                    new_clones=paste0(PT_ID,"_",combined_clonotype)) 
#take only Diseased, TEC, PB
T_cells_filt_clono<-subset(T_cells_filt_clono, Therapy=="TEC"&Tissue=="PB")
meta<-T_cells_filt_clono@meta.data
meta<-meta%>%select(c("PT_ID", "combined_clonotype", "Timepoint"))
meta<-meta%>%mutate(sample=paste0(PT_ID, "_", Timepoint))
total<-meta%>%group_by(sample)%>%summarise(total=n())
sub<-meta%>%group_by(PT_ID,sample, Timepoint, combined_clonotype)%>%summarise(n=n())

paired<-sub%>%group_by(PT_ID)%>%summarise(n=n_distinct(Timepoint))%>%filter(n>1)%>%pull(PT_ID)
sub<-sub%>%filter(PT_ID%in%paired)

patients<-unique(sub$PT_ID)
res_list <- list()


# Poisson test with sample sizes
poisson_test_with_size <- function(count1, count2, size1, size2) {
  result <- poisson.test(c(count1, count2), T = c(size1, size2))
  return(result$p.value)
}

# Get results over 2 lists
poisson_testing = function(counts1, counts2){
  
  size1=sum(counts1)
  size2=sum(counts2)
  
  p_values <- mapply(poisson_test_with_size, counts1, counts2, 
                     MoreArgs = list(size1 = size1, size2 = size2))
  
  as_tibble(data.frame(
    n1 = counts1,
    n2 = counts2,
    p1 = counts1 / size1,
    p2 = counts2 / size2,
    p = p_values
  )) %>% mutate(fdr = p.adjust(p, method = 'BH'), sig = fdr < 0.05)
}


for(pt in patients){
try<-sub%>%ungroup()%>%filter(PT_ID==pt)%>%select(-c("sample"))
sub_wide<-try%>%pivot_wider(names_from = Timepoint, values_from = n, values_fill = 0)
sub_wide<-sub_wide%>%arrange(combined_clonotype)

pre<-unlist(sub_wide$Pre)
post<-unlist(sub_wide$Post)

total<-meta%>%filter(PT_ID==pt)%>%group_by(Timepoint)%>%summarise(n=n())
tot_pre<-total$n[2]
tot_post<-total$n[1]

res<-poisson_testing(pre, post)%>%
  mutate(combined_clonotype = sub_wide$combined_clonotype)%>% arrange(fdr)
res_list[[pt]] <- res

}

res_df <- bind_rows(res_list, .id = "PT_ID")

#extract meta
output_dir_car<-"~/Immune_project/Paper_Figures_June/"
date_of_sampling <- read_xlsx(paste0(output_dir_car, "Swimmer_Plot_TEC.xlsx"))
# Load the dictionary used for de-identifying samples
deidentified <- read_csv("~/Immune_project/data_to_share_CART_TEC/Dictionary_for_deidentifying_seurat_object.csv")

# Keep only relevant therapy groups
deidentified <- deidentified %>%
  filter(Therapy %in% c("Healthy", "TEC", "TALQ"))

# Select only the columns required for matching and de-identification
#   Sample_ID = original library/sample name
#   ID / PT_ID = patient identifiers
#   Final_deID_sampleID = final de-identified ID used for sharing
deidentified <- deidentified %>%
  select(ID, PT_ID)

# Ensure each Final_deID_sampleID appears only once
# If duplicates exist, keep the first one
deidentified <- deidentified %>%
  distinct(PT_ID, .keep_all = TRUE)

date_of_sampling<-date_of_sampling%>%left_join(deidentified, by="ID")

paired_meta<-date_of_sampling%>%filter(PT_ID%in%paired)

res_df<-res_df%>%left_join(paired_meta[c("PT_ID", "Disease", "BOR_FC")], by="PT_ID")


outdir<-paste0(output_dir, "plots/")
# assicurati che le colonne siano coerenti
res_df <- res_df %>%
  rename(freq_pre = p1, freq_post = p2) %>%
  mutate(
    Disease = ifelse(is.na(Disease), "NA", Disease),
    BOR_FC  = ifelse(is.na(BOR_FC), "NA", BOR_FC)
  )

# --- user settings ---
fdr_thr <- 0.1
large_thr <- 0.01  

# res_df formatting
plot_df <- res_df %>%
  mutate(
    freq_pre  = as.numeric(freq_pre),
    freq_post = as.numeric(freq_post),
    
    # pseudocount for log scale 
    eps = 1e-6,
    freq_pre_pc  = pmax(freq_pre,  eps),
    freq_post_pc = pmax(freq_post, eps),
    
    # label facet: PT_ID + Disease + BOR
    pt_label = paste0(PT_ID, "\n", Disease, " | ", BOR_FC),
    
    # direction 
    change_dir = case_when(
      fdr < fdr_thr & freq_post > freq_pre ~ "expansion",
      fdr < fdr_thr & freq_post < freq_pre ~ "contraction",
      TRUE ~ "stable"
    ),
    
    # small vs large 
    size_class = case_when(
      pmax(freq_pre, freq_post) >= large_thr ~ "large",
      TRUE ~ "small"
    ),
    
    clonal_change = case_when(
      change_dir == "stable" ~ "stable",
      change_dir == "expansion"   & size_class == "small" ~ "small clone expansion",
      change_dir == "expansion"   & size_class == "large" ~ "large clone expansion",
      change_dir == "contraction" & size_class == "small" ~ "small clone contraction",
      change_dir == "contraction" & size_class == "large" ~ "large clone contraction"
    )
  ) %>%
  mutate(
    clonal_change = factor(
      clonal_change,
      levels = c(
        "large clone contraction",
        "small clone contraction",
        "stable",
        "small clone expansion",
        "large clone expansion"
      )
    )
  )

# colori simili alla figura
cols <- c(
  "large clone contraction" = "darkblue",  # magenta
  "small clone contraction" = "#00A6D6",  # cyan
  "stable"                  = "grey60",
  "small clone expansion"   = "#7E57C2",  # purple
  "large clone expansion"   = "#D81B60"   # dark blue
)

# breaks/labels come in figura (percentuali)
brks <- c(1e-4, 1e-3, 1e-2, 1e-1)
labs <- c("0.01%", "0.1%", "1%", "10%")

p_smm <- plot_df%>%filter(Disease=="HRSMM")%>%ggplot( aes(x = freq_pre_pc, y = freq_post_pc)) +
  # crosshair a 1%
  geom_vline(xintercept = large_thr, linewidth = 0.4) +
  geom_hline(yintercept = large_thr, linewidth = 0.4) +
  
  # diagonal reference
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50", linewidth = 0.4) +
  
  # points + jitter
  geom_jitter(
    aes(color = clonal_change),
    width  = 0.08,   
    height = 0.08,
    size = 1.6,
    alpha = 0.5
  ) +
  
  scale_x_log10(breaks = brks, labels = labs ) +
  scale_y_log10(breaks = brks, labels = labs) +
  scale_color_manual(values = cols, name = "Clonal change category:") +
  
  facet_wrap(~ pt_label, ncol = 4) +
  labs(
    x = "Relative abundance pre TEC (log10)",
    y = "Relative abundance post TEC (log10)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

p_smm

ggsave(p_smm, filename=paste0(outdir, "clonotype_tracking_smm_tec.pdf"), width=8, height=6)
p_rrmm <- plot_df%>%filter(Disease=="RRMM")%>%ggplot( aes(x = freq_pre_pc, y = freq_post_pc)) +
  # crosshair a 1%
  geom_vline(xintercept = large_thr, linewidth = 0.4) +
  geom_hline(yintercept = large_thr, linewidth = 0.4) +
  
  # diagonal reference
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50", linewidth = 0.4) +
  
  # points + jitter
  geom_jitter(
    aes(color = clonal_change),
    width  = 0.08,   
    height = 0.08,
    size = 1.6,
    alpha = 0.5
  ) +
  
  scale_x_log10(breaks = brks, labels = labs ) +
  scale_y_log10(breaks = brks, labels = labs) +
  scale_color_manual(values = cols, name = "Clonal change category:") +
  
  facet_wrap(~ pt_label, ncol = 4) +
  labs(
    x = "Relative abundance pre TEC (log10)",
    y = "Relative abundance post TEC (log10)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

p_rrmm

ggsave(p_rrmm, filename=paste0(outdir, "clonotype_tracking_rrmm_tec.pdf"), width=8, height=6)

#save res_df
write_csv(res_df, paste0(output_dir, "Poisson_TEC.csv"))

# -------------------------------------------------------------------------
#endogenous CAR-T
# -------------------------------------------------------------------------

endogenous_t_cells<-qread( "~/Immune_project/data_to_share_CART_TEC/endogenous_t_cells_deid.qs")

remove <- unique(rownames(endogenous_t_cells@meta.data[is.na(endogenous_t_cells@meta.data$clonotypeID), ]))
keep<-setdiff(unique(rownames(endogenous_t_cells@meta.data)), remove)
T_cells_filt_clono<-subset(endogenous_t_cells, cells = keep)
T_cells_filt_clono@meta.data<-T_cells_filt_clono@meta.data%>%mutate(combined_clonotype = paste(TRA_cdr3_aa, TRB_cdr3_aa, sep="_"),
                                                                    new_clones=paste0(PT_ID,"_",combined_clonotype)) 
#take only Diseased, TEC, PB
T_cells_filt_clono<-subset(T_cells_filt_clono, Product=="Cilta-cel"&Tissue=="PB"&Timepoint==c("Pre", "Day30"))
meta<-T_cells_filt_clono@meta.data
meta<-meta%>%select(c("PT_ID", "combined_clonotype", "Timepoint"))
meta<-meta%>%mutate(sample=paste0(PT_ID, "_", Timepoint))
total<-meta%>%group_by(sample)%>%summarise(total=n())
sub<-meta%>%group_by(PT_ID,sample, Timepoint, combined_clonotype)%>%summarise(n=n())

paired<-sub%>%group_by(PT_ID)%>%summarise(n=n_distinct(Timepoint))%>%filter(n>1)%>%pull(PT_ID)
sub<-sub%>%filter(PT_ID%in%paired)

patients<-unique(sub$PT_ID)
res_list <- list()


# Poisson test with sample sizes
poisson_test_with_size <- function(count1, count2, size1, size2) {
  result <- poisson.test(c(count1, count2), T = c(size1, size2))
  return(result$p.value)
}

# Get results over 2 lists
poisson_testing = function(counts1, counts2){
  
  size1=sum(counts1)
  size2=sum(counts2)
  
  p_values <- mapply(poisson_test_with_size, counts1, counts2, 
                     MoreArgs = list(size1 = size1, size2 = size2))
  
  as_tibble(data.frame(
    n1 = counts1,
    n2 = counts2,
    p1 = counts1 / size1,
    p2 = counts2 / size2,
    p = p_values
  )) %>% mutate(fdr = p.adjust(p, method = 'BH'), sig = fdr < 0.05)
}


for(pt in patients){
  try<-sub%>%ungroup()%>%filter(PT_ID==pt)%>%select(-c("sample"))
  sub_wide<-try%>%pivot_wider(names_from = Timepoint, values_from = n, values_fill = 0)
  sub_wide<-sub_wide%>%arrange(combined_clonotype)
  
  pre<-unlist(sub_wide$Pre)
  post<-unlist(sub_wide$Day30)
  
  total<-meta%>%filter(PT_ID==pt)%>%group_by(Timepoint)%>%summarise(n=n())
  tot_pre<-total$n[2]
  tot_post<-total$n[1]
  
  
  res<-poisson_testing( pre, post) %>%
    mutate(combined_clonotype = sub_wide$combined_clonotype)%>%
    arrange(fdr)
  res_list[[pt]] <- res
  
}

res_df <- bind_rows(res_list, .id = "PT_ID")

#extract meta
date_of_sampling <- T_cells_filt_clono@meta.data%>%distinct(PT_ID, .keep_all = TRUE)
res_df<-res_df%>%left_join(date_of_sampling[c("PT_ID", "Disease", "BOR_FC")], by="PT_ID")


outdir<-paste0(output_dir, "plots/")
# assicurati che le colonne siano coerenti
res_df <- res_df %>%
  rename(freq_pre = p1, freq_post = p2) %>%
  mutate(
    Disease = ifelse(is.na(Disease), "NA", Disease),
    BOR_FC  = ifelse(is.na(BOR_FC), "NA", BOR_FC)
  )

# --- user settings ---
fdr_thr <- 0.1
large_thr <- 0.01  

# res_df formatting
plot_df <- res_df %>%
  mutate(
    freq_pre  = as.numeric(freq_pre),
    freq_post = as.numeric(freq_post),
    
    # pseudocount for log scale 
    eps = 1e-6,
    freq_pre_pc  = pmax(freq_pre,  eps),
    freq_post_pc = pmax(freq_post, eps),
    
    # label facet: PT_ID + Disease + BOR
    pt_label = paste0(PT_ID, "\n", Disease, " | ", BOR_FC),
    
    # direction 
    change_dir = case_when(
      fdr < fdr_thr & freq_post > freq_pre ~ "expansion",
      fdr < fdr_thr & freq_post < freq_pre ~ "contraction",
      TRUE ~ "stable"
    ),
    
    # small vs large 
    size_class = case_when(
      pmax(freq_pre, freq_post) >= large_thr ~ "large",
      TRUE ~ "small"
    ),
    
    clonal_change = case_when(
      change_dir == "stable" ~ "stable",
      change_dir == "expansion"   & size_class == "small" ~ "small clone expansion",
      change_dir == "expansion"   & size_class == "large" ~ "large clone expansion",
      change_dir == "contraction" & size_class == "small" ~ "small clone contraction",
      change_dir == "contraction" & size_class == "large" ~ "large clone contraction"
    )
  ) %>%
  mutate(
    clonal_change = factor(
      clonal_change,
      levels = c(
        "large clone contraction",
        "small clone contraction",
        "stable",
        "small clone expansion",
        "large clone expansion"
      )
    )
  )

# colori simili alla figura
cols <- c(
  "large clone contraction" = "darkblue",  # magenta
  "small clone contraction" = "#00A6D6",  # cyan
  "stable"                  = "grey60",
  "small clone expansion"   = "#7E57C2",  # purple
  "large clone expansion"   = "#D81B60"   # dark blue
)

# breaks/labels come in figura (percentuali)
brks <- c(1e-4, 1e-3, 1e-2, 1e-1)
labs <- c("0.01%", "0.1%", "1%", "10%")

p_smm <- plot_df%>%filter(Disease=="HRSMM")%>%ggplot( aes(x = freq_pre_pc, y = freq_post_pc)) +
  # crosshair a 1%
  geom_vline(xintercept = large_thr, linewidth = 0.4) +
  geom_hline(yintercept = large_thr, linewidth = 0.4) +
  
  # diagonal reference
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50", linewidth = 0.4) +
  
  # points + jitter
  geom_jitter(
    aes(color = clonal_change),
    width  = 0.08,   
    height = 0.08,
    size = 1.6,
    alpha = 0.5
  ) +
  
  scale_x_log10(breaks = brks, labels = labs ) +
  scale_y_log10(breaks = brks, labels = labs) +
  scale_color_manual(values = cols, name = "Clonal change category:") +
  
  facet_wrap(~ pt_label, ncol = 4) +
  labs(
    x = "Relative abundance pre Cilta (log10)",
    y = "Relative abundance post Cilta (log10)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

p_smm

ggsave(p_smm, filename=paste0(outdir, "clonotype_tracking_smm_CAR.pdf"), width=8, height=6)
p_rrmm <- plot_df%>%filter(Disease=="RRMM")%>%ggplot( aes(x = freq_pre_pc, y = freq_post_pc)) +
  # crosshair a 1%
  geom_vline(xintercept = large_thr, linewidth = 0.4) +
  geom_hline(yintercept = large_thr, linewidth = 0.4) +
  
  # diagonal reference
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50", linewidth = 0.4) +
  
  # points + jitter
  geom_jitter(
    aes(color = clonal_change),
    width  = 0.08,   
    height = 0.08,
    size = 1.6,
    alpha = 0.5
  ) +
  
  scale_x_log10(breaks = brks, labels = labs ) +
  scale_y_log10(breaks = brks, labels = labs) +
  scale_color_manual(values = cols, name = "Clonal change category:") +
  
  facet_wrap(~ pt_label, ncol = 4) +
  labs(
    x = "Relative abundance pre Cilta (log10)",
    y = "Relative abundance post Cilta (log10)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

p_rrmm

ggsave(p_rrmm, filename=paste0(outdir, "clonotype_tracking_rrmm_car.pdf"), width=8, height=6)

#save res_df
write_csv(res_df, paste0(output_dir, "Poisson_CART.csv"))
res_df<-read_csv( paste0(output_dir, "Poisson_CART.csv"))
res_df2<-read_csv( paste0(output_dir, "Poisson_TEC.csv"))
res_df<-rbind(res_df, res_df2)
res_df <- res_df %>%
  mutate(
    category = case_when(
      sig == TRUE  & freq_post > freq_pre ~ "expanding",
      sig == TRUE  & freq_post < freq_pre ~ "contracting",
      sig == FALSE                        ~ "stable",
      TRUE                                ~ NA_character_
    )
  )


# -------------------------------------------------------------------------
#now to understand how many clones are expanding post-treatment
#left join by PT_ID and clonotype so that each clone pre and post is assigned a category
#then select only post-treatment samples and perform expanding clones/ total clones
# -------------------------------------------------------------------------
T_cells_filt2<-qread("~/Immune_project/data_to_share_CART_TEC/T_cells_nov25_deid.qs")

remove <- unique(rownames(T_cells_filt2@meta.data[is.na(T_cells_filt2@meta.data$clonotypeID), ]))
keep<-setdiff(unique(rownames(T_cells_filt2@meta.data)), remove)
T_cells_filt_clono<-subset(T_cells_filt2, cells = keep)
T_cells_filt_clono@meta.data<-T_cells_filt_clono@meta.data%>%mutate(combined_clonotype = paste(TRA_cdr3_aa, TRB_cdr3_aa, sep="_"),
                                                                    new_clones=paste0(PT_ID,"_",combined_clonotype)) 
#take only Diseased, TEC, PB
T_cells_filt_clono<-subset(T_cells_filt_clono, Therapy=="TEC"&Tissue=="PB")
paired_ids<-T_cells_filt_clono@meta.data%>%group_by(PT_ID)%>%summarise(n=n_distinct(Timepoint))%>%filter(n>1)%>%pull(PT_ID)
T_cells_filt_clono<-subset(T_cells_filt_clono, PT_ID%in%paired_ids)

#start with TEC
res_df<-res_df%>%select(-c("Disease"))
#assign the label "category" to each clonotype
T_cells_filt_clono@meta.data<-T_cells_filt_clono@meta.data%>%left_join(res_df, by=c("PT_ID", "combined_clonotype"))
#select Post only
post_tec<-T_cells_filt_clono@meta.data%>%filter(Timepoint=="Post")

#percentage of patients with at least one expanding clones
pt_exp <- post_tec %>%
  group_by(PT_ID, Disease) %>%
  summarise(
    has_expanding = any(category == "expanding" & n() > 0),
    .groups = "drop"
  )

# tbl 2x2
tab <- pt_exp %>%
  count(Disease, has_expanding) %>%
  pivot_wider(names_from = has_expanding, values_from = n, values_fill = 0) %>%
  arrange(Disease)

# matrice per fisher
mat <- as.matrix(tab[, c("FALSE", "TRUE")])
rownames(mat) <- tab$Disease
colnames(mat) <- c("no_expanding", "has_expanding")

mat
fisher_res <- fisher.test(mat)
fisher_res

df_mat <- as.data.frame(mat) %>%
  rownames_to_column("Disease") %>%
  pivot_longer(
    cols = c(no_expanding, has_expanding),
    names_to = "status",
    values_to = "n")

df_mat_perc <- df_mat %>%
  group_by(Disease) %>%
  mutate(
    perc = n / sum(n) * 100
  ) %>%
  ungroup()

df_mat_perc$status <- factor(
  df_mat_perc$status,
  levels = c("no_expanding", "has_expanding"),
  labels = c("Not Expanding", "Expanding")
)

p<-ggplot(df_mat_perc, aes(x = Disease, y = perc, fill = status)) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Not Expanding"    = "lightgrey",
      "Expanding" = "steelblue"   # oppure "firebrick"
    )
  ) +
  labs(
    y = "% of Patients",
    x = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 22) +
  theme(
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

ggsave(p, filename="~/TCR_project/plots/n_pts_expanding_TEC.pdf", width=6, height=6)

#number of CD4 and CD8 cells expanding in each post sample
sub <- post_tec %>%
  group_by(PT_ID, Disease,lineage, category) %>%
  summarise(n = n())

#fraction of expanding 
sub<-sub%>%group_by(PT_ID, lineage)%>%mutate(tot=sum(n), fraction=n/tot)

sub_expanding_cd8 <- sub %>%
  filter(category == "expanding"&lineage=="CD8") %>%
  mutate(Disease = factor(Disease, levels = c("HRSMM", "RRMM"))) %>%
  mutate(jitter_x = as.numeric(Disease) + runif(n(), -0.12, 0.12))

boxplot_plot_final <- ggplot(sub_expanding_cd8, aes(x = Disease, y = fraction * 100, fill = Disease)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.7  ) +
  geom_point(aes(x = jitter_x), fill = "darkgrey", alpha = 0.5,
             shape = 21, size = 4, stroke = 0.6, color = "black") +
  geom_pwc(method="wilcox.test", method.args = list(alternative="greater"), label.size=5) +
  scale_fill_manual(values = c("HRSMM" = "lightgrey", "RRMM" = "gold")) +
  labs(y = "% of CD8+ T cells", x = NULL, title="Cells in an expanding clone") +
  theme_classic(base_size = 18) +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 13),
        legend.position = "none")

ggsave(boxplot_plot_final, filename = paste0("~/TCR_project/plots/CD8_expanding.pdf"), width = 6, height = 4)

write_csv(sub, "~/Immune_project/data_to_share_CART_TEC/clonotype_tracking_TEC.csv")


# -------------------------------------------------------------------------
#select T cells in expanding clones 
# -------------------------------------------------------------------------

post_tec_exp <- post_tec %>%
  filter(category == "expanding") %>%
  count(PT_ID, final_annotation, name = "n")   # equivalente a group_by+summarise(n=n())

# lista completa delle annotazioni da includere (coerente su tutti i PT)
all_ann <- sort(unique(post_tec_exp$final_annotation))  # oppure unique(post_tec_exp$final_annotation)

post_tec_exp <- post_tec_exp %>%
  complete(
    PT_ID,
    final_annotation = all_ann,
    fill = list(n = 0)
  )
pt_disease <- post_tec %>%
  distinct(PT_ID, Disease)

post_tec_exp <- post_tec_exp %>%
  left_join(pt_disease, by = "PT_ID")

post_tec_exp <- post_tec_exp %>%
  group_by(PT_ID) %>%
  mutate(
    total = sum(n),
    fraction = ifelse(total > 0, n / total, 0)
  ) %>%
  ungroup()

mat <- post_tec_exp %>%
  select(PT_ID, final_annotation, fraction) %>%
  pivot_wider(names_from = final_annotation, values_from = fraction) %>%
  column_to_rownames("PT_ID") %>%
  as.matrix()

# 2) annotation_row con Disease
ann_row <- post_tec_exp %>%
  distinct(PT_ID, Disease) %>%
  column_to_rownames("PT_ID")

# allinea all'ordine delle righe di mat
ann_row <- ann_row[rownames(mat), , drop = FALSE]

# 3) ordina righe per Disease (es. HRSMM sopra, RRMM sotto)
ord <- order(ann_row$Disease)
mat <- mat[ord, , drop = FALSE]
ann_row <- ann_row[ord, , drop = FALSE]

# 4) colori Disease
ann_colors <- list(
  Disease = c(
    HRSMM = "grey70",
    RRMM  = "gold"
  )
)

# 5) heatmap
pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  border_color = NA,
  fontsize_row = 8,
  fontsize_col = 10
)


# -------------------------------------------------------------------------
#For each patient calculate the percentage of expanding cells in each phenotype
# -------------------------------------------------------------------------
ggplot(post_tec_exp, aes(x=final_annotation, y=fraction))+geom_point()+geom_boxplot()
