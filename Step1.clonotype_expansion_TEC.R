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


res<-poisson_testing( pre, post) %>% arrange(fdr)%>%
  mutate(combined_clonotype = sub_wide$combined_clonotype)
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
  
  
  res<-poisson_testing( pre, post) %>% arrange(fdr)%>%
    mutate(combined_clonotype = sub_wide$combined_clonotype)
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

