#install packages
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz",
  repos = NULL,
  type = "source"
)

remotes::install_github("HetzDra/turboGliph")

#import dependencies
library(qs)
library(turboGliph)
library(dplyr)
T_cells_filt2<-qread("~/Immune_project/data_to_share_CART_TEC/T_cells_nov25_deid.qs")
endogenous_t_cells<-qread( "~/Immune_project/data_to_share_CART_TEC/endogenous_t_cells_deid.qs")

df<-T_cells_filt2@meta.data%>%filter(!is.na(clonotypeID))
keep<-unique(rownames(df))
T_cells_filt_clono<-subset(T_cells_filt2, cells=keep, final_annotation%in%c("CD8+GZMB+TEM",
                                                                            "CD8+TNFRSF9+",
                                                                            "CD8+GZMK+TEM",
                                                                            "CD8+IFN+",
                                                                            "CD8+TCM",
                                                                            "CD8+KIR+TEM"
                                                                            ))

counts<-as.data.frame(T_cells_filt_clono@meta.data%>%group_by(Final_deID_sampleID)%>%summarise(counts=n_distinct(clonotypeID)))
T_cells_filt_clono@meta.data<-T_cells_filt_clono@meta.data%>%left_join(counts, by="Final_deID_sampleID")
sub<-T_cells_filt_clono@meta.data%>%select(c("Final_deID_sampleID", "TRB_cdr3_aa", "TRB_V_gene", "counts"))
sub<-sub%>%filter(!is.na(TRB_cdr3_aa))

colnames(sub)<-c("patient","CDR3b","TRBV", "counts")
sub <- sub %>% 
  dplyr::select(CDR3b, dplyr::everything())


# -------------------------------------------------------------------------
#do the same for CAR-T cells
# -------------------------------------------------------------------------
df<-endogenous_t_cells@meta.data%>%filter(!is.na(clonotypeID))
keep<-unique(rownames(df))
clono<-subset(endogenous_t_cells, cells=keep, final_annotation%in%c("CD8+GZMB+TEM",
                                                                            "CD8+GATA3+",
                                                                            "CD8+GZMK+TEM",
                                                                            "CD8+IFN+",
                                                                            "CD8+TCM",
                                                                            "CD8+KIR+TEM",
                                                                            "CD8+Cycling"
))

counts2<-as.data.frame(clono@meta.data%>%group_by(Final_deID_sampleID)%>%summarise(counts=n_distinct(clonotypeID)))
clono@meta.data<-clono@meta.data%>%left_join(counts2, by="Final_deID_sampleID")
sub2<-clono@meta.data%>%select(c("Final_deID_sampleID", "TRB_cdr3_aa", "TRB_V_gene", "counts"))
sub2<-sub2%>%filter(!is.na(TRB_cdr3_aa))

colnames(sub2)<-c("patient","CDR3b","TRBV", "counts")
sub2 <- sub2 %>% 
  dplyr::select(CDR3b, dplyr::everything())

final<-rbind(sub,sub2)
final_dedup <- final %>%
  distinct(patient, CDR3b, .keep_all = TRUE)

# -------------------------------------------------------------------------
# To use **turbo_gliph**, **gliph2** and **gliph_combined**, 
# CDR3 beta sequences are required. These must be specified either by a character vector 
# or in a dataframe in a column named *CDR3b*. 
# Additional information is not required for clustering, but is recommended for automatic cluster scoring. 
# This includes the following information:
#   
# * column *TRBV* : V gene information of the particular sequence
# * column *patient* : index number or similar of the patient from which the clone was isolated
# * column *counts* : number of clones in the sample
# * column *HLA* : HLA alleles of the particular patient, separated by commas
# 
# Notice: For **turbo_gliph**, **gliph2** and **gliph_combined** the column names are important, the order of the columns can be arbitrary. In this example, four additional columns (TRBJ, CDR3a, TRAV, TRAJ) are given but will not be used by the GLIPH algorithms. However, this additional information is stored by the functions and will be visible in the returned clusters.The additional columns do not affect the algorithms as long as all column names specified in the input and output in this tutorial are avoided.<br><br>
#   
#   
#   # Individual reference repertoire
#   By default, the GLIPH algorithms in this package use the naive reference repertoire from the original GLIPH publication. However, there are several other reference repertoires available on the GLIPH2 website\(^4\) for different T cell subtypes as well as different species, which can be used in this package as well as user-created reference repertoires. For this, the reference sequences must be passed to the functions as a data frame under the parameter *refdb_beta*. In any case, a column named *CDR3b* is required which contains the sequences. Optionally (but mandatory for some parameter settings), a column named *TRBV* is expected to contain the V genes of the sequences. Additional information is ignored by the functions.
# 
# # Function: turbo_gliph
# 
# The turbo_gliph function is a runtime-efficient implementation of the GLIPH algorithm provided by Glanville&nbsp;*et&nbsp;al.* in their publication as a Perl script.\(^1\)\(^,\)\(^2\) The choice of the R programming language, a restructuring of some parts of the original code, and the incorporated option of parallelization allow a significant reduction in runtime, especially with respect to larger sample input sizes. The ability to customize the algorithm using a variety of input parameters was retained and extended to include some parameters that are listed on the github website but not included in the Perl script. Details about the algorithm and the input parameters can be found in the original publication, the github website and the documentation of this package.
# 
# ## Brief overview of the algorithm
# 1. Preparation of the input and the reference repertoire
# + Exclusion of CDR3b sequences with characters outside the amino acid one-letter code.
# + Exclusion of CDR3b sequences with length below a customizable threshold.
# + If desired, exclusion of CDR3b sequences that do not start with the amino acid cysteine and end with the amino acid phenylalanine.
# 2. Identification of all desired motifs as local similarity (default: 2, 3 and 4 kmers)
# + Identification and summary of all motifs in the sample sequences
# + Repeated random sampling from the reference repertoire with sample set depth and identification and summary of all motifs in these subsamples.
# 3. Identification of significantly enriched motifs in the sample set by the following criteria:
#   + Falling below an adjustable p-value (= probability of encountering a motif in repeated random sampling at least as often as in the sample set).
# + Exceeding an adjustable fold change (= factor of the overrepresentation of the motif in the sample set compared to the mean value of the repeated random sampling)
# + Exceeding a minimum frequency to reduce contamination of significant motifs with rare motifs.
# 4. Connecting all sequences with a Hamming distance below a customizable threshold as an expression of global similarity
# 5. Finding all connections of the sample sequences via identified local and global similarities and clustering based on these connections
# 6. Scoring the clusters based on the following quantities:
#   + cluster size
# + CDR3b length enrichment
# + Clonal expansion enrichment (optional)
# + V gene enrichment (optional)
# + Enrichment of common HLA alleles (optional).
# + A summary score for all other scores (similar to the product of all scores and an additional factor)
# 
# 
# ## Calling the function and output explanation
# The sequences are simply analyzed by calling the function **turbo_gliph** with the input dataframe as *cdr3_sequences* parameter.


sub
# -------------------------------------------------------------------------
res_turbogliph <- turboGliph::turbo_gliph(cdr3_sequences = final_dedup,
                                          n_cores = 8)
saveRDS(res_turbogliph, "~/TCR_project/res_turbogliph.rds")
plot_network(clustering_output = res_turbogliph,
             n_cores = 1)

