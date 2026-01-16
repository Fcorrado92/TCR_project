library(readr)
vdjdb<-read_tsv("~/TCR_project/SearchTable-2026-01-16 00_26_28.282.tsv")


norm_hla <- function(x) {
  x %>%
    str_to_upper() %>%
    str_replace_all("\\s+", "") %>%
    str_replace("^HLA-", "") %>%        # remove HLA- if present
    str_replace_all("HLA", "") %>%      # sometimes "HLAA*02:01" weirdness
    str_replace_all("^-", "")           # safety
}
hla_long <- hla_tbl_patient %>%
  mutate(HLA = str_split(HLA, ",")) %>%
  unnest(HLA) %>%
  mutate(HLA = norm_hla(HLA)) %>%
  distinct(patient_key, HLA)


hla_long <- hla_long %>%
  mutate(
    HLA_4digit  = str_replace(HLA, "^(.*\\*[0-9]{2}:[0-9]{2}).*$", "\\1")  # -> C*07:02
  )


vdjdb2 <- vdjdb %>%
  mutate(mhc = norm_hla(vdjdb$`MHC A`),
             V = str_remove(V, "\\*.*$")
           )

colnames(final3)<-c("CDR3", "V")

vdjdb2<-vdjdb2%>%filter(Gene=="TRB")

tcr2 <- final3 %>%
  mutate(
    CDR3  = str_to_upper(str_replace_all(CDR3b, "\\s+", "")),
    V = str_to_upper(TRBV)  )
tcr2 <- tcr2 %>%
  mutate(patient_key = str_extract(patient, "P\\d+"))



hits <- tcr2 %>%
  left_join(
    vdjdb2,
    by = c("CDR3" = "CDR3", "V" = "V")
  )

hits<-hits%>%select(-c("HLA"))
