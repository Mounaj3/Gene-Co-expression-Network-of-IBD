
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
library(data.table)
library(GEOquery)
library(dplyr)



#---------------------------Metadata------------------------------------#

#meta <- read.csv("data/SraRunTable.csv")   # adaptez le chemin


eset = getGEO(filename = "data/GSE193677_series_matrix.txt", getGPL = FALSE)
meta= pData(eset)

#meta$title
names(meta)

#Remove "ch1$" 
names(meta) = sub(":ch1$", "", names(meta))
# Extract important info from title column

meta$samples_name=sub(",.*", "", meta$title) # replace what comes after 1st ',' by ''
meta$reGRID_patient= sub("MSCCR_reGRID_(\\d+)_Biopsy_\\d+", "\\1",
                           meta$samples_name) # in : MSCCR_reGRID_1_Biopsy_3, 1 is patient's ID 


keep= c("geo_accession", "samples_name", "reGRID_patient",
          "typere", "regionre", "ibd_disease", "diseasetypere", "diseasebi",
          "demographics_gender", "study_eligibility_age_at_endo",
          "crp_jjmgl_log2", "log2_fecalcalpro_mgperg",
          "harveybradshawindex_hbi_score", "colitisactivityindex_sccai",
          "ibd_clinicianmeasure_inactive_active",
          "ibdsescd_totalsescd", "ibdmesuc_mayo_score", "ibd_endoseverity_4levels", "endoremiss",
          "ghas_sum7", "nancyindex", "max_ghas_sum7", "max_nancy", "historemiss")


meta=meta[,keep]


dim(meta)        
colnames(meta)

# Exclude the controle 

meta=meta[meta$diseasebi =="IBD",]

# Resampling "NonI" group 
table(meta$typere)


# Statut inflammatoire : "I" ou "NonI"

pat= meta
#Create 3 regions groups
pat= pat %>%
  mutate(region_groups = case_when(
    regionre == "Ileum"  ~ "Ileum",
    regionre == "Rectum" ~ "Rectum",
    regionre %in% c("Cecum", "RightColon", "Transverse", "LeftColon", "Sigmoid") ~ "Colon"
  ))

table(pat$region_groups)

#Count the effective per group in the Inflamed biopsy samples
goal = pat %>% filter(typere == "I") %>% count(ibd_disease, region_groups, name = "n_cible")
goal

set.seed(123)
NonI_samples = pat %>%
  filter(typere == "NonI") %>%
  inner_join(goal, by = c("ibd_disease", "region_groups")) %>%
  group_by(ibd_disease, region_groups) %>%
  slice_sample(prop = 1) %>%
  filter(row_number() <= n_cible) %>%
  ungroup()

count(NonI_samples, ibd_disease, region_groups) 



I_samples= pat %>% filter(typere== "I")
final_samples= bind_rows(
  I_samples,
  NonI_samples
)

table(final_samples$typere)   
table(final_samples$typere, final_samples$region_groups)



#-----------------------Gene expression ------------------------------------#

#Link to seq_RNA data: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE193677
# We have two count files: all genes and filtred 
load_expr= function(path){
  
  counts=fread(path)
  genes_id=counts$V1
  samples=colnames(counts[,-1])  
  mat=t(as.matrix(counts[,-1]))
  rownames(mat)= samples
  colnames(mat)=genes_id
  return(list(
    expr=mat,
    samples=samples,
    genes_id=genes_id
  ))
}


#------------------------Prepare datasets --------------------------------#

loaded_counts=load_expr("data/GSE193677_MSCCR_Biopsy_counts.txt")
loaded_counts$expr[1:3,1:3]
dim(loaded_counts$expr)
length(loaded_counts$genes_id)
length(loaded_counts$samples)
counts=loaded_counts$expr

# Keep counts of selected samples 
counts= counts[final_samples$samples_name,]
dim(counts)
dim(final_samples)


# From the data.frame / data.table with the "gene" column

fwrite(as.data.table(counts, keep.rownames = "samples_name"),
       "data/counts_selected_samples.csv")
# Compressed (much smaller file, R reads it directly)
fwrite(final_samples, "data/meta_selected_samples.csv")


