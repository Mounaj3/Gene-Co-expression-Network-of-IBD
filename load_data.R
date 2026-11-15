
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
library(data.table)
library(GEOquery)
library(dplyr)



#---------------------------Metadata------------------------------------#

eset = getGEO( filename ="data/GSE193677_series_matrix.txt",getGPL = FALSE)

meta= pData(eset)      


#Remove "ch1$" 
names(meta) = sub(":ch1$", "", names(meta))
# Extract important info from title column

meta$samples_name=sub(",.*", "", meta$title) # replace what comes after 1st ',' by ''
meta$reGRID_patient= sub("MSCCR_reGRID_(\\d+)_Biopsy_\\d+", "\\1",
                         meta$samples_name) # in : MSCCR_reGRID_1_Biopsy_3, 1 is patient's ID 
colnames(meta)
keep <- c(
  
  "geo_accession", "samples_name", "reGRID_patient","typere",  "demographics_gender" ,       
  "diseasebi","regionre", "ibd_disease",
  "ibd_endoseverity_4levels", "ghas_sum7" ,       # sévérité endoscopique (patient)
  "crp_jjmgl_log2")               # CRP (inflammation systémique), ~4 % de NA


meta=meta[,keep]
dim(meta)        
colnames(meta)
str(meta)
meta[] <- lapply(meta, function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "n/a", "NaN", "Not Available", "not available")] <- NA
  x
})

meta = meta %>% 
  mutate_at(c("crp_jjmgl_log2","ghas_sum7"), as.numeric)

# Exclude the controle 
meta=meta[meta$diseasebi =="IBD",]

# Resampling "NonI" group 
table(meta$typere)
# Statut inflammatoire : "I" ou "NonI"

pat= meta
class(pat)
str(pat)
#Create 3 regions groups
pat= pat %>%
  mutate(region_groups = case_when(
    regionre == "Ileum"  ~ "Ileum",
    regionre == "Rectum" ~ "Rectum",
    regionre %in% c("Cecum", "RightColon", "Transverse", "LeftColon", "Sigmoid") ~ "Colon"
  ))

table(pat$region_groups)
goal <- pat %>%
  dplyr::filter(typere == "I") %>%
  dplyr::count(ibd_disease, region_groups, name = "n_cible")
#Count the samples per group in the Inflamed biopsy samples
goal = pat %>% dplyr::filter(typere == "I") %>% dplyr::count(ibd_disease, region_groups, name = "n_cible")
goal

set.seed(123)
NonI_samples = pat %>%
  filter(typere == "NonI") %>%
  inner_join(goal, by = c("ibd_disease", "region_groups")) %>%
  group_by(ibd_disease, region_groups) %>%
  slice_sample(prop = 1) %>%
  filter(row_number() <= n_cible) %>%
  ungroup()
NonI_samples=NonI_samples[,!names(NonI_samples) %in% "n_cible"]
dplyr::count(NonI_samples, ibd_disease, region_groups) 

I_samples= pat %>% dplyr::filter(typere== "I")

final_samples= bind_rows( I_samples,NonI_samples)

table(final_samples$typere)   
table(final_samples$typere, final_samples$region_groups)
str(final_samples)



#-----------------------Gene expression ------------------------------------#

#Link to seq_RNA data: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE193677
# We have two count files: row counts and filtred 
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


fwrite(as.data.table(counts, keep.rownames = "samples_name"),
       "data/counts_selected_samples.csv")
fwrite(final_samples, "data/meta_selected_samples.csv")


