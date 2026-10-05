#sdm data sheet prep

setwd("~/Desktop/CNR/SDM Bacteria/files for analysis/")
# #Raw Data Prep 
# #################
# #files are: genus.csv, genome_taxa.csv, metabolism_info.csv, KEGG_modules.csv, sample_metadata.csv, phys_chem.csv
# #Prep of genus dataframes 
abu<-read.csv2("genus.csv") #abundance of genera
library(dplyr)

abu <- abu %>%
  mutate(genus = sub(".*g__", "", gtdb_taxonomy)) %>%
  relocate(genus, .before = 1)

genome<-read.csv2("genome_taxa.csv")

genome <- genome %>%
  mutate(genus = sub(".*g__", "", gtdb_taxonomy)) %>%
  relocate(genus, .before = 1)


excluded_because_new_taxonomy <- anti_join(genome, abu, by = "genus") #I cannot find correspondance in the two datasets
abu <- abu %>%
    filter(genus%in%genome$genus)
genome <- genome %>%
    filter(genus%in%abu$genus)

metabolism<-read.csv2("metabolism_info.csv")
kegg<-read.csv2("KEGG_modules.csv") #info regarding metabolism modules in "metabolism"

summary_by_cat <- metabolism %>%
  group_by(kegg$category3) %>%
  dplyr::summarise(across(where(is.numeric), sum))
summary_by_cat<-summary_by_cat[-1,]
pw_names<-as.data.frame(cbind(summary_by_cat, paste("PW", seq_len(nrow(summary_by_cat)))))
colnames(pw_names)<-c("category3", "PW")
summary_by_cat<-summary_by_cat[,-1]
sum_metabo<-as.data.frame(t(summary_by_cat))
colnames(sum_metabo)<-paste(pw_names$PW)

sum_metabo <- sum_metabo %>%
  mutate(accession = rownames(.)) %>%   # use . instead of dataframe name
  relocate(accession, .before = 1)

all_traits<-inner_join(genome, sum_metabo, by="accession")


sample_metadata<-read.csv2("sample_metadata.csv") #metadata
phys_chem<-read.csv2("phys_chem_data.csv") #metadata
phys_chem$sample<-phys_chem$abbreviation
sample_metadata$sample%in%phys_chem$sample

sample_metadata <- sample_metadata %>%
  mutate(
    sample = gsub("[- ]", ".", sample),                # replace - and space with .
    sample = ifelse(grepl("^[0-9]", sample),          # if starts with a number
                    paste0("X", sample),              # add X
                    sample)
  )
phys_chem <- phys_chem %>%
  mutate(sample = gsub("[- ]", ".", sample),                # replace - and space with .
         sample = ifelse(grepl("^[0-9]", sample),          # if starts with a number
                         paste0("X", sample),              # add X
                         sample)
  )

sum(sample_metadata$sample%in%phys_chem$sample)
no_metabu<-anti_join(sample_metadata, phys_chem, by="sample")


sum(colnames(abu) %in% sample_metadata$sample)
abu2<-as.data.frame(t(abu))
colnames(abu2)<-abu$genus
abu2$sample<-rownames(abu2)

# ##no_metabu<-anti_join(sample_metadata, abu2, by="sample")
sample_metadata<-sample_metadata %>%
  filter(sample%in%abu2$sample)
abu2<-abu2 %>%
  filter(sample%in%sample_metadata$sample) %>%
  relocate(sample, .before = 1)

genus_site<-t(abu2)
genus_site<-as.data.frame(genus_site[-1,] ) #dataframe GENUS per SITE
genus_traits<-all_traits[,-c(2,14,15,16,17,18,19,20)] #dataframe Genus per Traits
taxonomy<-all_traits[,c(1,20)] #Taxonomy and accession
sample_genus<-abu2
sample_metadata<-sample_metadata
#
#
rm(abu, abu2, all_traits, excluded_because_new_taxonomy, genome,  metabolism, sum_metabo, summary_by_cat)
#



pw_all<-left_join(pw_names, kegg, by="category3")
#save.image(file="data_SDM.RData")

#####################

load(file="data_SDM_2.RData")
#prepared dataframes: 

head(genus_site[,1:5]) #genus per site abundances site è sample
head(genus_traits) #genus per traits PW is pw can be indentified in "pw_names"
head(phys_chem) #phys chemical factors, but there are about 400 samples missing which can be inspected in "no_metabu"
head(sample_genus[,1:6]) # transposed genus_site
head(sample_metadata)
head(pw_all) #pathways and their category in case you want to group
head(taxonomy) #taxonomy of genera


