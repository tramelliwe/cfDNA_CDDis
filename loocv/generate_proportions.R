#! /usr/bin/env Rscript

suppressPackageStartupMessages(if(!require(dplyr)){install.packages("dplyr")})

setwd("/home/cyril")

args <- commandArgs(trailingOnly=T)
output_file <- args[1]

files_groups <- read.delim("/home/cyril/atlas_uxm_02/groups.csv",sep=",")
groups <- unique(files_groups$group)
groups_enough <- c()

for (i in groups){
  files <- filter(files_groups,group==i)
  if (nrow(files)>1){ #only include groups that have more than one sample
    groups_enough <- append(groups_enough,i)
  }
}



total_bloodcells <- sample(seq(200,800),1) #in real samples, cfDNA deriving from 'blood cells' is the majority. Here, we reflect this by imposing a total proportion between 20 and 80%.
therest <- 1000-total_bloodcells # other tissues make up the rest.

x <- sample(seq(0,500),6)
x <- x/(sum(x)/total_bloodcells)
bloodcellgroup <- c("Erythrocyte_progenitors",
                    "Blood_T_cells",
                    "Blood_NK_cells",
                    "Monocytes_Macrophages",
                    "Blood_Granulocytes",
                    "Blood_B_cells")
df_props <- data.frame(tissue=bloodcellgroup,
                       proportion=x)

others <- groups_enough[-c(19, 20,21,22,23,24)] #remove the blood cells


#pick random tissue
#uniform distribution between 0 and therest

#give more chance to tissue of interest to be chosen first
tissue_interest <- c("Liver_Hepatocytes", "Heart_Cardiomyocytes","Heart_Fibroblasts","Neurons",
                     "Kidney_Epithelium","Lung_Alveolar_Epithelium","Lung_Bronchial_Epithelium",
                     "Oligodendrocytes")
others <- others[-which(others %in% tissue_interest)]

while (therest >= 1){ #while still 0.1% available 
  if (length(tissue_interest) > 0) {
    i <- sample(seq(1, length(tissue_interest)), 1)
    tissue <- tissue_interest[i]
    tissue_interest <- tissue_interest[-i]
  } else if (length(others) > 0) {
    i <- sample(seq(1, length(others)), 1)
    tissue <- others[i]
    others <- others[-i]
  } else {
    break
  }
  possible_prop <- seq(0, therest)
  prop <- sample(possible_prop, 1, prob = (possible_prop^2) / sum(possible_prop^2))  # higher probability for higher numbers
  therest <- therest - prop
  s <- data.frame(tissue = tissue, proportion = prop)
  df_props <- bind_rows(df_props, s)
}

#add all tissues not only others
others <- data.frame(tissue=c(others,tissue_interest),proportion=0)
df_props <- bind_rows(df_props,others)
df_props <- df_props[order(df_props$tissue),]
df_props$proportion <- df_props$proportion/1000
write.table(df_props, file = output_file, quote = F,sep = ",",row.names = F,col.names = T)
print(paste("Successfully written proportions v2.0",output_file))


