#! /usr/bin/env Rscript
if(!require(dplyr)){install.packages("dplyr")}
if(!require(purrr)){install.packages("purrr")}
if(!require(data.table)){install.packages("data.table")}
if(!require(tools)){install.packages("tools")}
if(!require(ggplot2)){install.packages("ggplot2")}



############
#PARAMETERS
############
args <- commandArgs(trailingOnly=T)
markers_per_tissue <- 500
min_markers <- 1 #number of markers that a tissue must have to be included
path_markers <- args[1]
id_atlas <- args[2]
############



############
#MAKE ATLAS
############
atlas <- data.frame()
excluded_types <- c()
#markers calculated with wgbstools find markers function
for (i in list.files(path = path_markers, pattern="Markers",full.names = T)){
  
  markers <- read.delim(i,header=T)
  if (nrow(markers) < 500){
    type <- strsplit(basename(i),split="\\.")[[1]][2]
    print(paste("Warning: ",type, "only contains",nrow(markers), "markers."))
  } else{}
  
  if (nrow(markers) >= min_markers){ #discarding tissues with not enough markers (here colon-fib)
    markers <- head(markers,n=markers_per_tissue)
    atlas <- bind_rows(atlas,markers)
  } else{
    type <- strsplit(basename(i),split="\\.")[[1]][2]
    print(paste(type, "was not included because it only contains",nrow(markers), "markers."))
    excluded_types <- append(excluded_types,type)
  }
}

atlas_nodups <- atlas %>%
  group_by(region) %>%
  filter(n() == 1 ) %>%
  ungroup() %>%
  group_by(startCpG) %>%
  arrange(.by_group = T)


write.table(atlas_nodups,paste0("markers/all_markers_",id_atlas,".tsv"),quote = F,row.names = F,col.names = T,sep = "\t")