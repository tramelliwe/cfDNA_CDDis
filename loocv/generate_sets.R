#! /usr/bin/env Rscript
if(!require(xlsx)){install.packages("xlsx")}
if(!require(stringr)){install.packages("stringr")}
if(!require(dplyr)){install.packages("dplyr")}

setwd("/home/cyril")

args <- commandArgs(trailingOnly=T)

files_groups <- read.delim("/home/cyril/workspace/atlas_uxm_02/groups.csv",sep=",") #links sample ids to one of the 39 groups
groups <- unique(files_groups$group)

atlas_start <- as.numeric(args[1]) #starting fold id
nr_atlas <- as.numeric(args[2]) #how many folds
output_folder <- args[3]

for (k in seq(atlas_start,nr_atlas+nr_atlas)){
  print(paste("Generating set",k, "on", atlas_start+nr_atlas))
  training_set <- data.frame(file_basename=c(),group=c())
  test_set <- data.frame(file_basename=c(),group=c())
  for (i in groups){
    files <- filter(files_groups,group==i)
    if (nrow(files)>1){ #if a group only has a single sample, then cross-validation of this group is not possible.
      test_file_id <- sample(seq(1,nrow(files)),1)
      test_file <- files[test_file_id,]
      training_files <- files[-test_file_id,]
      training_set <- bind_rows(training_set,training_files)
      test_set <- bind_rows(test_set, test_file)
    } 
  }
  colnames(training_set)[1] <- "name"
  training_set <- training_set[order(training_set$group),]
  test_set <- test_set[order(test_set$group),]  
  write.table(training_set, file = paste0(output_folder,"/sets/","training_files_",k,".csv"), quote = F,sep = ",",row.names = F,col.names = T)
  write.table(test_set, file = paste0(output_folder,"/sets/","test_files_",k,".csv"), quote = F,sep = ",",row.names = F,col.names = T)
}


