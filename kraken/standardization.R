library(xlsx)
library(tidyr)
library(dplyr)
library(rlang)


db <- read.delim("/Volumes/T7_Shield/Cyril/006_kraken/plufplf_db/inspect.txt", header=F)
db <- db %>%
  filter(V4=="S")%>%
  mutate(V6 = trimws(V6)) %>%
  select(c(V5,V6))
colnames(db) <- c("taxid","name")

meta <- read.xlsx("meta.xlsx",sheetIndex = 1)
ids <- meta$id

for (id in ids){
  lambda <- meta[meta$id==id,]$nr_lambda_reads
  data <- read.delim(paste0("kraken2_results/", id, "_report.kraken2"), header = F)
  data <- data %>%
    filter(V6 == "S") %>%
    mutate(V8 = trimws(V8)) %>%
    #filter(V5>V4/2) %>% #remove FP based on kmers_unique
    mutate(across(V2:V5, ~ .x*4000/lambda)) #normalization based on lambda count
  
  assign(id, data)
  
  data_clean <- data %>%
    select(c(V2,V7))
  colnames(data_clean) <- c(as.character(id),"taxid")
  
  assign(paste0(id,"_clean"), data_clean)
  
  db <- merge(db, data_clean, all=T, by="taxid")
  
}
db <- db %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) #%>%
# mutate(across(3:47, ~ ifelse(. < 3, 0, .))) 

db_t <- t(db)
db_t <- as.data.frame(db_t)
colnames(db_t) <- trimws(db_t["name",])
db_t <- db_t %>%
  filter(row_number()>2) 

db_t <- db_t %>%
  mutate(across(where(is.character), as.numeric))
# test <- db_t %>% summarize(across(everything(), max))

df_cleaned <- db_t %>%
  select(where(~ max(.) != 0))


colnames(df_cleaned) <- stringr::str_replace(colnames(df_cleaned), " ", "_")
df_cleaned$id <- rownames(df_cleaned)

species <- colnames(df_cleaned)
species <- species[-which(species=="id")]

df_cleaned <- merge(df_cleaned,meta, by="id")
rownames(df_cleaned) <- NULL
df_corradj <- df_cleaned
for (bug in species){
  for (n_batch in unique(df_cleaned$batch)){
    data <- filter(df_cleaned, batch == n_batch)
    n_obs <- nrow(data)
    
    if (max(data[[bug]])<50){
      df_corradj <- df_corradj %>%
        mutate(!!sym(bug) := ifelse(batch == n_batch, 0, !!sym(bug)))
      break
    }
    
    if (sum(data[[bug]]==0) < round(0.5*n_obs)){ #present in 50% of the samples
      corr <- cor(1/data$conc, data[[bug]])
      Q1 <- quantile(data[[bug]], 0.25)
      Q3 <- quantile(data[[bug]], 0.75)
      IQR <- Q3-Q1
      sd <- sd(data[[bug]])
      mean <- mean(data[[bug]])
      median <- median(data[[bug]])
      
      if (corr > 0.5){ #pre-dilution contaminant
        print(paste(bug, "has a correlation of ",corr, "so was removed in batch", n_batch))
        df_corradj <- df_corradj %>%
          mutate(!!sym(bug) := ifelse(batch == n_batch, 0, !!sym(bug)))
      }
      else if (max(data[[bug]]) >= Q3+1.5*IQR){ # if one sample is an outlier then proceed to next batch, we don't remove it
        print(paste("An outlier in bug", bug, "was detected: ", max(data[[bug]]), ">=", Q3+1.5*IQR))
        df_corradj <- df_corradj %>%
          mutate(!!sym(bug) := ifelse(batch == n_batch, !!sym(bug)-median, !!sym(bug)))
      }
      else if (sum(data[[bug]]!=0) == n_obs) { #post-dilution contaminant in all samples
        
        df_corradj <- df_corradj %>%
          mutate(!!sym(bug) := ifelse(batch == n_batch, 0, !!sym(bug)))
      }
    }
  }
}

if (!exists("backup")){
  backup <- df_corradj
}

all_species <- species

#major part of sensitivity/specificity is played here
df_corradj <- backup %>%
  mutate(across(all_species, ~ ifelse(. < 15, 0, .))) %>%
  select(where(~ max(.) > 49)) 


df_corradj$dilution <- ifelse(df_corradj$conc>0.5, df_corradj$conc/0.5, 1)

species <- all_species[all_species %in% names(df_corradj)]
test <- df_corradj %>%
  mutate(across(species, ~ . * dilution)) %>%
  mutate(across(species, ~ . * 500/plasma_volume))

a <- test %>%
  select(c(id,all_of(species))) %>%
  mutate(across(species, scale))


a <- as.data.frame(t(a))
colnames(a) <- paste0(a[1,],"_z")
a$species <- rownames(a)
a <- a[-1,]
rownames(a) <- NULL

test_t <- test %>%
  select(c(id,all_of(species)))
test_t <- as.data.frame(t(test_t))
colnames(test_t) <- paste0(test_t[1,])
test_t$species <- rownames(test_t)
test_t <- test_t[-1,]
rownames(test_t) <- NULL

a <- merge(a, test_t, by="species")

a <- a %>% select(order(colnames(a)))
a <- a %>%
  mutate(across(-species, ~ as.numeric(.)))

df <- data.frame(NULL)
for (id in ids){
  data <- a %>%
    select(c(starts_with(id), species)) %>%
    filter(.[[id]]>=500 | .[[paste0(id,"_z")]]>=4)
  if (nrow(data)>0){
    data$id <- id
    colnames(data) <- c("reads","z","specie", "id")
    df <- bind_rows(df,data)
  } 
  else{
    print(id)
  }
  
}
write.csv(df,file = "test.csv")








for (id in ids){
  lambda <- meta[meta$id==id,]$nr_lambda_reads
  data <- read.delim(paste0("kraken2_results/", id, "_report.kraken2"), header = F)
  data <- data %>%
    mutate(across(V2:V5, ~ ifelse(stringr::str_replace(trimws(V8), " ", "_") %in% species, as.numeric(a %>% filter(species == stringr::str_replace(trimws(V8), " ", "_")) %>% select(id)), 0)))
  
  write.table(data, file = paste0("kraken2_results/cleaned/",id,".tsv"), quote=F, sep="\t", row.names = F, col.names = F)
  
}


data <- data %>%
  mutate(V8 = stringr::str_replace(trimws(V8), " ", "_"))

merge(data, a, by.x="V8", by.y="species")

