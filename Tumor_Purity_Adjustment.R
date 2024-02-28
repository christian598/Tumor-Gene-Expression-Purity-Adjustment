library(GSVA)
library(estimate)
library(ggplot2)
library(readxl)

#LOAD DATA
load("data.Rdata")

## ESTIMATE PURITY OF SAMPLES

ESTIMATE_fun <- function(expr) {
  write.table(as.data.frame(expr), file = "rma.data.gct", quote = F, sep = '\t')
  
  filterCommonGenes(input.f = 'rma.data.gct', output.f = "RMA_10412.gct", id = "GeneSymbol")
  estimateScore("RMA_10412.gct", "estimate_score.gct", platform = "affymetrix")
  estimate <- read.table("estimate_score.gct", sep = '\t', row.names = 1, header = T, skip = 2); estimate <- estimate[,-1]
  
  # Clean up
  file.remove("rma.data.gct"); 
  file.remove("RMA_10412.gct"); 
  file.remove("estimate_score.gct")
  
  # Process scores for easy output
  estimate_scores <- as.data.frame(t(estimate))
  colnames(estimate_scores) <- c('StromalSignature', 'ImmuneSignature', 'ESTIMATEScore', 'TumorPurity')
  
  return(estimate_scores)
}

estimate_scores <- ESTIMATE_fun(CIT_full)


load("CITBCMST.rdata")

## Combine TumorPurity and CIT_full
resData <- data.frame(cbind(t(CIT_full),estimate_scores$TumorPurity))
colnames(resData)[colnames(resData) == "V20698"] <- "TumorPurity"


## Calculate linear relationship ##

calculate_linear_relationship <- function(data_frame) {
  # Get the column names of the genes
  gene_columns <- names(data_frame)[1:(ncol(data_frame)-1)]
  
  # Create an empty list to store the linear relationships
  linear_relationships <- list()
  
  # Iterate over each gene column
  for (gene in gene_columns) {
    # Get the gene expression values and purity values as vectors
    gene_expression <- data_frame[, gene]
    purity <- data_frame[, "TumorPurity"]
    
    # Calculate the linear relationship between the gene and purity
    linear_model <- lm(gene_expression ~ purity)
    
    # Store the linear model in the list
    linear_relationships[[gene]] <- linear_model
  }
  
  return(linear_relationships)
}

linear_relationships <- calculate_linear_relationship(resData)

## Solve Purity for 50% ##

results <- c()

y <- c() 

name_list <- list()

length(linear_relationships)

for (i in 1:length(linear_relationships)) { 
  
  print(i) 
  
  pred_val <- coef(linear_relationships[[i]])[1] * 0.5 + coef(linear_relationships[[i]])[2]
  
  adj_val <- data.frame("adj_val"  = c(pred_val + linear_relationships[[i]]$residuals))
  
  name_list <- append(name_list, names(linear_relationships)[i])
  
  results <- c(results, adj_val) 
}

#TURN IT INTO DATAFRAME AND ADD PROBE NAMES AND SAMPLE NAMES
results <- data.frame(results)
colnames(results) <- name_list
rownames(results) <- rownames(resData)
results <- t(results)

#REMOVES X FROM THE FRONT OF THE NAME IN THE COLOUMN
new_colnames <- sub("^.", "", colnames(results))
colnames(results) <- new_colnames

####################################################

article_genes <- read_excel("subset.xls")

gene_sub <- results[, colnames(results) %in% article_genes$`Probe sets`]

gene_sub <- results[rownames(results) %in% rownames(CIT_subtyping),]

write.csv2(gene_sub, "dtc_data.csv")

####################################################
################# Data Wrangle #####################
####################################################

article_genes <- read_excel("subset.xls")
gene_sub <- read.csv("gene_sub.csv")
rownames(gene_sub) <- gene_sub[[1]]
gene_sub <- gene_sub[, -1]

new_colnames <- sub("^.", "", colnames(gene_sub))
colnames(gene_sub) <- new_colnames

gene_sub <- cbind(gene_sub, CIT_classes)

means <- aggregate(gene_sub[, -ncol(gene_sub)], by = list(gene_sub$CIT_classes), FUN = mean)

colnames(means)[1] <- "CIT_classes"
rownames(means) <- means[[1]]
means <- means[, -1]
means <- t(means)

gene_order <- article_genes$"Probe sets"

gene_names <- rownames(means)

df2_reordered <- means[gene_order, ]

rownames(df2_reordered) <- gene_order

str(df2_reordered)
df2_reordered <- data.frame(df2_reordered)
df2_reordered[, 1:6] <- apply(df2_reordered[, 1:6], 2, as.numeric)

df2_reordered <- cbind(df2_reordered, article_genes$`Gene cluster`)
df2_reordered <- cbind(df2_reordered, article_genes$`num cluster`)

colnames(df2_reordered_correct)[7] <- "Gene_cluster"
colnames(df2_reordered)[7] <- "num_cluster"

####################################################
################## PCA #############################
####################################################
gene_sub <- gene_sub[, -ncol(gene_sub)]

df2_reordered <- df2_reordered[, -ncol(df2_reordered)]

pca_count <- prcomp(df2_reordered)
pca_count$x

summary(pca_count)

PCA <- data.frame(PC1 = pca_count$x[, 1], PC2 = pca_count$x[, 2])

ggplot(PCA, aes(x = PC1, y = PC2, color = article_genes$`Gene cluster`)) +
  geom_point()


####################################################
################## BOXPLOT #########################
####################################################
adj_genes <- read.csv2("adj_genes.csv")
colnames(adj_genes)[1] <- "CIT_classes"
rownames(adj_genes) <- adj_genes[[1]]
adj_genes <- adj_genes[, -1]

bplot_data <- data.frame(df2_reordered)
bplot_data

cluster1 <- bplot_data[bplot_data$`num_cluster` %in% "I",] 

cluster2 <- bplot_data[bplot_data$`num_cluster` %in% "II",] 

cluster3 <- bplot_data[bplot_data$`num_cluster` %in% "III",] 

cluster4 <- bplot_data[bplot_data$`num_cluster` %in% "IV",] 

cluster5 <- bplot_data[bplot_data$`num_cluster` %in% "V",] 

cluster6 <- bplot_data[bplot_data$`num_cluster` %in% "VI",] 

cluster7 <- bplot_data[bplot_data$`num_cluster` %in% "VII",] 

cluster8 <- bplot_data[bplot_data$`num_cluster` %in% "VIII",] 

cluster9 <- bplot_data[bplot_data$`num_cluster` %in% "IX",] 


cluster1_longer <-  cluster1 %>% pivot_longer(cols = -c("num_cluster"), values_transform = as.numeric) 
plot1 <- cluster1_longer %>% ggplot(aes(x = name, y = value))+ 
  geom_boxplot(aes(fill = name)) +
  labs(title = "basL")

cluster2_longer <-  cluster2 %>% pivot_longer(cols = -c("num_cluster"), values_transform = as.numeric) 
plot2 <- cluster2_longer %>% ggplot(aes(x = name, y = value))+ 
  geom_boxplot(aes(fill = name)) +
  labs(title = "cycle/proliferation")

cluster3_longer <-  cluster3 %>% pivot_longer(cols = -c("num_cluster"), values_transform = as.numeric) 
plot3 <- cluster3_longer %>% ggplot(aes(x = name, y = value))+ 
  geom_boxplot(aes(fill = name)) +
  labs(title = "lumA")

cluster4_longer <-  cluster4 %>% pivot_longer(cols = -c("num_cluster"), values_transform = as.numeric) 
plot4 <- cluster4_longer %>% ggplot(aes(x = name, y = value))+ 
  geom_boxplot(aes(fill = name)) +
  labs(title = "AR")

cluster5_longer <-  cluster5 %>% pivot_longer(cols = -c("num_cluster"), values_transform = as.numeric) 
plot5 <- cluster5_longer %>% ggplot(aes(x = name, y = value))+ 
  geom_boxplot(aes(fill = name)) +
  labs(title = "lumB")

cluster6_longer <-  cluster6 %>% pivot_longer(cols = -c("num_cluster"), values_transform = as.numeric) 
plot6 <- cluster6_longer %>% ggplot(aes(x = name, y = value))+ 
  geom_boxplot(aes(fill = name)) +
  labs(title = "ESR1")

cluster7_longer <-  cluster7 %>% pivot_longer(cols = -c("num_cluster"), values_transform = as.numeric) 
plot7 <- cluster7_longer %>% ggplot(aes(x = name, y = value))+ 
  geom_boxplot(aes(fill = name)) +
  labs(title = "lumA/normL")

cluster8_longer <-  cluster8 %>% pivot_longer(cols = -c("num_cluster"), values_transform = as.numeric) 
plot8 <- cluster8_longer %>% ggplot(aes(x = name, y = value))+ 
  geom_boxplot(aes(fill = name)) +
  labs(title = "mApo")

cluster9_longer <-  cluster9 %>% pivot_longer(cols = -c("num_cluster"), values_transform = as.numeric) 
plot9 <- cluster9_longer %>% ggplot(aes(x = name, y = value))+ 
  geom_boxplot(aes(fill = name)) +
  labs(title = "normL")

library(cowplot)

plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, ncol = 3)
