load("data.Rdata")
load("CITBCMST.rdata")

#########################################################
#########################################################
##################ADJUSTED DATA##########################
#########################################################
#########################################################
load("data.Rdata")
gene_sub <- read.csv("dtc_data.csv", sep=";")
rownames(gene_sub) <- gene_sub[[1]]
gene_sub <- gene_sub[, -1]


## leave-one-out classification of CIT subtyping (DTC)

pred_vector <- c()
for(i in 1:ncol(gene_sub)) {
  training <- gene_sub[,-i]
  training_classes <- CIT_classes[-i]
  test <- gene_sub[,i]
  test_class <- CIT_classes[i]
  
  centroids <- NULL
  for (class in unique(CIT_classes)) {
    class_centroid <- rowMeans(training[,training_classes==class])
    centroids <- cbind(centroids, class_centroid)
  }
  colnames(centroids) <- unique(CIT_classes)
  d <- as.matrix(dist(t(cbind(centroids, test))))
  class_pred <- names(which.min(d[1:6,7]))
  pred_vector <- c(pred_vector, test_class==class_pred)
}

table(pred_vector)

#########################################################
#########################################################
##################LUM C REMOVED ADJUSTED DATA############
#########################################################
#########################################################

gene_sub <- t(gene_sub)
gene_sub <- cbind(gene_sub, CIT_classes)
gene_sub <- data.frame(gene_sub)
gene_sub <- gene_sub[!gene_sub$CIT_classes == "lumC",]
gene_sub <- gene_sub[, -ncol(gene_sub)]
gene_sub <- t(gene_sub)
gene_sub <- data.frame(gene_sub)
gene_sub[, 1:307] <- apply(gene_sub[, 1:307], 2, as.numeric)

filtered_vector <- CIT_classes[names(CIT_classes) != "lumC"]
my_var <- gsub("lumC", NA, CIT_classes)
my_var <- na.omit(my_var)
CIT_classes <- my_var

table(my_var)
class(CIT_classes)
typeof(CIT_classes)
table(filtered_vector)
table(CIT_classes)
test <- data.frame(CIT_classes)

pred_vector <- c()
for(i in 1:ncol(gene_sub)) {
  training <- gene_sub[,-i]
  training_classes <- CIT_classes[-i]
  test <- gene_sub[,i]
  test_class <- CIT_classes[i]
  
  centroids <- NULL
  for (class in unique(CIT_classes)) {
    class_centroid <- rowMeans(training[,training_classes==class])
    centroids <- cbind(centroids, class_centroid)
  }
  colnames(centroids) <- unique(CIT_classes)
  d <- as.matrix(dist(t(cbind(centroids, test))))
  class_pred <- names(which.min(d[1:5,6]))
  pred_vector <- c(pred_vector, test_class==class_pred)
}

table(pred_vector)
