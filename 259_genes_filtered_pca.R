df_unique <- article_genes[!duplicated(article_genes$`Gene symbol`), ]
df_filtered <- CIT_core[rownames(CIT_core) %in% df_unique$`Probe sets`, ]


####Data wrangling###
gene_order <- df_unique$"Probe sets"

gene_names <- rownames(df_filtered)

df2_reordered <- df_filtered[gene_order, ]

rownames(df2_reordered) <- gene_order

df2_reordered <- cbind(df2_reordered, df_unique$`num cluster`)

df2_rendered <- data.frame(df2_reordered)

str(df2_rendered)

df2_rendered[, 1:355] <- apply(df2_rendered[, 1:355], 2, as.numeric)

###Col Mean###

means <- aggregate(df2_rendered[, -ncol(df2_rendered)], by = list(df2_rendered$V356), FUN = mean)

rownames(means) <- means[[1]]
means <- means[, -1]
means <- t(means)

means <- cbind(means, CIT_classes)
means <- means1


colnames(means)[1] <- "I"
rownames(means) <- means[[1]]
means <- means[, -1]

means <- means[, -ncol(means)]
means <- data.frame(means)
str(means)
means[] <- lapply(means, as.numeric)

###Perform PCA###

pca_count <- prcomp(means)
pca_count$x

summary(pca_count)

PCA <- data.frame(PC1 = pca_count$x[, 1], PC2 = pca_count$x[, 2])

ggplot(PCA, aes(x = PC1, y = PC2, color = means$CIT_classes)) +
  geom_point()















