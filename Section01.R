require(mp)
require(tidyverse)

## DATA IMPORT
DB <- read.csv("https://raw.githubusercontent.com/ProfNascimento/ECU/refs/heads/main/wine.csv")
names(DB)
DB$Type <- as.factor(DB$Type)

## DESCRIPTIVE
skim_summ <- skimr::skim_with(base = skimr::sfl())
skim_summ(DB)

## PAIRWISE FEATURE VISUALIZATION
library(GGally)
ggpairs(DB)

## CORRELATION
apply(DB[,-1],2,sd) # Tl sd=0
corr <- round(cor(DB[,-1]), 1)
ggcorrplot::ggcorrplot(corr)

######################################
## MATRIX SIMILARITY (DIST.)
library(vegan)
JAC.dist <- vegdist(DB[,-1], method = "jaccard")

## MULTIDIMENSIONAL PROJECTION (MP)
# Force Scheme
emb <- forceScheme(dist(DB[,-1], method = "euclidean"))
plot(emb, col=DB[,1], xlab ="", ylab ="", main="")

emb.FS <- as.data.frame(cbind(emb,DB[,1]))
colnames(emb.FS) <- c("X1","X2","Y")

ggplot( emb.FS, aes( x=X1, y=X2, color=as.factor(Y) ) )+
  geom_point() + theme(legend.position="top") + 
  scale_color_discrete(name = "WINE TYPE", labels = c("A", "B", "C"))

## Multidimensional Scaling (MDS)
emb1 <- cmdscale(dist(DB[,-1], method = "euclidean"))
plot(emb1, col=DB[,1])

## Local Affine Multidimensional Projection (LAMP)
emb2 <- lamp(DB[,-1])
plot(emb2, col=DB[,1])

# Least Square Projection (LSP)
emb3 <- lsp(DB[,-1])
plot(emb3, col=DB[,1])

# Pekalska’s approach to speeding up Sammon’s mapping
emb4 <- pekalska(dist(DB[,-1], method = "euclidean"), sample.indices = NULL, Ys = NULL)
plot(emb4, col=DB[,1])

# Part-Linear Multidimensional Projection (PLMP)
emb5 <- plmp(DB[,-1])
plot(emb5, col=DB[,1])

# t-Distributed Stochastic Neighbor Embedding (t-SNE)
emb6 <- tSNE(as.matrix(DB[,-1]))
plot(emb6, col=DB[,1])

# Uniform Manifold Approximation and Projection (UMAP)
library(umap)
emb7 <- umap(DB[,-1])
plot(emb7$layout, col=DB[,1])

# Polar Swarm (Pswarm)
emb8 <- ProjectionBasedClustering::PolarSwarm(as.matrix(DB[,-1]))
plot(emb8$ProjectedPoints, col=DB[,1])