#######################################
install.packages(c("tidyverse","fastICA","factoextra","skimr","Amelia","psych"))

library(tidyverse)
library(fastICA)
library(factoextra)
set.seed(21)
theme_set(theme_minimal())

## IMPORT DATA
DB <-read.csv("https://raw.githubusercontent.com/ProfNascimento/HighDim/main/harvesters.csv")
names(DB)
DB$Vehicle=as.factor(DB$Vehicle)

## DESCRIPTIVE 
skim_summ <- skimr::skim_with(base = skimr::sfl())
skim_summ(DB)
Amelia::missmap(DB)
dev.off()

## UNSUPERVISED LEARNING
## FEATURES ANALYSIS (Xs ONLY)
data_features <- DB[,-c(2:3,12)] %>% 
  dplyr::select(-Vehicle) %>% 
  mutate(across(everything(), ~ scale(.x)))

apply(data_features,2,summary)

cor(data_features) %>% 
  as_tibble() %>% 
  mutate(feature_y = names(.)) %>% 
  pivot_longer(cols = -feature_y, names_to = "feature_x", values_to = "correlation") %>% 
  mutate(feature_y = fct_rev(feature_y)) %>% 
  ggplot(aes(x = feature_x, y = feature_y, fill = correlation)) + 
  geom_tile() + 
  labs(x = NULL,y = NULL) +
  scico::scale_fill_scico(palette = "cork", limits = c(-1,1)) + 
  coord_fixed() + 
  theme(text = element_text(size = 30),
        axis.text.x = element_text(hjust = 1, angle = 30),
        legend.title=element_blank())

###############################
## ----------PCA------------ ##
###############################
pc_model <- prcomp(data_features,center = TRUE, scale. = TRUE)
(pc_model$variance <- pc_model$sdev^2)
#VARIANCE CONTRIBUTION
pc_model$variance/sum(pc_model$variance)
#CUMULATIVE CONTRIBUTION
cumsum(pc_model$variance/sum(pc_model$variance))

pc_model$variance %>% 
  as_tibble() %>%
  rename(eigenvalue = value) %>% 
  rownames_to_column("comp") %>% 
  mutate(comp = parse_number(comp),
         cum_variance = cumsum(eigenvalue)/sum(eigenvalue)) %>% 
  ggplot(aes(x = comp, y = eigenvalue)) + 
  geom_hline(yintercept = 1,lty=2) +
  geom_line(size = 1) + 
  geom_point(size = 3)

## NUMBERS OF COMPONENTS
n_comps=3

# PCA LOADING (EIGENVALUES)
pc_weight_matrix <- pc_model$rotation %>% 
  data.frame() %>% 
  rownames_to_column("variable") %>% 
  pivot_longer(starts_with("PC"), names_to = "prin_comp", values_to = "loading")

## PCA VISUALIZATION
fviz_pca_var(pc_model, col.var="contrib")+
  scale_color_gradient2(low="#56B4E9", mid="orange", 
                        high="#E69F00", midpoint=10)+theme_bw()

fviz_pca_biplot(pc_model,
                select.ind = list(contrib = 5))

fviz_pca_ind(pc_model, label="none", habillage=DB$Activity,
             addEllipses=TRUE, ellipse.level=0.7, alpha=0.05)

## Orthogonal (INDEPENDENT) SPACE generated from the PCA obtained dimensions
pc_model$x %>% 
  cor() %>% 
  data.frame() %>% 
  rownames_to_column("comp_x") %>% 
  pivot_longer(cols = starts_with("PC"), names_to = "comp_y", values_to = "correlation") %>% 
  filter(parse_number(comp_x) <= n_comps,
         parse_number(comp_y) <= n_comps) %>% 
  ggplot(aes(x = comp_x, y = comp_y, fill = correlation)) + 
  geom_tile() +
  geom_text(aes(label = round(correlation,4)), color = "white") +
  labs(title = "Correlation between PCs",
       x = NULL,
       y = NULL) +
  scico::scale_fill_scico(palette = "berlin", limits = c(-1,1)) + 
  coord_equal()

# SELECTING OBS (#4551)
as.matrix(data_features[4551,])%*%pc_model$rotation

###############################
## ----------ICA------------ ##
###############################
ica_model <- fastICA(data_features, n.comp = n_comps)
# A list containing the following components
# X	- pre-processed data matrix
# K	- pre-whitening matrix that projects data onto the first n.comp principal components.
# W	- estimated un-mixing matrix (see definition in details)
# A	- estimated mixing matrix
# S	- estimated source matrix

ica_weight_matrix <- data.frame(t(ica_model$A)) %>% 
  rename_with(~ str_glue("IC{seq(.)}")) %>%
  mutate(variable = names(data_features)) %>%
  pivot_longer(cols = starts_with("IC"), names_to = "ic", values_to = "loading")

# MIXED SIGNALS (ORIGINAL DATA)
dim(ica_model$X)
# ICA SOURCE ESTIMATION
dim(ica_model$S)

# VISUALIZING 
par(mfrow=c(3,1))
plot(ica_model$S[1:100,1], main="ICA components #1",type="l")
plot(ica_model$S[1:100,2], main="ICA components #2",type="l")
plot(ica_model$S[1:100,3], main="ICA components #3",type="l")

# PCA ESTIMATION
ica.pca=ica_model$X %*% ica_model$K
dim(ica.pca)
par(mfrow=c(3,1))
plot(ica.pca[1:100,1], main="PCA components",type="l")
plot(ica.pca[1:100,2], main="PCA components",type="l")
plot(ica.pca[1:100,3], main="PCA components",type="l")

## COMPARISON 3 -COMPONENTS (PCA vs ICA)
pairs(ica_model$S[1:100,], main="ICA Pre-processed DATA")
pairs(pc_model$x[1:100,1:3], main="PCA Pre-processed DATA")

# SELECTING OBS (#4551)
as.matrix(data_features[4551,])%*%t(ica_model$A)

###############################
## -----FACTOR ANALYSES----- ##
###############################
require(psych)
fa_model <- fa(data_features, nfactors = 3)
#Getting the factor loadings and model analysis
fa_model

fa_weight_matrix <- fa_model$loadings[] %>% 
  data.frame() %>% 
  rownames_to_column("variable") %>% 
  pivot_longer(starts_with("MR"), names_to = "factor", values_to = "loading")

fa_model$score.cor %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +
  geom_text(aes(label = round(value,4)), color = "white") +
  labs(title = "Correlation between FA scores",
       x = NULL,
       y = NULL) +
  coord_equal()

## NEW FEATURES BASED ON THE ADJUSTED FA
pairs(fa_model$scores[1:100,])

# LOADINGS COMPARISON (PCA-ICA-FA)
all_weight_matrices <- bind_rows(
  pc_weight_matrix %>% 
    rename(comp = prin_comp) %>% 
    mutate(alg = "PCA"), 
  ica_weight_matrix %>% 
    rename(comp = ic) %>% 
    mutate(alg = "ICA"),
  fa_weight_matrix %>% 
    rename(comp = factor) %>% 
    mutate(alg = "FA")
)

all_weight_matrices %>% 
  filter(parse_number(comp) <= n_comps) %>% 
  mutate(alg = str_glue("{alg} loadings"),
         alg = as_factor(alg)) %>% 
  ggplot(aes(x = comp, y = variable, fill = loading)) +
  geom_tile() + 
  labs(x = NULL,
       y = NULL) + 
  scico::scale_fill_scico(palette = "cork", limits = c(-1,1)) + 
  facet_wrap(~ alg, scales = "free_x")
