###Edward Gilbert R script GLMMTMB

#---------------------load and sort data--------------------#
library(dplyr)
TL<-read.csv("TL_temporal_transformed.csv")
TL$Year <- as.factor(TL$Year)
TL$Env.type <- as.factor(TL$Env.type)
TL$sex <- as.factor(TL$sex)
TL$tissue_type <- as.factor(TL$tissue_type)
TL$morphotype <- as.factor(TL$morphotype)
TL$shifted_zscore <- TL$zscore - min(TL$zscore, na.rm = TRUE) + 1
TL_clean <- TL %>% filter(!is.na(SVL) & !is.na(sex) & !is.na(shifted_zscore))

#---------------PCA for environmental variable---------------------#

##Select environmental variables only
env_vars <- TL_clean[, c(
  "elevation",
  "RH.mean",
  "wind.mean",
  "sol_rad",
  "TSKY.mean"
)]

#Run PCA (scale = TRUE is essential)
env_pca <- prcomp(env_vars, center = TRUE, scale. = TRUE)

#proportion of variance explained
summary(env_pca)

pca_variance <- data.frame(
  PC = paste0("PC", 1:length(env_pca$sdev)),
  Variance_Explained = env_pca$sdev^2 / sum(env_pca$sdev^2),
  Cumulative_Variance = cumsum(env_pca$sdev^2 / sum(env_pca$sdev^2))
)

pca_variance

plot(
  pca_variance$Variance_Explained,
  type = "b",
  xlab = "Principal Component",
  ylab = "Proportion of Variance Explained"
)

biplot(
  env_pca,
  choices = c(1, 2),
  scale = 0,
  cex = 0.7
)

library(ggplot2)

scores_df <- as.data.frame(env_pca$x[, 1:2])
scores_df$ID <- rownames(scores_df)

loadings_df <- as.data.frame(env_pca$rotation[, 1:2])
loadings_df$Variable <- rownames(loadings_df)


#Extract PC scores
TL_clean$Env_PC1 <- env_pca$x[, 1]
TL_clean$Env_PC2 <- env_pca$x[, 2]

var_explained <- (env_pca$sdev^2) / sum(env_pca$sdev^2)

pc1_lab <- paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")

ggplot(scores_df, aes(x = PC1, y = PC2)) +
  geom_point(
    size = 2,
    alpha = 0.6,
    colour = "grey30"
  ) +
  geom_segment(
    data = loadings_df,
    aes(
      x = 0, y = 0,
      xend = PC1 * 3,
      yend = PC2 * 3
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.7,
    colour = "black"
  ) +
  geom_text(
    data = loadings_df,
    aes(
      x = PC1 * 3.2,
      y = PC2 * 3.2,
      label = Variable
    ),
    size = 4,
    family = "Times"
  ) +
  xlab(pc1_lab) +
  ylab(pc2_lab) +
  theme_bw(base_size = 12, base_family = "Times") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.8),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

#extract and report loadings
env_pca$rotation

loadings_table <- data.frame(
  Variable = rownames(env_pca$rotation),
  PC1 = env_pca$rotation[, 1],
  PC2 = env_pca$rotation[, 2],
  PC3 = env_pca$rotation[, 3]
)

loadings_table


###-------------------GLMM--------------------------###

library(glmmTMB)

TLmod <- glmmTMB(
  shifted_zscore ~ Env.type + Env_PC1 + Env_PC2 + morphotype + tissue_type + sex + SVL +
    (1 | Year),
  family = Gamma(link = "log"),
  data = TL_clean
)

#install.packages("MuMIn")
library(MuMIn)

#set options so dredge works with glmmTMB
options(na.action = "na.fail")

#run dredge (random effect always included)
dredge_results <- dredge(
  TLmod,
  fixed = c("1", "(1 | Year)") # Ensures random effect is always present
)

head(dredge_results)

#model averaging
#close AICc and weights estimates effect sizes and confidence
#intervals for predictors across top models
#Gives averaged coefficients weighted by model support
avg_model <- model.avg(dredge_results, subset = delta < 2)

summary(avg_model)

###-----------------make plots and emmeans------------------###
library(emmeans)

#refit global model before dredging
#Figures show predictions from the global mixed-effects model; 
#inference is based on model-averaged coefficients
TLmod_global <- glmmTMB(
  shifted_zscore ~ Env_PC1 + Env_PC2 + SVL * sex +
    morphotype + tissue_type +
    (1 | Year),
  family = Gamma(link = "log"),
  data = TL_clean
)


##environmental PC1 — emmeans
emm_pc1 <- emmeans(
  TLmod_global,
  ~ Env_PC1,
  at = list(
    Env_PC1 = seq(
      min(TL_clean$Env_PC1),
      max(TL_clean$Env_PC1),
      length.out = 100
    )
  ),
  type = "response"
)

emm_pc1_df <- as.data.frame(emm_pc1)

#manually compute confidence limits
emm_pc1_df$asymp.LCL <- emm_pc1_df$response -
  qnorm(0.975) * emm_pc1_df$SE

emm_pc1_df$asymp.UCL <- emm_pc1_df$response +
  qnorm(0.975) * emm_pc1_df$SE

ggplot() +
  geom_jitter(
data = TL_clean, aes(x = Env_PC1, y = shifted_zscore, colour = Env.type), alpha = 0.6,width = 0.2,
    size = 4) + geom_line(data = emm_pc1_df,aes(x = Env_PC1, y = response),colour = "black",linewidth = 1) +
  geom_ribbon(data = emm_pc1_df,aes(x = Env_PC1,ymin = asymp.LCL,ymax = asymp.UCL),inherit.aes = FALSE,
    alpha = 0.2) +labs(
    x = "Environmental PC1",
    y = "Shifted z-score") +
  scale_color_brewer(palette = "Set2", name = "Environment type") +
  theme_classic(base_size = 16)

#extract model summary
mod_sum <- summary(TLmod_global)

#slope (log scale)
pc1_slope <- mod_sum$coefficients$cond["Env_PC1", "Estimate"]

#standard error
pc1_se <- mod_sum$coefficients$cond["Env_PC1", "Std. Error"]

#wald z-value
pc1_z <- mod_sum$coefficients$cond["Env_PC1", "z value"]

#p-value
pc1_p <- mod_sum$coefficients$cond["Env_PC1", "Pr(>|z|)"]
pc1_slope
pc1_p


###--------------sex and svl

#recompute emmeans
emm_svl <- emmeans(TLmod_global,~ SVL | sex, at = list(
    SVL = seq(min(TL_clean$SVL), max(TL_clean$SVL), length.out = 50)),type = "response")

emm_svl_df <- as.data.frame(emm_svl)

#compute confidence limits manually from SE
emm_svl_df$asymp.LCL <- emm_svl_df$response -
  qnorm(0.975) * emm_svl_df$SE

emm_svl_df$asymp.UCL <- emm_svl_df$response +
  qnorm(0.975) * emm_svl_df$SE

#plot using asymp.LCL / asymp.UCL
ggplot(TL_clean, aes(x = SVL, y = shifted_zscore, colour = sex)) +
  geom_jitter(alpha = 0.4, width = 0.2) +
  
  geom_line(
    data = emm_svl_df,
    inherit.aes = FALSE,
    aes(x = SVL, y = response, colour = sex),
    linewidth = 1
  ) +
  
  geom_ribbon(
    data = emm_svl_df,
    inherit.aes = FALSE,
    aes(x = SVL, ymin = asymp.LCL, ymax = asymp.UCL, fill = sex),
    alpha = 0.2,
    colour = NA
  ) +
  
  labs(
    x = "SVL",
    y = "Predicted shifted z-score"
  ) +
  theme_classic(base_size = 16)

mod_sum <- summary(TLmod_global)$coefficients$cond
svl_slope_ref <- mod_sum["SVL", "Estimate"]
svl_p_ref     <- mod_sum["SVL", "Pr(>|z|)"]

svl_slope_ref
svl_p_ref
#male
svl_slope_male <- svl_slope_ref + mod_sum["SVL:sexM", "Estimate"]
svl_p_male     <- mod_sum["SVL:sexM", "Pr(>|z|)"]

#female
svl_slope_female <- mod_sum["SVL", "Estimate"]
svl_p_female     <- mod_sum["SVL", "Pr(>|z|)"]

svl_slope_male
svl_p_male

svl_slope_female
svl_p_female

levels(TL_clean$sex)

pairs(emmeans(TLmod_global, ~ sex | SVL))

library(viridis)


##------------------Tissue type

#predicted means and SE
emm_tissue <- emmeans(
  TLmod_global,
  ~ tissue_type,
  type = "response")

#pairwise comparisons for manual annotation
pairs(emm_tissue, adjust = "tukey")

#convert to dataframe for plotting
emm_tissue_df <- as.data.frame(emm_tissue)

#plot predicted means with error bars
ggplot(TL_clean, aes(x = tissue_type, y = shifted_zscore, fill = tissue_type)) +
  geom_jitter(aes(color = Year, shape = sex), alpha = 0.6, width = 0.2, size = 3) +
  geom_errorbar(
    data = emm_tissue_df,
    inherit.aes = FALSE,
    aes(x = tissue_type, ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.2, color = "black"
  ) +
  geom_point(
    data = emm_tissue_df,
    inherit.aes = FALSE,
    aes(x = tissue_type, y = response),
    color = "black", size = 3
  ) +
  labs(x = "Tissue Type", y = "Predicted shifted z-score") +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  theme_classic(base_size = 16)

tissue_effects <- mod_sum[grep("^tissue_type", rownames(mod_sum)), ]

tissue_effects

###--------------------------Sex

#predicted means and SE
emm_sex <- emmeans(
  TLmod_global,
  ~ sex,
  type = "response")

#pairwise comparisons for manual annotation
pairs(emm_sex, adjust = "tukey")

#convert to dataframe for plotting
emm_sex_df <- as.data.frame(emm_sex)

#plot predicted means with error bars
ggplot(TL_clean, aes(x = sex, y = shifted_zscore, fill = sex)) +
  geom_jitter(aes(color = Year, shape = tissue_type), alpha = 0.6, width = 0.2, size = 3) +
  geom_errorbar(
    data = emm_sex_df,
    inherit.aes = FALSE,
    aes(x = sex, ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.2, color = "black"
  ) +
  geom_point(
    data = emm_sex_df,
    inherit.aes = FALSE,
    aes(x = sex, y = response),
    color = "black", size = 3
  ) +
  labs(x = "Sex", y = "Predicted shifted z-score") +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  theme_classic(base_size = 16)

sex_effects <- mod_sum[grep("^sex", rownames(mod_sum)), ]

sex_effects


###-----------intrinsic vs extrinsic variables--------------###
library(effectsize)
library(r2glmm)
library(hier.part)

#Intrinsic predictors
intrinsic_vars <- c("tissue_type", "sex", "SVL")

#Extrinsic predictors
extrinsic_vars <- c("Year", "Env_PC1", "Env_PC2")

#Full model (already fitted)
TLmod_global <- glmmTMB(shifted_zscore ~ tissue_type + sex*SVL + Year + Env_PC1 + Env_PC2 + (1 | Year),
                        family = Gamma(link = "log"), data = TL_clean)

##--------Standardised effect sizes
## "refit" ensures predictors are standardised before model fitting,
##which is the correct approach for GLMMs

std_effects <- standardize_parameters(
  TLmod_global,
  method = "refit"
)

#convert to data frame
std_df <- as.data.frame(std_effects)

#remove intercept
std_df <- std_df %>%
  filter(Parameter != "(Intercept)")

#clean labels for plotting
std_df$Parameter <- factor(
  std_df$Parameter,
  levels = std_df$Parameter[order(std_df$Std_Coefficient)]
)

#plot
ggplot(std_df, aes(x = Parameter, y = Std_Coefficient)) +geom_point(size = 3) +geom_errorbar(
    aes(ymin = CI_low, ymax = CI_high),width = 0.2) +coord_flip() +
  labs(
    x = "Predictor",
    y = "Standardised effect size",
    title = "Standardised effect sizes from GLMM",
    subtitle = "Gamma GLMM with environmental PCA and Year as random effect") +theme_bw(base_size = 16)

#fit intrinsic-only GLMM
int_mod <- glmmTMB(
  shifted_zscore ~ tissue_type + sex*SVL + (1 | Year),
  family = Gamma(link = "log"),
  data = TL_clean
)

#fit extrinsic-only GLMM
ext_mod <- glmmTMB(
  shifted_zscore ~ Env_PC1 + Env_PC2 + (1 | Year),
  family = Gamma(link = "log"),
  data = TL_clean
)

#Pseudo-R2 for GLMMs
library(performance)

full_r2 <- r2(TLmod_global)   #Marginal and conditional
int_r2  <- r2(int_mod)
ext_r2  <- r2(ext_mod)

full_r2
int_r2
ext_r2

#Deviance explained
#compute deviance explained relative to null model
null_model <- glmmTMB(shifted_zscore ~ 1 + (1 | Year),
                      family = Gamma(link = "log"), data = TL_clean)

dev_explained_full <- 1 - (deviance(TLmod_global) / deviance(null_model))
dev_explained_int  <- 1 - (deviance(int_mod) / deviance(null_model))
dev_explained_ext  <- 1 - (deviance(ext_mod) / deviance(null_model))

dev_explained_full
dev_explained_int
dev_explained_ext

###--------------------Hierarchical partitioning----------------------------###

#For GLMM, hier.part supports glm; for glmmTMB, we need to drop random effects
#Use fixed effects only (approximation)
predictors <- TL_clean[, c("Env_PC1", "Env_PC2", "tissue_type", "sex", "SVL")]

hp <- hier.part(TL_clean$shifted_zscore,
                xcan = predictors,
                family = Gamma(link = "log"))

effect_data <- data.frame(
  predictor = rownames(hp$I.perc),
  percent_effect = hp$I.perc[,1]
)

ggplot(effect_data, aes(x = reorder(predictor, percent_effect), y = percent_effect)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # To make it horizontal
  theme_classic(base_size = 16) +
  labs(x = "Predictor Variable", y = "% Independent Effect", 
       title = "Independent Effects of Predictors")

hp_percent <- data.frame(
  Predictor = rownames(hp$I.perc),
  Independent_contribution_percent = hp$I.perc[, 1]
)

#randomisation test (approximate)
set.seed(123)
rand_test <- rand.hp(TL_clean$shifted_zscore, predictors,
                     family = Gamma(link = "log"),
                     gof = "logLik", num.reps = 100)

rand_test

