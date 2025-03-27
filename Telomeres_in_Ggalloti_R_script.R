###Modelling telomere dynamics in Ggalloti###
###Edward Gilbert

###-----------raw CT-values-------------###
###Script by Megan Power
library(rptR)
data=read.table ("CT_values_All.csv", header=T, sep=',')
for(i in 1:nrow(data))
  
  #The following calculates the Pfaffl (2001) calculation for calculating relative
  #telomere length while taking into account any potential differences in plate
  #efficiences
  
{
  data$tel.cq[i] <- data$ref.tel[i]-data$tel[i]  #Gold tel value - sample tel value
  data$scg.cq[i] <- data$ref.scg[i]-data$scg[i]  #Gold scg value - sample scg value
  data$top[i] <- (data$eff.tel[i])^data$tel.cq[i] #Tel plate efficiency to the power of above
  data$bottom[i] <- (data$eff.scg[i])^data$scg.cq[i] #scg plate efficiency to the power of above
  data$ratio[i] <- ((data$top[i])/(data$bottom[i])) #Relative telomere length
}

#remove triplicates
data_unique <- data %>%
  group_by(Sample) %>% 
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE)))

#write csv file
write.csv(data_unique, file = "RTL_all.csv")

#ICC
data2=read.table ("CT_values_v2.csv", header=T, sep=',')

gold_data <- subset(data2, Sample == "Gold")

gold_tel_data <- subset(gold_data, !is.na(eff.tel))
gold_scg_data <- subset(gold_data, !is.na(eff.scg))

icc_gold_random_batch <- rpt(eff.tel ~ (1|plate), 
                             grname = "plate", 
                             data = gold_tel_data, 
                             datatype = "Gaussian", 
                             nboot = 1000, 
                             npermut = 0)


print(icc_gold_random_batch)

gold_tel_avg <- gold_tel_data %>%
  group_by(plate) %>%
  summarise(eff.tel = mean(eff.tel, na.rm = TRUE))

rpt(eff.tel ~ (1|plate), grname = "plate", data = gold_tel_avg, datatype = "Gaussian", 
    nboot = 1000, npermut=0)

###-----------model selection----------------###

library(glmulti)
library(mice)
library(dplyr)
library(flextable)
library(tidyverse)
library(performance)
library(emmeans)
library(car)
library(effectsize)
library(hier.part)

TL<-read.csv("TL_temporal_transformed.csv")
TL$Year <- as.factor(TL$Year)
TL$Env.type <- as.factor(TL$Env.type)
TL$sex <- as.factor(TL$sex)
TL$tissue_type <- as.factor(TL$tissue_type)
TL$morphotype <- as.factor(TL$morphotype)

#add constant to zscores to shift it into positive, to allow gamma link log function
TL$shifted_zscore <- TL$zscore - min(TL$zscore, na.rm = TRUE) + 1

TL_clean <- TL %>% filter(!is.na(SVL) & !is.na(sex) & !is.na(shifted_zscore))


TLmod <- glm(shifted_zscore ~ Year + Env.type + elevation + RH.mean + wind.mean + 
               sol_rad + morphotype + tissue_type + sex + SVL+TSKY.mean,
             family = Gamma(link = "log"), data = complete_data)


TL_clean$Year <- as.factor(TL_clean$Year)
#g = genetic algorithm
#h = exhaustive run through all models
#d = number of models

#run exhaustive algorithm models
h_model <- glmulti(shifted_zscore ~ Year + Env.type + elevation + RH.mean + wind.mean + 
                     sol_rad + morphotype + tissue_type + sex * SVL + TSKY.mean,
                   data = TL_clean,
                   crit = aicc, 
                   level = 1, 
                   method = 'h',
                   family = Gamma(link ="log"), 
                   fitfunction = glm, 
                   confsetsize = 100)

#best model selected
fmod<-glm(shifted_zscore ~ Year+tissue_type+sex+SVL+TSKY.mean,
          family = Gamma(link = "log"), data = TL_clean)
summary(fmod)

#shows the 100 best models
plot(h_model)

#compare the weights and aicc of top 6 models
weightable(h_model)[1:10,] %>% 
  regulartable() %>% 
  autofit()

#validate the best model by checking which variables reoccurs in the models
plot(h_model, type = "s")

#post-hoc
emm <- emmeans(fmod, ~ Year, type = "response")
pairs(emm, adjust = "tukey")
emm <- emmeans(fmod, ~ tissue_type, type = "response")
pairs(emm, adjust = "tukey")

#effect sizes
std_model <- standardize(fmod)

#extract standardized coefficients
coeff <- summary(std_model)$coefficients

#convert to df
effect_sizes <- data.frame(
  Predictor = rownames(coeff),
  Estimate = coeff[, "Estimate"],
  SE = coeff[, "Std. Error"]
) %>%
  filter(Predictor != "(Intercept)")

#intrinsic vs extrinsic models
#intrinsic variable model
int_mod<-glm(shifted_zscore ~ tissue_type+sex*SVL,
             family = Gamma(link = "log"), data = TL_clean)
#extrinsic variables model
ext_mod<-glm(shifted_zscore ~ Year+TSKY.mean,
             family = Gamma(link = "log"), data = TL_clean)

full_r2<-r2_nagelkerke(fmod)
int_r2<-r2_nagelkerke(int_mod)
ext_r2<-r2_nagelkerke(ext_mod)

# ext/int Deviance explained
dev_explained_model1 <- 1 - (deviance(fmod) / deviance(null_model))
dev_explained_model2 <- 1 - (deviance(int_mod) / deviance(null_model))
dev_explained_model3 <- 1 - (deviance(ext_mod) / deviance(null_model))

#Hierarchical partitioning
hp <- hier.part(TL_clean$shifted_zscore, 
                xcan = predictors, 
                family = Gamma(link = "log"))
                
predictors <- TL_clean[, c("Year", "tissue_type", "sex", "SVL", "TSKY.mean")]
effect_data$predictor <- rownames(effect_data)
effect_data <- hp$I.perc
colnames(effect_data)[1] <- "percent_effect"

#randomisation test of HP
set.seed(123)
rand_test <- rand.hp(TL_clean$shifted_zscore, predictors, family = Gamma(link = "log"), gof = "logLik", num.reps = 100)
print(rand_test)

###-----------------cross-correlation-coefficient analysis with droughts and heatwaves--------------###
library(ggplot2)
library(tidyr)
library(dplyr)

TL<-read.csv("TL_temporal_transformed_hw_dr.csv")

#min-max normalisation function
min_max_norm <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

#normalise considering NAs
TL_normalized <- TL %>%
  mutate(
    heatwave_norm = min_max_norm(Total_Heatwaves),
    drought_norm = min_max_norm(Total_Droughts),
    zscore_norm = min_max_norm(zscore) # Already normalized across all years
  )

#calculate mean
mean_zscore <- TL %>%
  group_by(Year) %>%
  summarise(mean_zscore = mean(zscore_norm, na.rm = TRUE), .groups = "drop")

#extract columns and remove missing data
zscore <- TL$zscore
heatwaves <- TL$Total_Heatwaves
droughts <- TL$Total_Droughts
complete_data <- na.omit(data.frame(zscore, heatwaves, droughts))

#CCF
ccf_zscore_heatwaves <- ccf(complete_data$zscore, complete_data$heatwaves, 
                            lag.max = 5, plot = TRUE,
                            main = "Cross-Correlation: zscore vs Heatwaves")
ccf_zscore_droughts <- ccf(complete_data$zscore, complete_data$droughts, 
                           lag.max = 5, plot = TRUE,
                           main = "Cross-Correlation: zscore vs Droughts")

#extract values and put in df
heatwave_ccf <- as.data.frame(ccf_zscore_heatwaves$acf) %>%
  mutate(Lag = ccf_zscore_heatwaves$lag)
drought_ccf <- as.data.frame(ccf_zscore_droughts$acf) %>%
  mutate(Lag = ccf_zscore_droughts$lag)

ccf_data <- rbind(
  data.frame(Lag = heatwave_ccf$Lag, Correlation = heatwave_ccf$V1, Variable = "Heatwaves"),
  data.frame(Lag = drought_ccf$Lag, Correlation = drought_ccf$V1, Variable = "Droughts")
)

#extract CI
heatwave_ci <- qnorm((1 + 0.95) / 2) / sqrt(nrow(complete_data)) # 95% CI
drought_ci <- qnorm((1 + 0.95) / 2) / sqrt(nrow(complete_data))

#filter for positive lags
ccf_data_positive_lags <- ccf_data %>%
  filter(Lag >= 0)

#Permutation testing
#Function
permutation_test_ccf <- function(zscore, climate_data, lag_max = 5, num_permutations = 5000) {
  observed_ccf <- ccf(zscore, climate_data, lag.max = lag_max, plot = FALSE)$acf
  permuted_ccf_values <- numeric(num_permutations)
  for (i in 1:num_permutations) {
    permuted_zscore <- sample(zscore)
    permuted_ccf <- ccf(permuted_zscore, climate_data, lag.max = lag_max, plot = FALSE)$acf
    permuted_ccf_values[i] <- permuted_ccf[lag_max + 1]
  }
  p_value <- mean(abs(permuted_ccf_values) >= abs(observed_ccf[lag_max + 1]))
  return(list(observed_ccf = observed_ccf[lag_max + 1], p_value = p_value, permuted_ccf_values = permuted_ccf_values))
}

#Perform the tests
heatwave_permutation_results <- permutation_test_ccf(complete_data$zscore, complete_data$heatwaves)
drought_permutation_results <- permutation_test_ccf(complete_data$zscore, complete_data$droughts)

###--------------------Multi-taxa comparison------------------###

library(dplyr)
library(ggplot2)
library(tidyr)
library(rotl)
library(dplyr)
library(tidyr)
library(ape)
library(phytools)

mt<-read.csv("Multi-taxa_zscore_by_species.csv")

#data organisation
mt_long <- pivot_longer(mt, cols = everything(), 
                        names_to = "Dataset", 
                        values_to = "Zscore")

mt_long_f <- mt_long %>%
  filter(!is.na(Zscore)) %>% 
  group_by(Dataset) %>%
  filter(
    Zscore >= quantile(Zscore, 0.25) - 1.5 * IQR(Zscore) & 
      Zscore <= quantile(Zscore, 0.75) + 1.5 * IQR(Zscore)
  ) %>%
  ungroup()

mean_values <- mt_long_f %>%
  group_by(Dataset) %>%
  summarize(
    mean_Zscore = mean(Zscore, na.rm = TRUE),
    sd_Zscore = sd(Zscore, na.rm = TRUE)
  )

#order level
mt_long_f$Dataset <- factor(mt_long_f$Dataset, levels = c(
  "Cocellatus", "Palgirus", "Ggalloti_Gilbert", 
  "Zvivipara", "Lagilis", "Ckingii"
))

mean_values$Dataset <- factor(mean_values$Dataset, levels = c(
  "Cocellatus", "Palgirus", "Ggalloti_Gilbert", 
  "Zvivipara", "Lagilis", "Ckingii"
))

#generate tree
tl<-read.csv("Species_summary.csv")
rownames(tl) <- tl$Species
species_names <- tl$Species

#force the tree to be ultrametric
tree <- compute.brlen(tree, method = "Grafen")

#match data to tree
species_in_tree <- tree$tip.label
species_in_data <- rownames(mean_values)
# Reorder
mean_values <- mean_values[match(tree$tip.label, rownames(mean_values), ]

#pagels lambda
lambda_mean <- phylosig(tree, mean_vec, method = "lambda", test = TRUE)
lambda_SD <- phylosig(tree, sd_vec, method = "lambda", test = TRUE)

#plot trees
phenogram(tree, mean_vec, spread.labels = TRUE, fsize = 0.8, col = "blue", main = "Mean Trait")
phenogram(tree, sd_vec, spread.labels = TRUE, fsize = 0.8, col = "red", main = "SD Trait")