# RNG Prelim Analysis
# Check BGC, DNA data by sites, samples, depths
# by Cliff Bueno de Mesquita, JGI, April 2023

#### 1. Setup ####
library(plyr)
library(tidyverse)
library(car)
library(emmeans)
library(multcomp)
setwd("~/Documents/GitHub/RNG/")
d <- read.csv("data/RNG_metadata.csv") %>%
  filter(DNA_Extraction != "No")



#### 2. Biogeochemistry ####
# Just soils
soil <- d %>%
  filter(Type == "Soil") %>%
  select(-Buffer_pH)
soil_long <- soil %>%
  pivot_longer(cols = c(16:27),
               names_to = "variable") %>%
  mutate(Sample_Site = factor(Sample_Site,
                              levels = c("DS3", "RGO", "HRA", "RGW", "HRC", "RGB"))) %>%
  mutate(variable = factor(variable,
                           levels = c("Clay_perc", "Silt_perc", "Sand_perc",
                                      "OM_perc", "TC_perc", "TN_perc",
                                      "pH", "NO3_mg_kg", "NH4_mg_kg", 
                                      "S_mg_kg", "P_mg_kg", "K_mg_kg")))

bgc_only <- soil %>%
  dplyr::select(16:27)
m <- list()
t <- list()
for (i in 1:ncol(bgc_only)) {
  m[[i]] <- aov(bgc_only[,i] ~ soil$Sample_Site)
  t[[i]] <- emmeans(object = m[[i]], specs = "Sample_Site") %>%
    cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
    mutate(name = names(bgc_only)[i],
           y = max(bgc_only[,i], na.rm = T)+
             (max(bgc_only[,i], na.rm = T)-min(bgc_only[,i], na.rm = T))/10) %>%
    rename(variable = name)
}
t_df <- do.call(rbind.data.frame, t) %>%
  mutate(variable = factor(variable, levels = levels(soil_long$variable)))

facet_names <- c("Clay_perc" = "Clay~('%')",
                 "Silt_perc" = "Silt~('%')",
                 "Sand_perc" = "Sand~('%')",
                 "OM_perc" = "OM~('%')",
                 "TC_perc" = "TC~('%')",
                 "TN_perc" = "TN~('%')",
                 "pH" = "pH",
                 "NO3_mg_kg" = "NO[3]^'-'~(mg/kg)",
                 "NH4_mg_kg" = "NH[4]^'+'~(mg/kg)",
                 "S_mg_kg" = "S~(mg/kg)",
                 "P_mg_kg" = "P~(mg/kg)",
                 "K_mg_kg" = "K~(mg/kg)")

png("InitialFigs/BGC.png", width = 7, height = 7, units = "in", res = 300)
ggplot(soil_long, aes(x = Sample_Site,
              y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75) +
  geom_text(data = t_df, aes(Sample_Site, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Site",
       y = NULL) +
  facet_wrap(~ variable, ncol = 3, scales = "free_y",
             labeller = as_labeller(facet_names, default = label_parsed)) +
  theme_bw()
dev.off()


#### 3. DNA Concentration ####
m <- aov(DNA_conc_ng_uL ~ Sample_Site + Depth, data = d)
hist(m$residuals)
plot(m$fitted.values, m$residuals)
Anova(m, type = "II")

ggplot(d, aes(Sample_Site, DNA_conc_ng_uL)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  labs(x = "Site",
       y = "DNA concentration (ng/µL)") +
  facet_grid(~ Type, scales = "free_x", space = "free") +
  theme_bw()

png("InitialFigs/DNA_conc.png", width = 6, height = 4, units = "in", res = 300)
ggplot(d, aes(x = reorder(Sample_Site, DNA_conc_ng_uL, mean),
              y = DNA_conc_ng_uL, 
              colour = Type)) +
  geom_hline(yintercept = 10, linetype = "dotted") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.75, position = position_jitterdodge()) +
  scale_colour_manual(values = c("grey50", "tan", "brown", "skyblue")) +
  labs(x = "Site",
       y = "DNA concentration (ng/µL)") +
  theme_bw() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank())
dev.off()
