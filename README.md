# Schisto burden paper
---
title: "Schisto burden"
author: "Sidy"
date: "2025-05"
output: html_document
---
#####This analysis is from baseline parasitological. I uploaded firt the individual level data, and the site level.
```{r}
####Individual prevalence DAta loarding
Parasitobase= read.csv2("~/Desktop/Developpement/DISES/DISES-analysis/Parasitology_Data/Parasitobase.csv")
```
### 1. Needed packages
```{r}
###Needed packages
library(binom)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(lme4)
library(tidyverse)
library(MASS)
library(car)
library(dunn.test)
library(broom)###extract t-test results and present them in a clean table
library(broom.mixed)
```

### 2. Data organization
```{r}
##### Data organization

# Convert ALL character columns that look numeric
IPB <- Parasitobase####Copy of the data to clean and use for analysis
# Columns conversion that should be numeric
num_cols <- c("S.mansoni.EPG","S.mansoni.mean.EPG", "Total.living.S.h.eggs", 
               "S.haematobium.infection..0.1.", "S.mansoni.infection..0.1.",
              "Latitude", "Longitude", "Actual.age","S.haematobium.living.eggs.round")

IPB[num_cols] <- lapply(IPB[num_cols], function(x) as.numeric(as.character(x)))
colSums(is.na(IPB))
IPB<- IPB[!is.na(IPB$S.mansoni.EPG) & !is.na(IPB$Total.living.S.h.eggs), ]
IPB <- IPB[!duplicated(IPB), ]
IPB$Sexe <- as.factor(IPB$Sexe)
IPB$SRB <- as.factor(IPB$SRB)
IPB$Hydrology.source <- as.factor(IPB$Hydrology.source)
IPB$scale_Lat <- scale(IPB$Latitude)####For GLMM model
IPB$scale_Long <- scale(IPB$Longitude)
IPB$scale_Age <- scale(IPB$Actual.age)

IPB <- IPB %>%
  dplyr::filter(
    !is.na(S.mansoni.infection..0.1.),
    !is.na(Hydrology.source),
    !is.na(Longitude),
    !is.na(Latitude),
    !is.na(Sexe),
    !is.na(Actual.age),
    !is.na(S.haematobium.infection..0.1.),
    !is.na(S.haematobium.living.eggs.round),
    !is.na(S.mansoni.mean.EPG))

IPB$S.haematobium.eggs.round[IPB$S.haematobium.eggs.round=="-"]<-0
IPB$S.haematobium.infection.intensity[IPB$S.haematobium.infection.intensity=="-"]<-0
IPB$S.haematobium.infection..0.1.[IPB$S.haematobium.infection..0.1.=="-"]<-0
IPB$Total.living.S.h.eggs[IPB$Total.living.S.h.eggs=="-"]<-0
```

```{r}
###### Data Summury
summary(IPB)
```

### 3. Prevalence comparison
```{r}
#####T-test: comparison male vs female infection prevalence 
t.test(S.mansoni.infection..0.1.~Sexe, data=IPB) ####mansoni
t.test(S.haematobium.infection..0.1.~Sexe,data = IPB)####Haematobium
#####Comparison males vs females Infection intensity  
t.test(S.mansoni.mean.EPG~Sexe, data=IPB)####mansoni
t.test(Total.living.S.h.eggs~Sexe, data=IPB)######haematobium

#####Comparison Delta vs Valley infection prevalence
t.test(S.mansoni.infection..0.1.~SRB, data = IPB)#####mansoni
t.test(S.haematobium.infection..0.1.~SRB, data = IPB)#####haematobium
######Comparison Delta vs Valley infection intensity
t.test(S.mansoni.mean.EPG~SRB, data=IPB)###mansoni
t.test(Total.living.S.h.eggs~SRB, data=IPB) ####Haematobium


```


```{r}
######Global prevalence
#####Global mansoni prevalence
S_mansoni_prevalence <- IPB %>%
  summarise(
    N = n(),
    Positive = sum(S.mansoni.infection..0.1. == 1, na.rm = TRUE),
    Prevalence = round(100 * Positive / N, 1))
S_mansoni_prevalence
S_mansoni_prevalence_CI <-S_mansoni_prevalence %>%
  rowwise() %>%
  mutate(
    CI = list(binom.confint(Positive, N, method = "wilson")),
    CI_lower = CI$lower * 100,
    CI_upper = CI$upper * 100)
S_mansoni_prevalence_CI

##### Global haematobium prevalence
S_haematobium_prevalence <- IPB %>%
  summarise(
    N = n(),
    Positive = sum(S.haematobium.infection..0.1. == 1, na.rm = TRUE),
    Prevalence = round(100 * Positive / N, 1))
S_haematobium_prevalence
S_haematobium_prevalence_CI <-S_haematobium_prevalence %>%
  rowwise() %>%
  mutate(
    CI = list(binom.confint(Positive, N, method = "wilson")),
    CI_lower = CI$lower * 100,
    CI_upper = CI$upper * 100)
S_haematobium_prevalence_CI

########Regional prevalence
####S h region
S_mansoni_prevalence_by_region <- IPB %>%
  group_by(SRB) %>%
  summarise(
    N = n(),
    Positive = sum(S.mansoni.infection..0.1. == 1, na.rm = TRUE),
    Prevalence = round(100 * Positive / N, 1))
S_mansoni_prevalence_by_region
####Sm P CI by region
S_mansoni_prevalence_by_region_CI <-S_mansoni_prevalence_by_region %>%
  rowwise() %>%
  mutate(
    CI = list(binom.confint(Positive, N, method = "wilson")),
    CI_lower = CI$lower * 100,
    CI_upper = CI$upper * 100)
S_mansoni_prevalence_by_region_CI
#####S h region
S_haematobium_prevalence_by_region <- IPB %>%
  group_by(SRB) %>%
  summarise(
    N = n(),
    Positive = sum(S.haematobium.infection..0.1. == 1, na.rm = TRUE),
    Prevalence = round(100 * Positive / N, 1))
S_haematobium_prevalence_by_region
####Sh P CI by region
S_haematobium_prevalence_by_region_CI <-S_haematobium_prevalence_by_region %>%
  rowwise() %>%
  mutate(
    CI = list(binom.confint(Positive, N, method = "wilson")),
    CI_lower = CI$lower * 100,
    CI_upper = CI$upper * 100)
S_haematobium_prevalence_by_region_CI
######Sites level prevaelnces
######Sites mansoni prevalence
S_mansoni_prevalence_by_site <- IPB %>%
  group_by(Sites) %>%
  summarise(
    N = n(),
    Positive = sum(S.mansoni.infection..0.1. == 1, na.rm = TRUE),
    Prevalence = round(100 * Positive / N, 1))
S_mansoni_prevalence_by_site
######S mansoni prevalence with Confidence Interval (CI)
S_mansoni_prevalence_by_site_CI <-S_mansoni_prevalence_by_site %>%
  rowwise() %>%
  mutate(
    CI = list(binom.confint(Positive, N, method = "wilson")),
    CI_lower = CI$lower * 100,
    CI_upper = CI$upper * 100)
S_mansoni_prevalence_by_site_CI

#####Sites haematobium prevalence
S_haematobium_prevalence_by_site <- IPB %>%
  group_by(Sites) %>%
  summarise(
    N = n(),
    Positive = sum(S.haematobium.infection..0.1. == 1, na.rm = TRUE),
    Prevalence = round(100 * Positive / N, 1))
S_haematobium_prevalence_by_site
###### S haematobium preva with  CI
S_haematobium_prevalence_by_site_CI <-S_haematobium_prevalence_by_site %>%
  rowwise() %>%
  mutate(
    CI = list(binom.confint(Positive, N, method = "wilson")),
    CI_lower = CI$lower * 100,
    CI_upper = CI$upper * 100)
S_haematobium_prevalence_by_site_CI

```
##### 4. Sites Mean eggs
```{r}
###### Total mean eggs by specy
####Mansoni
S_mansoni_EPG <- IPB %>%
  summarise(
    N = n(),
    Mean_mansoni_EPG = mean(S.mansoni.EPG, na.rm = TRUE),
    SD_EPG = sd(S.mansoni.EPG, na.rm = TRUE))
S_mansoni_EPG
####Haematobium
S_haematobium_eggs <- IPB %>%
  summarise(
    N = n(),
    Mean_haematobium_eggs = mean(S.haematobium.living.eggs.round, na.rm = TRUE),
    SD_haematobium_eggs = sd(S.haematobium.living.eggs.round, na.rm = TRUE))
S_haematobium_eggs
#####S mansoni with standard deviation by site
S_mansoni_eggs<-IPB %>%
  group_by(Sites) %>%
  summarise(
    N = n(),
    Mean_mansoni_EPG = mean(S.mansoni.EPG, na.rm = TRUE),
    SD_EPG = sd(S.mansoni.EPG, na.rm = TRUE))
S_mansoni_eggs
######S. haematobium mean eggs and standard dev

S_haematobium_eggs<-IPB %>%
  group_by(Sites) %>%
  summarise(
    N = n(),
    Mean_haematobium_eggs = mean(S.haematobium.living.eggs.round, na.rm = TRUE),
    SD_haem_Eggs = sd(S.haematobium.living.eggs.round, na.rm = TRUE))
S_haematobium_eggs

```

###### 5. Coinfection status

```{r}
######Coinfection status
S_sp_coinfection <- with(IPB, ifelse(S.mansoni.infection..0.1. == 1 & S.haematobium.infection..0.1. == 1, "Both",
                                  ifelse(S.mansoni.infection..0.1. == 1, "S_mansoni_only",
                                  ifelse(S.haematobium.infection..0.1. == 1, "S_haematobium_only", "None"))))
#####prevalence calcul
table(S_sp_coinfection)
prop.table(table(S_sp_coinfection)) * 100  # in %


```
#### 6. Schistosoma intensity Classifification 
```{r}
# S haematobium Classifification intensity
Intensityhae <- cut(IPB$Total.living.S.h.eggs,breaks = c(-1, 0, 49, Inf),
                     labels = c("Uninfected", "Light", "Heavy"))
table(Intensityhae)### %
prop.table(table(Intensityhae)) * 100  # percentages
IPB <- IPB %>% filter(!is.na(Intensityhae))

# S mansoni Classification intensity 
Intensityman <- cut(IPB$S.mansoni.EPG,
                        breaks =c(-1, 0, 99, 399, Inf),
                        labels = c("Uninfected", "Light", "Moderate", "Heavy"))

table(Intensityman) #### %
prop.table(table(Intensityman)) * 100  # percentages
IPB <- IPB %>% filter(!is.na(Intensityman))
#################PLot


Im=ggplot(IPB, aes(x =Intensityman, fill = Intensityman)) +
  geom_bar() +
  scale_fill_manual(values = c("gray", "lightblue","orange3","tomato")) +
  labs(title = "C",
       x = expression(italic("S. mansoni")),
       y = "Number of tested children",
       fill="Intensity") +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank())+ 
  theme(text = element_text(family = "Times New Roman"),
        axis.line = element_line(color = "black", linewidth = 0.5),   # draw axis lines
        axis.ticks = element_line(color = "black"),                   # draw tick marks
        axis.ticks.length = unit(0.2, "cm")  )


Ih=ggplot(IPB, aes(x = Intensityhae, fill = Intensityhae)) +
  geom_bar() +
  scale_fill_manual(values = c("gray", "lightblue", "tomato")) +
  labs(title = "D",
       x = expression(italic("S. haematobium")),
       y = "Number of tested children",
       fill="Intensity") +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")+
  theme(text = element_text(family = "Times New Roman"),
        axis.line = element_line(color = "black", linewidth = 0.5),   
        axis.ticks = element_line(color = "black"),                   
        axis.ticks.length = unit(0.2, "cm"))
Im+Ih



```
######### 7. Prevalence, Intensity Plotting 
```{r}
####Site data level loading
Baseline.Prevalence <- read.csv("~/Desktop/Developpement/DISES/DISES-analysis/Parasitology_Data/Baseline Prevalence.csv", sep=";")
attach(Baseline.Prevalence)
Baseline.Prevalence$Latitude <- as.numeric(gsub(",", ".", Baseline.Prevalence$Latitude))
Baseline.Prevalence$Longitude <- as.numeric(gsub(",", ".", Baseline.Prevalence$Longitude))


Ph=ggplot(data = Baseline.Prevalence, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(size = S..h_P, color = S..h_P), alpha = 0.9) +
  guides(size = guide_legend(title = "S.h P (%)"),
         color = guide_colorbar(title = "S.h P (%)"))+
  scale_color_viridis_c(option = "E",name = "Eggs per 10_ml") +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  labs(title = "S.h",
       x = "Longitude", y = "Latitude",
       fill="Intensity")+
  theme(legend.position = "none")+
  theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),   # draw axis lines
        axis.ticks = element_line(color = "black"),                   # draw tick marks
        axis.ticks.length = unit(0.2, "cm"))
Ph
Pm=ggplot(data = Baseline.Prevalence, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(size = S..m_P, color = S..m_P), alpha = 0.9) +
  guides(size = guide_legend(title = "S. sp P (%)"),
         color = guide_colorbar(title = "S.sp P (%)"))+
  scale_color_viridis_c(option = "E", name = "EPG") +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  labs(title = "S.m",
       x = "Longitude", y = "Latitude",
       fill="Intensity")+
  theme(legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),  
        axis.ticks = element_line(color = "black"),                   
        axis.ticks.length = unit(0.2, "cm"))

Pm

Ih=ggplot(data = Baseline.Prevalence, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(size = S.h_I, color = S.h_I), alpha = 0.9) +
  guides(size = guide_legend(title = "S.h P (%)"),
         color = guide_colorbar(title = "S.h P (%)"))+
  scale_color_viridis_c(option = "C", name = "Eggs per 10_ml") +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  labs(title = "S.h",
       x = "Longitude", y = "Latitude",
       fill="Intensity")+
  theme(legend.key.size = unit(0.3, "cm"))+
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),   # draw axis lines
    axis.ticks = element_line(color = "black"),                   # draw tick marks
    axis.ticks.length = unit(0.2, "cm"))
Ih
Im=ggplot(data = Baseline.Prevalence, aes(x = Longitude, y = Latitude)) +
  geom_point(aes(size = S.m_I, color = S.m_I), alpha = 0.9) +
  guides(size = guide_legend(title = "S. sp I (Eggs)"),
         color = guide_colorbar(title = "S.sp I (Eggs)"))+
  scale_color_viridis_c(option = "C", name = "EPG") +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  labs(title = "S.m",
       x = "Longitude", y = "Latitude",
       fill="Intensity")+
  theme(legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        text = element_text(family = "Times New Roman"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),  
        axis.ticks = element_line(color = "black"),                   
        axis.ticks.length = unit(0.2, "cm"))
IM

Ph+Pm+Ih+Im


```



### 8. S. sp modelling
##### 8.1.S mansoni infection modelling
```{r}
#######Statistical test
IPB$Actual.age=as.numeric(IPB$Actual.age)
chisq.test(table(SRB,S.mansoni.infection..0.1.))
fisher.test(table(SRB,S.mansoni.infection..0.1.))
####S mansoniLogistical model
modelmansoni<-glm(S.mansoni.infection..0.1. ~Hydrology.source+Longitude+Latitude+ Sexe+Actual.age, data=IPB, family = "binomial")
summary(modelmansoni)

modelmansonilong<-glm(S.mansoni.infection..0.1. ~Hydrology.source*Longitude+Latitude+ Sexe+Actual.age, data=IPB, family = "binomial")###Interaction hydro source*long
summary(modelmansonilong)

modelmansonilat<-glm(S.mansoni.infection..0.1. ~Hydrology.source*Latitude+Longitude+ Sexe+Actual.age, data=IPB, family = "binomial")###Interaction source and lat
summary(modelmansonilat)

modelmansonisex<-glm(S.mansoni.infection..0.1. ~Hydrology.source*Sexe+Longitude+Latitude+ Actual.age, data=IPB, family = "binomial")###Interaction source and sexe
summary(modelmansonisex)

modelmansoniage<-glm(S.mansoni.infection..0.1. ~Hydrology.source* Actual.age+Sexe+Longitude+Latitude, data=IPB, family = "binomial")###Interaction source and age
summary(modelmansoniage)
AIC(modelmansoni,modelmansonilong,modelmansonilat,modelmansonisex,modelmansoniage)

vif(modelmansonilat)####Smallest AIC but Strong colinearity detected in interaction term Hydro source*Lat.NOT BEST MODEL
vif(modelmansonilong)##### 2nd Not best model. Severe multicolinearity
vif(modelmansoniage) ##### 3rd smallest AIC. Colinearity Drop all model interaction

vif(modelmansoni)###4 smallest AIC ####Best model 
```
###### 8.2. Best mansoni infection model
```{r}
modelmansoni<-glm(S.mansoni.infection..0.1. ~Hydrology.source+Longitude+Latitude+ Sexe+Actual.age, data=IPB, family = "binomial")
summary(modelmansoni)

coefmodelmansoni=coef(modelmansoni)
ORmansoni<-exp(coefmodelmansoni)#####)Odd Ratio (OR)
ORmansoni
ORmansoni_decimal <- format(ORmansoni, scientific = FALSE)
ORmansoni_decimal

# Get summary to extract p-values
summary_model <- summary(modelmansoni)

# Combine ORs, CIs, and p-values
OR_table <- cbind(
  OR = exp(coef(modelmansoni)),
  Lower_CI = exp(confint(modelmansoni)[, 1]),
  Upper_CI = exp(confint(modelmansoni)[, 2]),
  p_value = summary_model$coefficients[, 4])

# View rounded results
round(OR_table, 3)

```

##### 8.3. S haematobium infection modelling

```{r}
####S haematobium infection logistical model
modelhaematobium<-glm(S.haematobium.infection..0.1. ~SRB+Hydrology.source+Latitude+Longitude+Sexe+Actual.age, data=IPB, family = "binomial")
summary(modelhaematobium)

modelhaematobiumlong<-glm(S.haematobium.infection..0.1.  ~SRB+Hydrology.source*Longitude+Latitude+ Sexe+Actual.age, data=IPB, family = "binomial")###Interaction hydro source*long
summary(modelhaematobiumlong)


modelhaematobiumlat<-glm(S.haematobium.infection..0.1.  ~SRB+Hydrology.source*Latitude+Longitude+Sexe+Actual.age, data=IPB, family = "binomial")###Interaction hydro source*lat
summary(modelhaematobiumlat)


modelhaematobiumsex<-glm(S.haematobium.infection..0.1.  ~SRB+Hydrology.source*Sexe+Longitude+Latitude+ Actual.age, data=IPB, family = "binomial")###Interaction source*sexe
summary(modelhaematobiumsex)


modelhaematobiumage<-glm(S.haematobium.infection..0.1.  ~SRB+Hydrology.source* Actual.age+Sexe+Longitude+Latitude, data=IPB, family = "binomial")###Interaction source vs age
summary(modelhaematobiumage)

AIC(modelhaematobium,modelhaematobiumlong,modelhaematobiumlat,modelhaematobiumsex,modelhaematobiumage)

vif(modelhaematobiumlat)####Smallest AIC but Strong colinearity detected in interaction term Hydro source*Lat.NOT BEST MODEL

vif(modelhaematobium)####Best model 2nd smallest AIC model, but with no interaction.
```

##### 8.4. Haematobium infection best model
```{r}
modelhaematobium<-glm(S.haematobium.infection..0.1. ~SRB+Hydrology.source+Latitude+Longitude+Sexe+Actual.age, data=IPB, family = "binomial")
summary(modelhaematobium)

coefmodelhaematobium=coef(modelhaematobium)
ORhaem<-exp(coefmodelhaematobium)#####)Odd Ratio (OR)
ORhaem
ORhaem_decimal <- format(ORhaem, scientific = FALSE)
ORhaem_decimal
# 1. Get confidence intervals on the log-odds scale
confint_log <- confint(modelhaematobium)

# 2. Exponentiate to get confidence intervals on the odds ratio scale
confint_odds <- exp(confint_log)

# 3. Get the odds ratios
odds_ratios <- exp(coef(modelhaematobium))
odds_ratios_dec<-format(odds_ratios,scientific=FALSE)

# Combine ORs, CIs, and p-values
OR_t <- cbind(
  OR = exp(coef(modelhaematobium)),
  Lower_CI = exp(confint(modelhaematobium)[, 1]),
  Upper_CI = exp(confint(modelhaematobium)[, 2]),
  p_value = summary_model$coefficients[, 4])
 
print(OR_t)
```

###### 9. Infection intensity modeling
```{r}
######Intensity infection modeling
library(glmmTMB)
####Mansoni intensity
Mansoni_intensity <- glmmTMB(
  S.mansoni.EPG ~ Hydrology.source + scale(Longitude) + scale(Latitude) + Sexe + scale(Actual.age) + (1 | Sites),
  family = nbinom2,
  data = IPB)
Mansoni_intensity
######INteraction between hydrology and other factors
Mansoni_intensitylong <- glmmTMB(
  S.mansoni.EPG ~ Hydrology.source*scale(Longitude) + scale(Latitude) + Sexe + scale(Actual.age) + (1 | Sites),
  family = nbinom2,
  data = IPB)
Mansoni_intensitylong

Mansoni_intensitylat <- glmmTMB(
  S.mansoni.EPG ~ Hydrology.source*scale(Latitude)+scale(Longitude) + Sexe + scale(Actual.age) + (1 | Sites),
  family = nbinom2,
  data = IPB)
Mansoni_intensitylat

Mansoni_intensitysex <- glmmTMB(
  S.mansoni.EPG ~ Hydrology.source*Sexe +scale(Latitude)+scale(Longitude)+ scale(Actual.age) + (1 | Sites),
  family = nbinom2,
  data = IPB)
Mansoni_intensitysex#####best model according to AIC

Mansoni_intensityage <- glmmTMB(
  S.mansoni.EPG ~ Hydrology.source*scale(Actual.age)+Sexe +scale(Latitude)+scale(Longitude) +(1 | Sites),
  family = nbinom2,
  data = IPB)
Mansoni_intensityage

AIC(Mansoni_intensity,Mansoni_intensitylong,Mansoni_intensitylat,Mansoni_intensitysex,Mansoni_intensityage)

exp(fixef(Mansoni_intensity_SRB)$cond)####Ratio
confint(Mansoni_intensity_SRB, level = 0.95)###CI
confint(Mansoni_intensity, method = "wald")

########Haematobium intensity
Haematobium_intensity <- glmmTMB(
  S.haematobium.living.eggs.round ~ SRB + Hydrology.source + scale(Longitude) + scale(Latitude) + Sexe + scale(Actual.age) + (1 | Sites),
  family = nbinom2,
  data = IPB)
Haematobium_intensity

#####Interaction between hydrology source and other factors
Haematobium_intensitylong <- glmmTMB(
  S.haematobium.living.eggs.round ~ SRB + Hydrology.source*scale(Longitude) + scale(Latitude) + Sexe + scale(Actual.age) + (1 | Sites),
  family = nbinom2,
  data = IPB)
Haematobium_intensitylong

Haematobium_intensitylat <- glmmTMB(
  S.haematobium.living.eggs.round ~ SRB + Hydrology.source*scale(Latitude)+scale(Longitude) + Sexe + scale(Actual.age) + (1 | Sites),
  family = nbinom2,
  data = IPB)
Haematobium_intensitylat

Haematobium_intensitysex <- glmmTMB(
  S.haematobium.living.eggs.round ~ SRB + Hydrology.source*Sexe +scale(Latitude)+scale(Longitude)+ scale(Actual.age) + (1 | Sites),
  family = nbinom2,
  data = IPB)
Haematobium_intensitysex#####best model according to AIC

Haematobium_intensityage <- glmmTMB(
  S.haematobium.living.eggs.round ~ SRB + Hydrology.source*scale(Actual.age)+Sexe +scale(Latitude)+scale(Longitude) +(1 | Sites),
  family = nbinom2,
  data = IPB)
Haematobium_intensityage


AIC(Haematobium_intensity,Haematobium_intensitylong,Haematobium_intensitylat,Haematobium_intensitysex,Haematobium_intensityage)
```

##### 10. Best intensity infection model extraction
```{r}
########Haematobium
Haematobium_intensitysex <- glmmTMB(
  S.haematobium.living.eggs.round ~ SRB + Hydrology.source*Sexe +scale(Latitude)+scale(Longitude)+ scale(Actual.age) + (1 | Sites),
  family = nbinom2,
  data = IPB)
Haematobium_intensitysex#####best model according to AIC

tidy(Haematobium_intensitysex, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)

######Mansoni
Mansoni_intensitysex <- glmmTMB(
  S.mansoni.EPG ~ Hydrology.source*Sexe +scale(Latitude)+scale(Longitude)+ scale(Actual.age) + (1 | Sites),
  family = nbinom2,
  data = IPB)
Mansoni_intensitysex#####best model according to AIC

tidy(Mansoni_intensitysex, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)



######Coef calculating 

#######Haematobium intensity coef
confint(Haematobium_intensitysex, method = "profile")  # Better with small samples
coefsHaematobium_intensitysex <- summary(Haematobium_intensitysex)$coefficients$cond
coefsHaematobium_intensitysex
# SE, z, p, CI Estimations
resultsHaematobium_intensitysex <- broom.mixed::tidy(Haematobium_intensitysex, effects = "fixed", conf.int = TRUE)
resultsHaematobium_intensitysex
resultsHaematobium_intensitysex$rate_ratio <- exp(resultsHaematobium_intensitysex$estimate)
resultsHaematobium_intensitysex$CI_low <- exp(resultsHaematobium_intensitysex$conf.low)
resultsHaematobium_intensitysex$CI_high <- exp(resultsHaematobium_intensitysex$conf.high)
resultsHaematobium_intensitysex

######Mansoni infection coef
confint(Mansoni_intensitysex, method = "profile")  # Better with small samples
coefsMansoni_intensitysex <- summary(Mansoni_intensitysex)$coefficients$cond
coefsMansoni_intensitysex
# , SE, z, p, CI Estimation
resultsMansoni_intensitysex <- broom.mixed::tidy(Mansoni_intensitysex, effects = "fixed", conf.int = TRUE)
resultsMansoni_intensitysex
resultsMansoni_intensitysex$rate_ratio <- exp(resultsMansoni_intensitysex$estimate)
resultsMansoni_intensitysex$CI_low <- exp(resultsMansoni_intensitysex$conf.low)
resultsMansoni_intensitysex$CI_high <- exp(resultsMansoni_intensitysex$conf.high)
resultsMansoni_intensitysex

```




