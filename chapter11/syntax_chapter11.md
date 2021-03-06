---
output:
  html_document: default
  pdf_document: default
---
# Syntax Chapter 11  

**This chapter needs to be exectuted in R**

## How to get the HRS data
```
install.packages("devtools")
library(devtools)
install_github("ajdamico/lodown")
library(lodown)
```

```
# If you would like to download all of the HRS data that you see on the website, this would be the command:
lodown( "hrs" , output_dir = getwd() , 
        your_username = "your username" , 
        your_password = "your password" )
```
 
```
# For this exercise, we do not need all of the data available from HRS, but can focus on specific versions and subsets, which we can download like this:
# Create a list of available data sets and define the download # destination using output_dir:
hrs_cat <- get_catalog( "hrs" ,
               output_dir = getwd() , 
               your_username = "your username" , 
               your_password = "your password" )
               
# If you execute the command, in RStudio, you can see the available datasets in parallel to the website:
view(hrs_cat)

# Define the first dataset we wish to download, namely the PGSs: PGENSCORE3.zip:
PGS <- subset( hrs_cat , grepl( 'PGENSCORE3.zip' , file_name ) )

# Download the data:
PGS <- lodown( "hrs" , PGS,
               your_username = "your username" , 
               your_password = "your password" )
 
# Next define the second dataset we wish to download, namely the RAND longitudinal data: randhrs1992_2014v2_STATA.zip:
Pheno <- subset( hrs_cat , grepl( 'randhrs1992_2014v2_STATA.zip' , file_name ))

# Download the data:
Pheno <- lodown( "hrs" , Pheno,
               your_username = "your username" , 
               your_password = "your password" )              
```

```
# Read in the data frame for the PGS:
hrs_PGS <- readRDS("PGENSCORE3/pgenscore3e_r.rds")

# Read in the data frame for the Phenos:
hrs_Phenos <- readRDS("randhrs1992_2014v2_STATA/randhrs1992_2014v2.rds")

# Merge data frames
hrs_GePhen <- merge(hrs_PGS,hrs_Phenos,by=c("hhid","pn"))
```

```
install.packages("dplyr")
library(dplyr)

# Select one individual per household
hrs_GePhen_ci <-hrs_GePhen %>% group_by(hhid) %>% dplyr::mutate(whh_count = row_number())

# We keep arbitrarily always the first person in the household
hrs_GePhen_uni <- hrs_GePhen_ci[ which(hrs_GePhen_ci$wf_count =="1"),]

# Finally, we clean the working space in RStudio
rm(hrs_PGS,hrs_GePhen_ci,hrs_GePhen,hrs_Phenos)

# The data file we are using from here is named hrs_GePhen_uni and contains 8,451 individuals and 11,523 variables. Be sure before embarking on the analyses in Chapter 11 that you also have the same number. 
```

## 11.2 Polygenic Score (PGS) Applications
```
# Select variables
vars_prediction <- c("raedyrs","ea_pgs3_edu2_ssgac16", 
                     "ea_pgs3_edu3_ssgac18","rabyear","ragender","pc1_5a", 
                     "pc1_5b","pc1_5c","pc1_5d","pc1_5e","pc6_10a", "pc6_10b", 
                     "pc6_10c","pc6_10d","pc6_10e")

# Summarize variables
data_prediction <- hrs_GePhen_uni[vars_prediction]
summary(data_prediction$raedyrs)

# Select and summarise score
d_PGS3 <- density(data_prediction$ea_pgs3_edu3_ssgac18)
plot(d_PGS3, main= "Kernel Density PGS3")
```

```
# Estimate models
# Model 1: Edu = PGS2 + Cov
predic_control_PGS2 <- lm(raedyrs~ ea_pgs3_edu2_ssgac16+rabyear+
                            ragender +pc1_5a+pc1_5b+pc1_5c+pc1_5d +pc1_5e+ +
                            pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e, data = data_prediction)

# Model 2: Edu = PGS3 + Cov
predic_control_PGS3 <- lm(raedyrs~ ea_pgs3_edu3_ssgac18+rabyear+ragender 
                          +pc1_5a+pc1_5b+pc1_5c+pc1_5d +pc1_5e+pc6_10a+pc6_10b+pc6_10c+
                            pc6_10d+pc6_10e, data = data_prediction)

# Model 2: Edu = Cov
predic_control <- lm(raedyrs~ rabyear+ragender +pc1_5a+ pc1_5b+pc1_5c+
                       pc1_5d+pc1_5e+pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data=data_prediction)

summary(predic_control_PGS2)
summary(predic_control_PGS2)$r.square-summary(predic_control)$r.square
summary(predic_control_PGS3)$r.square-summary(predic_control)$r.square

# Clean R
rm(vars_prediction,predic_control,predic_control_PGS2,
   predic_control_PGS3,data_prediction,d_EduYears, d_PGS2,d_PGS3)
```

```
# Average BMI across waves: BMI_AV
bmi <-c("r1bmi", "r2bmi", "r3bmi", "r4bmi", "r5bmi", 
        "r6bmi", "r7bmi", "r8bmi", "r9bmi", "r10bmi", 
        "r11bmi", "r12bmi")

hrs_GePhen_uni$BMI_AV = rowMeans(hrs_GePhen_uni[,bmi], na.rm = TRUE) 
summary(hrs_GePhen_uni$BMI_AV)

# Average height across waves: Height_AV
height <-c("r1height", "r2height", "r3height", "r4height", 
           "r5height", "r6height", "r7height", "r8height", 
           "r9height", "r10height", "r11height", "r12height")
           
hrs_GePhen_uni$Height_AV = rowMeans(hrs_GePhen_uni[,height], na.rm = TRUE)

# Average Age across waves: Age_AV
age <-c("r1agey_b", "r2agey_b","r3agey_b","r4agey_b","r5agey_b",
        "r6agey_b","r7agey_b","r8agey_b","r9agey_b","r10agey_b",
        "r11agey_b", "r12agey_b")
        
hrs_GePhen_uni$Age_AV = rowMeans(hrs_GePhen_uni[,age], na.rm = TRUE)

# Average Health across waves: Health_AV
health <- c("r1conde", "r2conde","r3conde", "r4conde", "r5conde", 
            "r6conde", "r7conde", "r8conde", "r9conde", "r10conde", 
            "r11conde", "r12conde")
            
hrs_GePhen_uni$Health_AV = rowMeans(hrs_GePhen_uni[,health], na.rm = TRUE)
summary(hrs_GePhen_uni$Health_AV)
```

```
# Select variables for analyses
vars_rG <- c("Height_AV", "BMI_AV", "Health_AV", "ea_pgs3_edu3_ssgac18", 
             "Age_AV", "raedyrs", "raevbrn", "rabyear", "ragender", 
             "pc1_5a", "pc1_5b","pc1_5c","pc1_5d", "pc1_5e", 
             "pc6_10a", "pc6_10b", "pc6_10c", "pc6_10d", "pc6_10e")

data_rG <- hrs_GePhen_uni[vars_rG]
```

```
# Center variables
data_rG$Height_AV_scaled =scale(data_rG$Height_AV) 
data_rG$BMI_AV_scaled = scale(data_rG$BMI_AV)
data_rG$raedyrs_scaled =scale(data_rG$raedyrs) 
data_rG$Health_AV_scaled =scale(data_rG$Health_AV) 
data_rG$raevbrn_scaled =scale(data_rG$raevbrn)
```

```
# Estimate models
M_rGHeight <- lm(Height_AV_scaled~ ea_pgs3_edu3_ssgac18+rabyear+ragender 
                 +pc1_5a+pc1_5b+pc1_5c+pc1_5d +pc1_5e+
                   pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data =data_rG)

M_rHeight <- lm(Height_AV_scaled~ raedyrs_scaled +rabyear+ragender +
                  pc1_5a+pc1_5b+pc1_5c+pc1_5d+pc1_5e+
                  pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data=data_rG)

M_rGBMI <- lm(BMI_AV_scaled~ ea_pgs3_edu3_ssgac18+rabyear+ragender 
              +pc1_5a+pc1_5b+pc1_5c+pc1_5d +pc1_5e+
                pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data =data_rG)
 
M_rBMI<- lm(BMI_AV_scaled~ raedyrs_scaled+rabyear+ragender +
              pc1_5a+pc1_5b+pc1_5c+pc1_5d+ pc1_5e+
              pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data= data_rG)
   
M_rGFert <- lm(raevbrn_scaled~ ea_pgs3_edu3_ssgac18+rabyear+ragender +
                 pc1_5a+pc1_5b+pc1_5c+pc1_5d+pc1_5e+
                 pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data =data_rG)
 
M_rFert<- lm(raevbrn_scaled~ raedyrs_scaled+rabyear+ragender +
               pc1_5a+pc1_5b+pc1_5c+pc1_5d+ pc1_5e+
               pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data= data_rG)

M_rGHealth <- lm(Health_AV_scaled~ Age_AV+ea_pgs3_edu3_ssgac18+rabyear+ragender 
                 +pc1_5a+pc1_5b+pc1_5c+pc1_5d+pc1_5e+
                   pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e, data=data_rG)

M_rHealth<- lm(Health_AV_scaled~ raedyrs_scaled+rabyear+ragender +
                 pc1_5a+pc1_5b+pc1_5c+pc1_5d+ pc1_5e+
                 pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data= data_rG)
```

```
install.packages("tidyverse") 
library(tidyverse) 
install.packages("broom") 
library(broom) 
install.packages("stargazer") 
library(stargazer)

# Plot results
stargazer(M_rGHeight, M_rHeight, M_rGBMI, M_rBMI, M_rGFert,M_rFert, 
          M_rGHealth, M_rHealth, type = "text", single.row=T)

m1 <- tidy(M_rGHeight) %>% mutate(var = "height", var2 = "genetic") %>% 
  filter(term == "ea_pgs3_edu3_ssgac18") %>%  
  mutate(up=estimate+std.error*1.96, low = estimate - std.error*1.96)

m2 <- tidy(M_rHeight) %>% mutate(var = "height", var2 = 'social') %>% 
  filter(term == 'raedyrs_scaled') %>% 
  mutate(up=estimate+std.error*1.96, low = estimate - std.error*1.96)

m3 <- tidy(M_rGBMI) %>% mutate(var = "bmi", var2 = 'gene- tic') %>% 
  filter(term == 'ea_pgs3_edu3_ssgac18')%>% 
  mutate(up=estimate+std.error*1.96, low = estimate - std.error*1.96)
 
m4 <- tidy(M_rBMI) %>% mutate(var = 'bmi', var2 = 'social') %>% 
  filter(term == 'raedyrs_scaled') %>% 
  mutate(up=estimate+std.error*1.96, low = estimate - std.error*1.96)
  
m5 <- tidy(M_rGFert) %>% mutate(var = 'fertility', var2 = 'genetic') %>% 
  filter(term == 'ea_pgs3_edu3_ssgac18') %>% 
  mutate(up=estimate+std.error*1.96, low = estimate - std.error*1.96)
                                                                                                   
m6 <- tidy(M_rFert) %>% mutate(var = 'fertility', var2 = 'social') %>% 
  filter(term == 'raedyrs_scaled') %>% 
  mutate(up=estimate+std.error*1.96,low = estimate - std.error*1.96)

m7 <- tidy(M_rGHealth) %>% mutate(var = 'health', var2 = 'genetic') %>% 
  filter(term == 'ea_pgs3_edu3_ssgac18') %>%
  mutate(up=estimate+std.error*1.96, low = estimate - std.error*1.96)
                                                                                                   
m8 <- tidy(M_rHealth) %>% mutate(var = 'health', var2 = 'social') %>%
  filter(term == 'raedyrs_scaled') %>% 
  mutate(up=estimate+std.error*1.96, low = estimate - std.error*1.96)

df <- rbind(m1,m2,m3,m4,m5,m6,m7,m8) 
df
df %>% ggplot(aes(x = var2, y = estimate)) +
  geom_bar(stat = 'identity') + facet_grid(~var) + theme_bw() +
  geom_errorbar(aes(ymin=low, ymax=up), 
                width=.1, position=position_dodge(0.1))
```

```
# Define the variables for subset
vars_confounding <- c("raedyrs","rafeduc",
"ea_pgs3_edu2_ssgac16", "ea_pgs3_edu3_ssgac18", "rabyear", "ragender",
"pc1_5a", "pc1_5b","pc1_5c","pc1_5d", "pc1_5e", 
"pc6_10a", "pc6_10b", "pc6_10c", "pc6_10d", "pc6_10e")

#  Create a subset of data
data_confounding <- hrs_GePhen_uni[vars_confounding]
                                                        
#  Intergenerational transmission educational attainment - jointly fitting paternal effects
M_intergen <- lm(raedyrs~ rafeduc+ rabyear+ragender +
pc1_5a+pc1_5b+pc1_5c+pc1_5d+pc1_5e +
pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data=data_confounding)

summary(M_intergen)

M_intergen_PGS3 <- lm(raedyrs~ rafeduc+ea_pgs3_edu3_ssgac18+ 
rabyear+ragender +pc1_5a+pc1_5b+pc1_5c+pc1_5d+pc1_5e+ 
  pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data=data_confounding)
summary(M_intergen_PGS3)
```

## 11.3 Gene-Environment Interaction

```
# Define and extract variables for analysis
bmi <-c("r1bmi","r2bmi","r3bmi","r4bmi","r5bmi",
        "r6bmi","r7bmi", "r8bmi","r9bmi","r10bmi","r11bmi","r12bmi")

hrs_GePhen_uni$BMI_AV = rowMeans(hrs_GePhen_uni[,bmi], na.rm = TRUE)
summary(hrs_GePhen_uni$BMI_AV)

age <-c("r1agey_b","r2agey_b","r3agey_b","r4agey_b","r5agey_b","r6agey_b",
        "r7agey_b","r8agey_b","r9agey_b","r10agey_b","r11agey_b","r12agey_b")

hrs_GePhen_uni$Age_AV = rowMeans(hrs_GePhen_uni[,age], na.rm = TRUE)
summary(hrs_GePhen_uni$Age_AV)

vars_GxE <- c("BMI_AV","Age_AV","ea_pgs3_bmi_giant15","rabyear","ragender",
              "pc1_5a","pc1_5b","pc1_5c","pc1_5d","pc1_5e",
              "pc6_10a","pc6_10b","pc6_10c","pc6_10d","pc6_10e")

data_GxE <- hrs_GePhen_uni[vars_GxE]
```

```
#  Estimate additive model of birth cohort and BMI PGS
M_G_bc <- lm(BMI_AV~ ea_pgs3_bmi_giant15 +rabyear + Age_AV + ragender +
                 pc1_5a+pc1_5b+pc1_5c+pc1_5d+pc1_5e+pc6_10a+pc6_10b+
              pc6_10c+pc6_10d+pc6_10e,data=data_GxE)
summary(M_G_bc)
```

```
# Estimate model of multiplicative interaction term of Birth year x BMI-PGS

M_Gxbc <- lm(BMI_AV~ ea_pgs3_bmi_giant15 +rabyear + ea_pgs3_bmi_giant15*rabyear +ragender+
               Age_AV+ +pc1_5a+pc1_5b+pc1_5c+pc1_5d+pc1_5e+
               pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data=data_GxE)
summary(M_Gxbc)
```

```
# Center interaction variables and re-estimate model 
install.packages('jtools') 
library(jtools)

data_GxE$rabyear_center <- center(data_GxE$rabyear) 
data_GxE$Age_AV_center <- center(data_GxE$Age_AV) 
data_GxE$Gender_center <- center(data_GxE$ragender)

M_Gxbc_center <- lm(BMI_AV~ ea_pgs3_bmi_giant15 +rabyear_center + 
                      ea_pgs3_bmi_giant15*rabyear_center + Age_AV_center 
                    + Gender_center +pc1_5a+pc1_5b+pc1_5c+pc1_5d+pc1_5e+
                    pc6_10a+pc6_10b+pc6_10c+pc6_10d+pc6_10e,data=data_GxE) 
summary(M_Gxbc_center)
```

```
#  Produce plot of GxE interaction
install.packages('arm')
library(arm) 
install.packages('lme4') 
library(lme4)
install.packages('interplot') 
library(interplot)

interplot(m = M_Gxbc_center, var1 = "ea_pgs3_bmi_giant15", 
          var2 = "rabyear_center")+
  xlab("Year born") +
  ylab("Estimated Coefficient for BMI genetic score") + theme_bw() +
  ggtitle("GxE interaction between PGS for BMI and birth
year") +
  theme(plot.title = element_text(face="bold"))+
  geom_hline(yintercept = 0, linetype = "dashed")
```

```
# Create interaction plot of predicted BMI by BMI-PGS x Birth Year
install.packages('jtools')
library(jtools) 
install.packages("interactions") 
library(interactions)

interact_plot(M_Gxbc, pred = "ea_pgs3_bmi_giant15", 
              modx = "rabyear", interval = TRUE,int.width = 0.95, 
              x.label="BMI-PGS", y.label="Predicted value BMI", 
              main.title = "Predicted BMI by birth year", legend.main="Birth year")
```

```
# Create interaction plot of predicted BMI by Birth Year x BMI-PGS
interact_plot(M_Gxbc, pred = "rabyear", modx = "ea_pgs3_bmi_giant15", 
              interval = TRUE, int.width = 0.95, x.label="Birth year", 
              y.label="Predicted value BMI", 
              main.title = "Predicted birth year by BMI-PGS", 
              legend.main="BMI-PGS")
```