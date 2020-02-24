###############################################################################################
###################                Syntax Chapter 10       ####################################
###############################################################################################


./plink 	--bfile 1kg_hm3_qc \
 		--score score_rs9930506.txt 1 2 3 \
 	--pheno BMI_pheno.txt \
 	--out 1kg_FTOscore


head 1kg_FTOscore.profile


Rscript PRSice.R --dir . \
    --prsice ./PRSice \
    --base BMI.txt \
    --target 1kg_hm3_qc \
    --snp MarkerName \
    --A1 A1 \
    --A2 A2 \
    --stat Beta \
    --pvalue Pval \
    --pheno-file BMI_pheno.txt \
    --bar-levels 1 \
    --fastscore \
    --binary-target F \
    --out BMI_score_all



Rscript PRSice.R --dir . \
    --prsice ./PRSice \
    --base BMI.txt \
    --target 1kg_hm3_qc \
    --thread 1 \
    --snp MarkerName \
    --A1 A1 \
    --A2 A2 \
    --stat Beta \
    --pvalue Pval \
    --bar-levels 5e-08,5e-07,5e-06,5e-05,5e-04,5e-03,5e-02,5e-01 \
    --fastscore \
    --all-score \
    --no-regress \
    --binary-target F \
    --out BMIscore_thresholds  


head BMIscore_thresholds.all.score

Rscript PRSice.R --dir . \
    --prsice ./PRSice \
    --base BMI.txt \
    --target 1kg_hm3_qc \
    --thread 1 \
    --snp MarkerName \
    --A1 A1 \
    --A2 A2 \
    --stat Beta \
    --pvalue Pval \
    --pheno-file BMI_pheno.txt \
    --interval 0.00005 \
    --lower 0.0001 \
    --quantile 20 \
    --all-score \
    --binary-target F \
    --out BMIscore_graphics

Rscript PRSice.R --dir . \
    --prsice ./PRSice \
    --base BMI.txt \
    --target 1kg_hm3_qc \
    --thread 1 \
    --snp MarkerName \
    --A1 A1 \
    --A2 A2 \
    --stat Beta \
    --pvalue Pval \
    --no-clump F \
    --pheno-file Obesity_pheno.txt \
    --interval 0.00005 \
    --lower 0.0001 \
    --quantile 5 \
    --all-score \
    --binary-target T \
    --out Obesity_score_graphics





###############################################################################################
###################     This part needs to be exectuted in R     ##############################
###############################################################################################


#import external data
data<-read.table("BMIscore_thresholds.all.score", header=T)

# show first rows of data
head(data)

# calculate new standardised variable
data$PGS=(data$X1-mean(data$X1))/sd(data$X1, na.rm=T)

#plot histogram of the polygenic score 
 hist(data$PGS)

# import external data with phenotype
pheno_BMI<-read.table("BMI_pheno.txt", header=T)

# merge the two datasets
data.with.pheno<-merge(data,pheno_BMI, by="IID")



# run a linear regression model
mod1<-(lm(BMI~PGS, data=data.with.pheno))
summary(mod1)


summary(mod1)$r.square




# create a vector with column names
columns=c("FID", "IID", "pca1", "pca2", "pca3", "pca4", "pca5", "pca6", "pca7", "pca8", "pca9", "pca10")

# read external data with PCAs
pca<- read.table(file="pca.eigenvec", sep = "", header=F, col.names=columns)[,c(2:12)]

# merge file with covariates with the rest of the file
data.with.covars<-merge(data.with.pheno,pca, by="IID")

# Estimate new linear regression model
mod2<-lm(BMI~PGS+pca1+pca2+pca3+pca3+pca5+pca6+pca7+pca8+pca9+pca10, data=data.with.covars)

# Calculate R-Square
summary(mod2)$r.square

# Estimate linear regression model
mod2.no.pgs<-update(mod2,  ~ . - PGS)

# Estimate difference in R-squared
print(summary(mod2)$r.square-summary(mod2.no.pgs)$r.square)
# install new package in R to perform bootstrap

# install new package in R to perform bootstrap

install.packages("boot")
library(boot)
set.seed(12345)

# define new function to calculate differential r-squared.
# The following commands define the new function
rsq <- function(formula, PGS=PGS, data, indices) {
  d <- data[indices,] 
  fit1 <- lm(formula,  data=d)
  fit2 <- update(fit1,  ~ . + PGS)
  
  return(summary(fit2)$r.square-summary(fit1)$r.square)
} 


# bootstrapping with 1000 replications 
results <- boot(data=data.with.covars, statistic=rsq, 
  	R=1000, formula=BMI~pca1+pca2+pca3+pca3+pca5+pca6+pca7+pca8+pca9+pca10, PGS=PGS)
boot.ci(results, type="norm")      


# import external data with multi-ancestry information

data<-read.table("BMIscore_MULTIANCESTRY.best", header=T)
head(data)

# Calculate  standardized PGS
data$PGS=(data$PRS-mean(data$PRS))/sd(data$PRS, na.rm=T)

data.with.pheno<-merge(data,pheno_BMI, by="IID")

#import dataset with geographical information (ancestral populations)
geo<-read.table(file="1kg_samples.txt", sep = "\t", header=T)[, c(1,4,5,6,7)]
names(geo)[1]<-"IID"

# Create a vector with new columns names
columns=c("FID", "IID", "pca1", "pca2", "pca3", "pca4", "pca5", "pca6", "pca7", "pca8", "pca9", "pca10")

# import exteral data with PCAs
pca<- read.table(file="pca.eigenvec", sep = "", header=F, col.names=columns)[,c(2:12)]

# merge information with covariates and ancestral information
data.with.covars<-merge(data.with.pheno,pca, by="IID")
data.with.geo<-merge(data.with.covars,geo, by="IID")

# fit a linear regression model
mod.AFR<-lm(BMI~PGS, data=subset(data.with.geo, Superpopulation.name=="African"))

summary(mod.AFR)$r.square


# run bootstrap function to  estimate increased R squared
results <- boot(data=subset(data.with.geo, Superpopulation.name=="African"), statistic=rsq, 
  	R=1000, formula=BMI~pca1+pca2+pca3+pca3+pca5+pca6+pca7+pca8+pca9+pca10, PGS=PGS)

boot.ci(results, type="norm")      



###############################################################################################
###################     This part needs to be exectuted in the bash    ##############################
###############################################################################################


python --version

head BMI_LDpred.txt
python ldpred/LDpred.py coord --help

python ldpred/LDpred.py coord --gf=1kg_EU_qc --ssf=BMI_LDpred.txt --ssf-format=STANDARD --N 500000 --out BMI_ldpred_coord    

python ldpred/LDpred.py gibbs --help 

python ldpred/LDpred.py gibbs --cf BMI_ldpred_coord --ldr 100 --ldf 1kg_bmi --f 1 --N=500000 --out BMI_weights_LD 



head BMI_weights_LD_LDpred_p1.0000e+00.txt

python ldpred/LDpred.py score --help


python ldpred/LDpred.py score \
--gf 1kg_EU_qc \
--rf BMI_weights_LD \
--pf BMI_pheno_LDpred.txt \
--rf-format LDPRED \
--out BMI_score_LD

head BMI_score_LD_LDpred_p1.0000e+00.txt


./plink --bfile 1kg_EU_qc \
        --score BMI_weights_LD_LDpred_p1.0000e+00.txt 3 4 7 header \
  --out scores_BMI_LD
   

