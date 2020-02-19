# Syntax Chapter 9     

## Data used in this chapter:
* aaa
* bbb
* ccc


## 9.2 Association Analysis
```
./plink    	 --bfile 1kg_EU_BMI \
        	 --snps rs9674439 \
       	 --assoc \
      	 --linear \
      	 --out BMIrs9674439
```

```
./plink    	 --bfile 1kg_EU_Overweight \
        	 --snps rs9674439 \
       	 --assoc \
      	 --logistic \
      	 --out Overweight_rs9674439
```


```
./plink    	 --bfile 1kg_EU_BMI \
        	 --snps rs9674439 \
       	 --assoc \
      	 --linear dominant \
      	 --out BMIrs9674439
```		 

```		 
./plink    	 --bfile 1kg_EU_BMI \
       	 --assoc \
      	 --linear \
      	 --out BMIgwas
```		 


## 9.3 Linkage disequilibrium
```
./plink --bfile hapmap-ceu --ld rs2883059 rs2777888 --out ld_example
```

```
./plink 	  --bfile 1kg_hm3_qc --maf 0.01 \
        --indep-pairwise 50 5 0.2 \
        --out 1kg_hm3_qc_pruned
```

```
./plink	--bfile 1kg_hm3_qc \
	--extract 1kg_hm3_qc_pruned.prune.in \
	--make-bed \
 		--out 1kg_hm3_pruned
```

## 9.4 Population Stratification

```
./plink --bfile 1kg_hm3_pruned --pca 10 --out 1kg_pca
```



**This part needs to be exectuted in R**

### Panel A
```
library(ggplot2)
columns=c("fid", "Sample.name", "pca1", "pca2", "pca3", "pca4",   "pca5", "pca6", "pca7", "pca8", "pca9", "pca10")

pca<- read.table(file="1kg_pca.eigenvec", sep = "", header=F, col.names=columns)[,c(2:12)]

ggplot(pca, aes(x=pca1, y=pca2))+ geom_point()+ theme_bw()+  xlab("PC1") + ylab("PC2")
```


#### Panel B
```
geo <- read.table(file="1kg_samples.txt", sep = "\t",header=T)[,c(1,4,5,6,7)]


data <- merge(geo, pca, by="Sample.name")
ggplot(data, aes(x=pca1, y=pca2, col=Superpopulation.name))+
     	 geom_point()+ theme_bw()+
     	 xlab("PC1") + ylab("PC2")+
 	 labs(col = "")
```




## 9.5 Genetic Relatedness
```
./plink --bfile 1kg_hm3_pruned --keep 1kg_samples_EUR.txt  --distance --out ibs_matrix

./plink --bfile 1kg_hm3_pruned --keep 1kg_samples_EUR.txt --make-rel --out rel_matrix
```

```
./gcta64	--bfile 1kg_EU_BMI \
 	 	--autosome \
 	--maf 0.01 \
 	--make-grm \
 	--out 1kg_gcta 
```

```
./gcta64 --grm 1kg_gcta --grm-cutoff 0.025 --make-grm --out 1kg_rm025
./gcta64 --grm 1kg_rm025 --pheno BMI_pheno.txt --reml --out 1kg_BMI_h2
```