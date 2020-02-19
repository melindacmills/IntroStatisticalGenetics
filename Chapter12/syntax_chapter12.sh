###############################################################################################
###################                Syntax Chapter 12      ####################################
###############################################################################################


##### IN R



# Load the library
install.packages("qqman")
library(qqman)


# Use R to download the summary statistics
download.file("http://ssgac.org/documents/EduYears_Main.txt.gz", dest="EA2_results.txt.gz")

# Import the summary statistics in R
gwasResults<-read.table("EA2_results.txt.gz", header=T)


# save the figure into an external png file format
png(file="manhattan_without_highlights.png" , width = 1200, height = 600)

manhattan(gwasResults, chr="CHR", bp="POS", snp="MarkerName", p="Pval",suggestiveline=F)

 dev.off()




 gwasResults<-subset(gwasResults, Pval<0.05)
hits<-gwasResults[gwasResults$Pval<5e-08,]

# create new variable with SNPs that we want to highlight
gwasResults$highlight.snps<-0

 for ( i in 1:  dim(hits)[1]){
 
   chr<-hits[i,2]
   loc_min<- hits[i, 3]-5000
   loc_max<- hits[i, 3]+5000
  
neighbours.snps<-gwasResults$MarkerName[gwasResults$CHR==chr &
                                      gwasResults$POS>loc_min & 
                                       gwasResults$POS<loc_max]

gwasResults$highlight.snps[gwasResults$MarkerName %in%
                            neighbours.snps]  <- 1
 }




 png(file="manhattan_with_highlights.png" , width = 1200, height = 600)

 # add highlight command to the Manhattan plot 
  manhattan(gwasResults, chr="CHR", bp="POS", snp="MarkerName", p="Pval", highlight=gwasResults$MarkerName[gwasResults$highlight.snps==1],suggestiveline=F)
  dev.off()


qq(gwasResults$Pval)