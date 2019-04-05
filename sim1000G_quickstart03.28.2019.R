####  QUICKSTART OF S1000G SIMULATOR  #####
#
# PURPOSE: to generate variants in the region of the example vcf file for 200 unrelated
#          individuals
#
##CREATED: 3.28.2019
##LAST UPDATED: 3.28.2019



library("sim1000G")
library("gplots")



# Read the example file included in sim1000G

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir,"region.vcf.gz")


#run the simulation 
vcf = readVCF( vcf_file, maxNumberOfVariants = 600 , min_maf = 0.01, max_maf = 1)

startSimulation(vcf, totalNumberOfIndividuals = 1000)
ids = generateUnrelatedIndividuals(200)

genotype = retrieveGenotypes(ids)


#look at genotypes and compare allele frequencies with ones in original vcf file 
plot( apply(genotype,2,mean)/2 ,  apply(vcf$gt1+vcf$gt2,1,mean)/2 )
abline(0,1,lty=1,lwd=9,col=rgb(0,0,1,0.3))

#show generated genotypes 
gplots::heatmap.2(genotype,col=c("white","orange","red"),Colv=F, trace="none")

#compute correlation between the markers to show a linkage disequilibrium plot of the region

gplots::heatmap.2( cor(genotype)^2 , trace="none", col=rev(heat.colors(200)) ,Rowv=F,Colv=F )


