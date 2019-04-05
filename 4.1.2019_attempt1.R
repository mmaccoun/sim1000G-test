####  sim1000G-miniProjects  #####

#      *****OBJECTIVE*****      
#     gain better understanding of functions available in the simulator 
#   
# PROCEDURE : (1) generate a new nuclear family, using vcf data file included in sim1000G package to generate 
#                       a simulated genomic region of the genome 
#             (2) simulate 3 generations of the family, with 2 offspring for each generation
#             (3) compute the exact IBD1 & IBD2 proportions for siblings in in one generation, and offspring with 
#                       both parents  
#
# CONCLUSION : based on the computed IBD1 and IBD2 states, I am not confident of the abilities of this simulation
#              as seen in the data frame fam_dat, there are a total of 10 individuals in the 3 generation of this 
#              pedigree, and IBD1 and 2 states are not consistant, nor does it seem plausible for pairs of siblings
#              to have 0 IBD1 and IBD2 proportionsn, and the same is observed with offspring-parent relationships. 
#
#             Perhaps using a VCF file of our own will yield more promising results.   
#             Used own data file, also failure -- will DIY 
#

##CREATED: 4.1.2019
##LAST UPDATED: 4.1.2019 

library(sim1000G)

setwd("~/Desktop/IBDstates")

##from documentation, get input and needed map 

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path("~/Desktop/1000genomes", "chr22_100.vcf")
vcf = readVCF( vcf_file, maxNumberOfVariants = 100, min_maf = 0.01, max_maf = NA)

genetic_map_of_region = downloadGeneticMap(22, dir = NA)
readGeneticMapFromFile(genetic_map_of_region)

#will read data from sample vfc file provided, will allocate memory for max # 
            # of indivudals in this simulation to be 2000
startSimulation(vcf, totalNumberOfIndividuals = 2000) 

#generate a new nuclear family in the beginning of the simulation 
#code can be found in the documentation for function newNuclearFamily(family_id)
fam1 = newNuclearFamily(1)


#use this family and generate genotype data foe a family of 3 generations with 2 offspring
  # in generation 2 and 2 offspring in generation 3


#sex 1 = male, sex 2 = female, can track mother/father 

fam1_gen3 = newFamily3generations(1, 2, c(2,2))
print(fam1_gen3)

#store data in data frame 
fam_dat = data.frame(fam1_gen3)
#assign family members letters to 
fam_dat$letterID = c("a", "b", "d", "e", "i", "j", "c", "f", "g", "h")


# will compute the mean IBD12(haploid identical with maternal/paternal AND identical  
  # state) computed from shared haplotypes for the different individuals 

# ***PARENTAL **** 
#ca, offspring paternal comparison from generation 2 to 1 
ca = computePairIBD12(7,1)

#cb, offspring maternal comparison form generation 2 to 1 
cb = computePairIBD12(7,2)

#da, offspring paternal comparison from generation 2 to 1 
da = computePairIBD12(3,1)

#db, offspring maternal comparison form generation 2 to 1 
db = computePairIBD12(3,2)

#gc, offspring paternal comparison from generation 3 to 2 
gc = computePairIBD12(9 ,7)

#gf, offspring maternal comparison from generation 3 to 2 
gf = computePairIBD12(9, 8)

#hc, offspring paternal comparison from generation 3 to 2 
hc = computePairIBD12(10, 7)

#hf, offspring maternal comparison from generation 3 to 2 
hf = computePairIBD12(10, 8)

#id, offspring paternal comparison from generation 3 to 2 
id = computePairIBD12(5, 3)

#ie, offspring maternal comparison from generation 3 to 2
ie = computePairIBD12(5,4)

#jd, offspring paternal comparison from generation 3 to 2 
jd = computePairIBD12(6,3)

#je, offspring maternal comparison from generation 3 to 2 
je = computePairIBD12(6,4)

IBD12_offspringdat = data.frame(ca, cb, da, db, gc, gf, hc, hf, id, ie, jd, je)
 
# *** SIBLING ***

#cd, sibling pair from generation 2 
cd = computePairIBD12(7,3)

#hg, sibling pair from generation 3 
hg = computePairIBD12(10,9)
gh = computePairIBD12(9, 10) #see rearranging index wont change mean computed 

#ij, sibling pair from generation 3
ij = computePairIBD12(5,6)

IBD12_sibdat = data.frame(cd, hg, gh, ij)


