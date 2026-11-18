install.packages("poppr", repos = "https://cloud.r-project.org/")
library(poppr)

###Data import and manipulation  ####

find.package("poppr") 
# "C:/Users/tallu/AppData/Local/R/win-library/4.3/poppr"

#Run code below, go to location of poppr package, select files and within files, select rootrot.csv, it will be saved as x
x <- getfile()

myData  <- read.genalex(x$files)
myData


popdata <- poppr(myData)

#Look at the data 
popdata


#Example of GenAIEx format 
data(microbov)
strata(microbov) <- data.frame(other(microbov)) # set the strata
microbov
genind2genalex(microbov, file = "C:/Users/tallu/OneDrive/Desktop/microbov.csv")

#Import data in GENAIEX format 

rootrot <- read.genalex(system.file("files/rootrot.csv", package="poppr"))
rootrot

#Utilising the other slot example 

data(nancycats) # Load the data
other(nancycats) # geographical coordinates


repeats <- rep(2, nLoc(nancycats)) #nLoc = number of loci
repeats

other(nancycats)$repeat_lengths <- repeats
other(nancycats) # two items named xy and repeat_lengths

#Genclone example 

data(Aeut)
Aeut


#Set stratification 
strata(Aeut) <- other(Aeut)$population_hierarchy[-1]
Aeut
#Convert to genclone object 
agc <- as.genclone(Aeut)
agc
#Can also convert back to genind object

#Missing data example
library("poppr")
data(nancycats)
info_table(nancycats, plot = TRUE)

#locus fca8 has a lot of missing data. for an example, will use first 5 individuals at the first locus
tab(nancycats)[1:5, 8:13]

#Exclude loci and genotypes with >5% missing data

nanloci <-  missingno(nancycats, "loci")
nangeno <-  missingno(nancycats, "geno")
tab(nanloci)[1:5, 8:13]

#Example for extracting populations 

data("H3N2", package = "adegenet")
strata(H3N2) <- data.frame(other(H3N2)$x)
H3N2

#Set the population factor to country 

setPop(H3N2) <- ~country
popNames(H3N2) # Only two countries from North America.
v_na <- popsub(H3N2, sublist = c("USA", "Canada"))
popNames(v_na)

#If we want to see the population size, we can use the adegenet function nInd():

c(NorthAmerica = nInd(v_na), Total = nInd(H3N2))

#Dataset without North American countries 

v_na_minus <- popsub(H3N2, exclude = c("USA", "Canada"))
popNames(v_na_minus)

#Number of individuals in both sets equal to original data set

(nInd(v_na_minus) + nInd(v_na)) == nInd(H3N2)

#First 10 populations in alphabetical order, but no countries in North America

vsort <- sort(popNames(H3N2))[1:10]
vsort
valph <- popsub(H3N2, sublist = vsort, exclude = c("USA", "Canada"))
popNames(valph)

#Example for clone correcting data 

data(Aeut)
strata(Aeut) <- data.frame(other(Aeut)$population_hierarchy[-1])
Aeut

#when read in data with read.genalex(), default to represent it in gencline object. clonecorrect() works on multilocus genotype - more efficient to convert to genclone object 1st. 

aphan <- as.genclone(Aeut)
nameStrata(Aeut) <- ~field/sample #rename strata

clonecorrect(aphan,  strata = ~Pop/Subpop) # correct by sample w.r.t field and keep field as the population 

clonecorrect(aphan, strata = ~Pop) # correcting by field - MLG much closer to census

clonecorrect(aphan, strata = NA) #correcting over whole data set. MLG equal to census (diploid individuals)

#Shufflepop Bruvo's distance example

data(nancycats)
nan1 <- popsub(nancycats, 1)
reps <- rep(2, 9) # Assuming dinucleotide repeats.
observed <- mean(bruvo.dist(nan1, replen = reps))
observed


#Replicate resampling 999 times & create histogram to visualise distribution what would happen under different panmiia assumps

set.seed(9999)
bd.test <- replicate(999, mean(bruvo.dist(shufflepop(nan1, method = 2), replen = reps)))
hist(bd.test, xlab = "Bruvo's Distance", main = "Average Bruvo's distance over 999 randomizations")
abline(v = observed, col = "red")
legend('topleft', legend="observed", col="red", lty = 1)

#Uninformative microsat loci example

data(nancycats)
naninform <- informloci(nancycats, cutoff = 0.05)


#Determining how many MLGs in a dataset example

H3N2 #information displayed when view object
H3N2_mlg <- mlg(H3N2) #store number of MLGs

#MLGs across populations 

setPop(H3N2) <- ~country
v.dup <- mlg.crosspop(H3N2, quiet=TRUE) # quiet = FALSE gives output of MLG name, total individuals of MLG and populations where found 


#Count number of populations each MLG crosses 
v.num <- sapply(v.dup, length) # count the number of populations each MLG crosses.
head(v.num)

#MLG tables and graphs 

mlg.table(H3N2, sublist = "Norway", plot = TRUE)


setPop(H3N2) <- ~year
summary(H3N2) # Check the data to make sure it's correct.
library("vegan")
H.year <- mlg.table(H3N2, plot = FALSE)
rarecurve(H.year, ylab="Number of expected MLGs", sample=min(rowSums(H.year)),
          border = NA, fill = NA, font = 2, cex = 1, col = "blue")

#Combining MLG functions 

v.tab <- mlg.table(H3N2, plot = FALSE)
v.tab[1:10, 1:10] # Showing the first 10 columns and rows of the table.

setPop(H3N2) <- ~country
UGNN.list <- c("United Kingdom", "Germany", "Netherlands", "Norway")
UGNN <- mlg.crosspop(H3N2, sublist=UGNN.list, indexreturn=TRUE)


UGNN # Note that we have three numbers here. This will index the columns for us.
UGNN.list # And let's not forget that we have the population names.
v.tab[UGNN.list, UGNN]





####

###Analysis of Multilocus Genotypes and Lineages ####

library("poppr")
data(monpop)
head(mll(monpop, "original"), 20) # Showing the definitions for the first 20 samples
#Get a warning because monopop data was created for mll() was implemented 

mll(monpop) <- "original"
monpop
# 264 original multilocus genotypes (current MLG number and definition)

mll(monpop) <- "custom"
monpop
# Now 264 custom MLGs.

#Demonstration of how to contrust each MLG definition:

#Naive ("original")

mll(monpop) <- "original"
monpop
#There are 694 samples, but only 264 MLGs over 13 loci. 
#Naive definition takes strings of alleles and compares them for equality - doesn't take inot account genotyping error, hyper-variable loci and missing data - contribute to unique genotpe but not truly unique


#Example of 5 samples with alleles typed at 2 loci

grid_example <- matrix(c(1, 4,
                         1, 1,
                         5, 1,
                         9, 1,
                         9, 4),
                       ncol = 2,
                       byrow = TRUE)
rownames(grid_example) <- LETTERS[1:5]
colnames(grid_example) <- c("x", "y")
grid_example


x <- as.genclone(df2genind(grid_example, ploidy = 1))
tab(x)  # Look at the multilocus genotype table


nmll(x) # count the number of multilocus genotypes
mll(x)  # show the multilocus genotype definitions
#nmll() function counts the number of MLGs in sample
#mll() function displays the assignment of ML lineages in sample 

#Include another sample with partial missing data , missing y locus
x <- as.genclone(df2genind(rbind(grid_example, new = c(5, NA)), ploidy = 1))
tab(x)  # Note the missing data at locus 2.

nmll(x) #6 MLGs even though one contains missing data
mll(x)
#Include another sample with low-frequency allele at locus y "2", but allele "5" at locus x
grid_new <- rbind(grid_example,
                  new = c(5, NA),
                  mut = c(5, 2)
)
x <- as.genclone(df2genind(grid_new, ploidy = 1))
tab(x)
nmll(x)
mll(x)

#genotype added is closer to genotype that shares the 5 allele - happens because naive algorithm reports the rank of concatenated alleles

(xt <- apply(tab(x), 1, paste, collapse = ""))
rank(xt, ties.method = "first")

#Calculate pairwise euclidean distances between samples:
library("phangorn")
library("ape")
raw_dist <- function(x){
  dist(genind2df(x, usepop = FALSE))
}
(xdis <- raw_dist(x)) # new, mut and C are very similar 

plot.phylo(upgma(xdis)) # new, mut, C may be same multilocus lineage (MLL)

# Filtered (contracted)
# Fix issues with naive MLGs, can utlise genetic distance - allow us to collapse MLGs that are under a specific distance threshold 
#mlg.filter() creates a dissilimarity matrix from data and then filter based off matrix 
#need to define a threshold : minimum distance to consider 2 genotypes unique (anything below cluster into same MLG)

x # We have 7 MLGs before filtering

mlg.filter(x, distance = xdis) <- 1 + .Machine$double.eps^0.5
x # Now we have 5 MLGs

mll(x) <- "original" # We'll reset to the naive definition

#THRESHOLDS 
#To collapse all 3 new samples into the original 5 MLGs, had to set a threshold just above 1
#Use function to return the MLGs it created and the thresholds that were passed to clister MLLS:

mll(x, "original")
mlg.filter(x, distance = xdis, threshold = 1, stats = c("mlg", "thresholds"))
#MLG 5 has collapsed into MLG3 
#What happens when add very small number to 1:

(e <- .Machine$double.eps^0.5) # A very tiny number
mlg.filter(x, distance = xdis, threshold = 1 + e, stats = c("mlg", "thresholds"))

#Using threshold of 0 allows 2 MLGs seperated by zero genetic distance to remain mutually unique 
mlg.filter(x, distance = xdis, threshold = 0, stats = c("mlg", "thresholds"))
mll(x, "original")


#Setting the genotypes 
#Want the reassign sample genotypes to the collapsed version:
#This code is not working for some reason.....

x
mlg.table(x) # Before: 7 MLGs
mlg.filter(x, distance = xdis) <- 1 + e
x
mlg.table(x) # After: 5 MLGs
#See https://zkamvar.github.io/poppr/articles/mlg.html for this section above if needed 

#TIE BREAKERS (algorithms)
# Underlying mlg.filter are three algorithms that decide what genotypes go together- nearest neighbour, farthest neighbour, average neighbout (UPGMA)
# filter_stats() plots output of filtering algorithm 

data(Pinf)
Pinf
pinfreps <- fix_replen(Pinf, c(2, 2, 6, 2, 2, 2, 2, 2, 3, 3, 2))
pinf_filtered <- filter_stats(Pinf, distance = bruvo.dist, replen = pinfreps, plot = TRUE)
Pinf <- as.genclone(Pinf)
#Choosing a threshold 
#After choosing genetic distance and a filtering algorithm, must decide on threshold to represent min genetictance at which 2 individuals would be considered from different clonal lineages
#One method = look for an initial small peak in the histogram of pairwise genetic distances & set threshold between that peak and larger peak. Initial peak likely clones differentiated by a small set of random mutations. 
#If peak not obvious, then look for largest gap between all putative thresholds; use cutoff_predictor() with the output of filter_stats(). 

print(farthest_thresh <- cutoff_predictor(pinf_filtered$farthest$THRESHOLDS))
print(average_thresh  <- cutoff_predictor(pinf_filtered$average$THRESHOLDS))
print(nearest_thresh  <- cutoff_predictor(pinf_filtered$nearest$THRESHOLDS))

#Then can define multilocus lineages with following criteari: t = threshold, d= distance/Bruvo's Distance, a = algothithm (farthest neighbout)

mlg.filter(Pinf, distance = bruvo.dist, replen = pinfreps, algorithm = "f") <- farthest_thresh
Pinf

#Custom ("custom"), make custom MLGs
#Example:
data(partial_clone)
pc <- as.genclone(partial_clone)
mll(pc)

#Change MLGs to letters:
LETTERS[mll(pc)]  # The new MLGs
library(ggplot2)
mll.custom(pc) <- LETTERS[mll(pc)]
mlg.table(pc) #not working again, ggplot2 theme error

#Minimum spanning networks:
pcpal <- colorRampPalette(c("blue", "gold"))
set.seed(9001)
pcmsn <- bruvo.msn(pc, replen = rep(1, nLoc(pc)), palette = pcpal,
                   vertex.label.color = "firebrick", vertex.label.font = 2,
                   vertex.label.cex = 1.5)

#If we found Q is the same is M and wanted to change it:
mll.levels(pc)[mll.levels(pc) == "Q"] <- "M"

set.seed(9001)
pcmsn <- bruvo.msn(pc, replen = rep(1, nLoc(pc)), palette = pcpal,
                   vertex.label.color = "firebrick", vertex.label.font = 2,
                   vertex.label.cex = 1.5)


#DIVERSITY ANALYSIS
#Look at distribution of multilocus genotypes 
data(monpop)
splitStrata(monpop) <- ~Tree/Year/Symptom
montab <- mlg.table(monpop, strata = ~Symptom/Year, plot = FALSE) #again doesn't plot, turned off splot to do downstram stats

#Basic Statistics 
#The function diversity_stats() is used to get basic diversity statistics:

(monstat <- diversity_stats(montab))

#Confidence intervals 

diversity_ci(montab, n = 100L, raw = FALSE) #again, won't plot

