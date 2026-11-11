install.packages("poppr", repos = "https://cloud.r-project.org/")
library(poppr)

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


