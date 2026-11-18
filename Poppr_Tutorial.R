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
monpop
head(mll(monpop, "original"), 20) # Showing the definitions for the first 20 samples
#Get a warning because monopop data was created for mll() was implemented 

mll(monpop) <- "original"
monpop
# 264 original multilocus genotypes (current MLG number and definition)

mll(monpop) <- "custom"
monpop
# Now 264 custom MLGs.