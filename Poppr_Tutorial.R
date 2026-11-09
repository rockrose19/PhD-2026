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
