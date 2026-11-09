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
