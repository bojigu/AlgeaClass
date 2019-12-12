# Creation of dataset for machine learning analysis of microbial flow cytometry data #

directory <- "D:\\Ronny_Proj\\Algae FACS raw data\\live_dead_staining_defined_ratios_of_living_and_dead_algae_FACS\\fac_data"
setwd(directory)
library(plyr)
library(reshape)
library(flowCore)
library(flowFP)
library(MASS)
library(h2o)
library(ggplot2)


if(file.exists("fset.RData")) {
  load("fset.RData")
} else {
  fset <- read.flowSet(path=directory, pattern="*.fcs")
  save(fset, file="fset.RData")
}

parse_name <- function(name) {
  csource <- substr(name, 1,3)
  day <- substr(name, 8,8)
  rep <- substr(name, 10,10)
  c(csource, day, rep, name)
}

if(!file.exists("samples.csv")) {
  samples <- sampleNames(fset)
  
  
  sample_table <- ldply(samples, parse_name)
  names(sample_table) <- c("csource", "day", "rep", "name")  
  write.csv(sample_table, "samples.csv", row.names=FALSE)
} else {
  sample_table <- read.csv("samples.csv")
}


#if(!file.exists("ValidationSet.RData")) {
split_flow_frame <- function(ff) {
  ## take a flow frame, extract the relevant columns, and 
  ## return it as a dataframe of 10 observations, each having 
  ## 1000 points of 3 readings each. 
  dat <- log(exprs(ff)[,c("FSC-A", "SSC-A", "PerCP-Cy5-5-A")])
  dat <- data.frame(dat)
  dat$mod <- 1:nrow(dat) %% 100
  
  split <- ddply(dat, .(mod), function(d) {
    c(d$'FSC.A', d$'SSC.A', d$'PerCP-Cy5-5.A')
  })
  
  split
}


bigset <- fsApply(fset, function(ff) {
  split <- split_flow_frame(ff) 
  split <- subset(split, select=-c(mod)) ## drop that column
  split$label <- parse_name(identifier(ff))[1]
  split$labelT <- paste0(parse_name(identifier(ff))[1], parse_name(identifier(ff))[2])
  split
})

bigset <- ldply(bigset, function(d) {d})
bigset <- bigset[, !(names(bigset) %in% c(".id"))]

datasplitter <- 1:nrow(bigset) # total number of rows .
datasplitter <- datasplitter %% 4; 

training_set <- bigset[which(datasplitter %in% c(0,1,2,3)),] ## Training set
validation_set <- bigset[which(datasplitter %in% c(2)),]
test_set <- bigset[which(datasplitter %in% c(3)),]

write.csv(training_set, file="trainingset.csv", row.names=F)
write.csv(validation_set, file="validationset.csv", row.names=F)
write.csv(test_set, file="testset.csv", row.names=F)
save(training_set, file="TrainingSet.RData")
save(validation_set, file="ValidationSet.RData")
save(test_set, file="TestingSet.RData")
#}










####PLAY WITH FLOWFP
mod <- flowFPModel(fset, name="FSC/SSC Model", parameters=c(1,4), nRecursions=7)
show(mod)
plot(mod)
fset
fp <- flowFP (fset, mod)
plot(fp, type="stack")
p <- flowFP (fset, param=c("FSC-A", "SSC-A"), nRecursions=8)
plex <- flowFPPlex()
for (levels in 8:5) {
  nRecursions(fp) <- levels
  plex <- append (plex, fp)
}
p <- flowFP (fset, param=c("FSC-A", "SSC-A"), nRecursions=7)
plex <- flowFPPlex()
for (levels in 7:5) {
  nRecursions(fp) <- levels
  plex <- append (plex, fp)
}
plot (plex, type="tangle", transformation="norm")
fp1 <- flowFP (fset, parameters=c("FSC-A", "SSC-A"), name="self model: fs1", nRecursions=7)
plot (fp, type="qc", main="Gate QC for Sample fs1")
fp <- flowFP (fset, parameters=c("FSC-A","SSC-A","PerCP-Cy5-5-A"), nRecursions=5)
plot (fp, type='plate')
nRecursions(fp)
counts(fp)
sampleNames(fp)
sampleClasses(fp)
parameters(fp)
tags(fp)
binBoundary(fp)
counts(fp)
fp <- flowFP (fset, parameters=c("FSC-A","SSC-A","PerCP-Cy5-5-A"), nRecursions=7)
counts(fp)
plot (fp, type='plate')
binBoundary(fp)
plot(fp)
plot(counts(fp))
a=counts(fp)
plot(a)
fp <- flowFP (fset, parameters=c("FSC-A","SSC-A","PerCP-Cy5-5-A"), nRecursions=5)
a=counts(fp)
plot(a)
matrix(a)
counts(fp)
