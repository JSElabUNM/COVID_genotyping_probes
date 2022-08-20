while(is.null(dev.list())==F) dev.off()
rm(list=ls())

###### run these if you need to install packages######
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install(version = "3.12")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("genefilter")
# 
# install.packages("xlsx")
# install.packages("writexl")
# install.packages("readxl")
# install.packages("rlist")
# install.packages("devtools")
# install.packages("Rtools")
# # install_github("easyGgplot2","kassambara")
# install.packages("tidyverse")
# install.packages("jsonlite")
# install.packages("futile.logger")
# install.packages("VennDiagram")`
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")
# install.packages("data.table")

# library(SummarizedExperiment)
library(data.table)
# library(writexl)
# library(readxl)
library(rlist)
library(devtools)
library(gtools)
library(ggplot2)
library(futile.logger)
library(stringr)
library(VennDiagram)

date <- Sys.Date()

#USAGE#
#cd /d H: or whatever drive/path it is all in
#need to make PCRXX_WYZZ.txt PCRXX_WY36_TV.txt first
#can highlight and ctrl+alt+enter to run from this window!

#in windows command prompt, quickly make list of all files using:
#dir /b >file.names.txt && Rscript replacement_B3.R file.names.txt

#or

#Rscript replacement_B3.R all_calls_11-Mar-2021_012_R_19_Clade.csv
#Rscript replacement_B3.R PCR32.csv
#setwd("C:/Users/Tara/Desktop/replacement_B3")
#raw = read.table("PCR32.tsv",header = T, sep="\t")
#raw = read.table("PCR32_WY26.csv",header = T, sep="\t")

#use file.names.txt for input to run all at once
#otherwise, input 1 file at a time



#don't add arguments for TV list or Variants list, may not get result you are looking for. 

####use if doing it within rstudio, with example with WY24
#x <- "PCR32_WY26.csv"
#file.name <- str_sub(x,1L,end = str_length(name.of.file)-4)
#x <- "PCR59_WY59.csv"
#x <- "PCR64_WY64.csv"
# raw = read.table(x,header = T, sep=",")
# file.name <- "WY24"
# name.of.file ="WY24.txt"
#test command line arguments


# setwd("H:Ryota School") #if running line by line, uncomment



export.fasta <- function(n,file.name.handle) 
{
  # n <- raw   #test for 
  # file.name.handle <- test.file.name
  
  # output.path <- ".\\"
  # full.output.path <- paste(output.path,"\\",file.name.handle,".fasta", collapse="",sep="")
  full.output.path <- paste(file.name.handle,".fasta", collapse="",sep="")
  header <- paste(">", file.name.handle, date)
  cat(header, file = full.output.path, append = FALSE)
  sink(file = full.output.path, append = TRUE)
  cat("\n",n, fill = 70,sep="")
  sink(file=NULL)
  
  # write.table(n, full.output.path, row.names = FALSE, quote = FALSE, append = TRUE)
  
}

export.data <- function(n,file.name.handle, notes = '', type) 
{ #get rid of column/row names if necessary
  #2/2/2020 working as intended
  #usage #export.data(data.frame_to_print, set.id = 'subset name',notes = 'any other notes')
  # #uncomment for test locally
  # n <- fixQ.highFIX.nonref #test for table #fix2$FIXCALL   #test for table
  # # n <- paste(n,collapse ="")
  # set.id <- "fixQ.highFIX.nonref"
  # notes <- "this is a test"  
  #easy to sort via ">" in word or text or ctrl+f
  
  #working on new version 3/12/2021
  #print out all position, ref, call, Qscore from MA.final
  
  
  # file.name.handle <- file.name
  # header <- c(">Sample\n", file.name.handle,"\nsubset.id", set.id, "\nread-count", dim(n)[1], "\nnotes",notes)
  
  header <- paste(">",date, file.name.handle, type,notes,"Calls: ", dim(n)[1],"\n") #  
  
  # output.path <- ".\\"
  # full.output.path <- paste(output.path,"\\",date,"_",file.name.handle,"_",type,".txt", collapse="",sep="")
  full.output.path <- paste(file.name.handle,"_",date,"_",type,".txt",collapse="",sep="")#
  sink(file = full.output.path)
  cat(header, file = full.output.path, append = FALSE)
  sink(file=NULL)
  write.table(n, full.output.path, row.names = FALSE, quote = FALSE, append = TRUE)
  
}

# x: the vector
# n: the number of samples
# centered: if FALSE, then average current sample and previous (n-1) samples
#           if TRUE, then average symmetrically in past and future. (If n is even, use one more sample from future.)
movingAverage <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}

sequence.sample <- function(x){
  
  #name.of.file <- "all_calls_11-Mar-2021_012_R_19_Clade.csv"
  #name.of.file <- input.argument
  name.of.file <- x #CHANGED so that it takes in single input for function sake
  file.name <- str_sub(name.of.file,1L,end = str_length(name.of.file)-4)
  cat("\nRunning sample:",file.name, sep="\t")
  # args <- commandArgs(trailingOnly = TRUE) #CHANGED
  args <- file.name
  # cat("\n",args, sep="\t")

if(length(args)==0){
  args <- c("/Users/jeremyedwards/Dropbox/JeremyEdwards/COVID/ARTIC2.26/all_calls_11-Mar-2021_012_R_19_Clade") #default runs this sample
}

if(length(args)>=2){
  name.of.TV.file <- paste(args[2],".csv",collapse="",sep="")
  TV = read.table(name.of.TV.file, header = T, sep="\t")
}

if(length(args)>=3){
  name.of.variants <- args[3]
} else{
  name.of.variants <- "VarStrains.csv"  #default opens pre-made VarStrains.csv
}
# name.of.TV <- paste(args[1], "_TV.txt",collapse="",sep="")
  # name.of.variants <- "variants_positions.txt"
raw = read.table(name.of.file,header = T, sep=",") #change sep="\t" or ","
vars = read.table(name.of.variants,header = T, sep=",") # or sep = ","
vars <- vars[,2:5]




if(length(args)>=2){
  TV <- TV$True.Variants
}
###start analyzing core0 data


colnames(raw)[1:6] = c("location","ref","core0","CALL12","Q","Q12") #correct column names as necessary

total <- dim(raw)[1]
Qth = 19.999999  #set lower bound of Q-score to use

#######table to quickly count contents

#core0
raw$Category <- 0
raw$Qscore <- 0
raw.call.reads <- sum(ifelse(raw$Q > Qth,1,0))    #Q > Qth, all calls
raw.temp <- raw[raw$Q > Qth,]
raw.temp$Category <- "Call"
raw <- raw[!(raw$location %in% raw.temp$location),]
raw <- rbind  (raw, raw.temp)

raw.noncall.reads <- sum(ifelse(raw$Q < Qth,1,0)) #number of low Q(noncalls) #226
raw.temp <- raw[raw$Q < Qth,]
raw.temp$Category <- "noncall"
raw <- raw[!(raw$location %in% raw.temp$location),]
raw <- rbind  (raw, raw.temp)
# unique(raw$Category)

raw.ref.call <- sum(ifelse(raw$Q > Qth & raw$core0 == raw$ref,1,0)) #high Q, refernce calls
raw.temp <- raw[raw$core0 ==  raw$ref & raw$Category == "Call",]
raw.temp$Category <- "Reference Call"
raw <- raw[!(raw$location %in% raw.temp$location),]
raw <- rbind  (raw, raw.temp)
# unique(raw$Category)

raw.nonref.calls <- sum(ifelse(raw$Q > Qth & raw$ref != raw$core0,1,0)) #highQ non reference (variant calls?)
raw.temp <- raw[raw$core0 !=  raw$ref & raw$Category == "Call",]
raw.temp$Category <- "Non-reference Call"
raw <- raw[!(raw$location %in% raw.temp$location),]
raw <- rbind  (raw, raw.temp)
# unique(raw$Category)

rm(raw.temp)

## output comment as necessary
# cat("\nCore0 data, table")
# table(raw$Category)
# 
# cat("\nJSE-fix data, table")
# table(raw$Category2)


#make function for generating fasta

#####recreate JSE-fix core1/2 replacement from spreadsheet in R  --> fix to replace lower Q with higher Q12 if core  --> new.data
#
######## replace non-ref core0 with ref core1/2 that is Q>Qth
Fix.list <- raw[which(raw$Q < 50 & raw$ref != raw$core0),]  #calls are non-ref and Q<50 (WY64 is 59) 
Fix.list <- Fix.list[which(Fix.list$ref == Fix.list$CALL12),] #19 CALL12 is reference reads
Fix.list <- Fix.list[which(Fix.list$Q12 > Fix.list$Q),] #of 15, 31 have Q12>Q (how many are calls)
Fix.list <- Fix.list[,1:6] #remove old fix columns
Fix.list <- cbind(Fix.list, Fix.list$CALL12,Fix.list$Q12) #generate new columns with calls to use (ie from core12, Q12)
names(Fix.list)[c(7,8)] <- c("FIXCALL","FIX")
number.replaced <- dim(Fix.list)[1] #how many core0 were replaced by core1/2

#fix is supposed to replicate jse/spreadsheet version
fix <- raw[-which(raw$location %in% Fix.list$location),c(1:6)] #remove Fix.list locations from raw data into fix
# length(unique(fix$location))
fix <- cbind(fix,fix[,c(3,5)]) #replace old fixcall/fixQ with the core0,Q of raw
names(fix)[c(7,8)] <- c("FIXCALL","FIX") #change name of core0,Q of [,7:8] 
fix <- rbind(fix, Fix.list) #append Fix.list to fix, ready for using known variants data

o <- order(fix$location) #order by location
fix <- fix[o,]  #fix should equal exactly the same as spreadsheet version
# sum(ifelse(fix$FIXCALL == raw$FIXCALL && fix$FIX == raw$FIX, 1,0)) #check for identical dataframes
rm(o,Fix.list)

#####fix/replace additional verified calls via core1/2 --> new.data
#replace low Q corresponding to core0 = ref with correct fixcall and highQ when core12 = fixcall, was correct call so was not filtered out in first pass
fix.new <- fix[(fix$Q12 > Qth & fix$Q12 > fix$FIX),] #Q12 > Qfix after initial pass of fixing using core1/2 

fix.new <- fix[(fix$ref == fix$CALL12 & #call12 = ref, reference reads
                  fix$Q12 > Qth & #Q12 > Q, how many are call
                  fix$Q12 > fix$FIX & #Q12 > FIX, how many have higher Q12
                  fix$FIXCALL == fix$ref),] #FIXCALL = ref, how many FIXCALL are reference
fix.new$FIX <- fix.new$Q12 #replace FIX with Q12 which is larger, all FIXCALL = reference

new.data <- fix[!(fix$location %in% fix.new$location),]
new.data <- rbind(new.data,fix.new) #append Fix.new.list to new.data, ready for using known variants data

o <- order(new.data$location) #order by location
new.data <- new.data[o,]   #replaced all calls where CALL12 = FIXCALL and reference, Q12 > FIX > Q12 
#new.data ready for categories

#test new.data for variant calls and higher Q
temp <-  new.data[new.data$core0 == new.data$CALL12,]
# sum(ifelse(new.data$core0 != new.data$ref,1,0))
temp2 <- temp[temp$core0 != temp$ref,] #core0=core12 but not reference
rm (temp,temp2)

#make category/Q with new.data
new.data$Category <- 0
new.data$Qscore <- 0
new.temp <- new.data[new.data$FIX > Qth,]
new.temp$Category <- "Call"  #name all with FIX > Qth as "CALLS"
new.data<- new.data[!(new.data$location %in% new.temp$location),]
new.data <- rbind  (new.data, new.temp)
# unique(new.data$Category)

new.temp <- new.data[new.data$FIX < Qth,]
new.temp$Category <- "noncall" #make FIX < Qth as "noncall"
new.data <- new.data[!(new.data$location %in% new.temp$location),]
new.data <- rbind  (new.data, new.temp)
# unique(new.data$Category)

new.temp <- new.data[new.data$FIXCALL == new.data$ref & new.data$Category == "Call",]
new.temp$Category <- "Reference Call" #name all reference calls 
new.data <- new.data[!(new.data$location %in% new.temp$location),]
new.data <- rbind  (new.data, new.temp)
# unique(new.data$Category)

new.temp <- new.data[new.data$FIXCALL != new.data$ref & new.data$Category == "Call",]
new.temp$Category <- "Non-reference Call" #name all non-reference calls 
new.data <- new.data[!(new.data$location %in% new.temp$location),]
new.data <- rbind  (new.data, new.temp)
# unique(new.data$Category)

rm(fix.new,new.temp)

o <- order(new.data$location) #order by location
new.data <- new.data[o,]   #replaced all calls where CALL12 = FIXCALL and reference, Q12 > FIX > Q12 

# ## output comment as necessary
# cat("\nJSE fix + replacement data, table")
# table(new.data$Category)
# table(raw$Category)

####################new.data analysis/check section
# length(unique(new.data$location)) #check for no duplicates during replacement
# length(unique(raw$location)) #check for no duplicates during replacement

# sum(ifelse(new.data$FIX < Qth ,1,0)) #all non-calls FIX < Q1th
# sum(ifelse(new.data$Category == "noncall",1,0))
# sum(ifelse(new.data$FIXCALL == new.data$ref & new.data$FIX > Qth,1,0)) #all reference calls
# sum(ifelse(new.data$Category == "Reference Call",1,0))
# sum(ifelse(new.data$ref != new.data$FIXCALL & new.data$FIX > Qth,1,0)) #all non-ref calls
# sum(ifelse(new.data$Category == "Non-reference Call",1,0))

##### looking at TVs of new.data
##### additional replacement of FIX by reference calls of core12 and Q12 > Q > Qth
#############end looking at FIX/FIXCALL



#make sure TV_JSE are correct variants, not miscalling




# new.data[new.data$core0 == new.data$CALL12 & new.data$FIXCALL != new.data$ref & new.data$FIX > Qth,]


######################################## try with and without this next section
#check for Q12 > FIX
new.data.temp <-new.data[new.data$Q12 > new.data$FIX,] #Q12> Q, 6683

new.data.temp.1 <- new.data.temp[new.data.temp$core0 == new.data.temp$ref,] #core0 = ref
new.data.temp.1.1 <- new.data.temp.1[new.data.temp.1$core0 == new.data.temp.1$CALL12,] #same as above but also core0=CALL12, ref=core0=core12
new.data.temp.1.1$FIX <- new.data.temp.1.1$Q12  #1.1 done, core12= ref, higher Q12 than Q 

new.data.temp.2 <- new.data.temp[new.data.temp$core0 != new.data.temp$ref,] #core0 != ref
new.data.temp.2.1 <- new.data.temp.2[new.data.temp.2$CALL12 == new.data.temp.2$ref,] #if CALL12=ref
new.data.temp.2.2 <- new.data.temp.2[new.data.temp.2$CALL12 == new.data.temp.2$core0,] #if CALL12=core0
new.data.temp.2.2$FIX <- new.data.temp.2.2$Q12 #2.2 done, CALL12 confirm core0 and higher Q12 than Q

new.data <- new.data[!(new.data$location %in% new.data.temp.1.1$location),]
new.data <- rbind(new.data, new.data.temp.1.1)

new.data <- new.data[!(new.data$location %in% new.data.temp.2.2$location),]
new.data <- rbind(new.data, new.data.temp.2.2)

rm(new.data.temp,new.data.temp.1, new.data.temp.2,new.data.temp.2.1,new.data.temp.1.1,new.data.temp.2.2)

#check for Q12> FIX and CALL12 == ref
new.data.temp<- new.data[new.data$Q12 > new.data$FIX & new.data$Q12 > new.data$Q,] #there should be no non-reference calls
new.data.temp<- new.data.temp[new.data.temp$CALL12 == new.data.temp$ref & new.data.temp$location >0,]
rm(new.data.temp)
#data check section for new.data data.frame
# unique(new.data$Category)
# sum(ifelse(new.data$FIX < Qth ,1,0)) #all non-calls FIX < Q1th
# sum(ifelse(new.data$Category == "noncall",1,0))
# sum(ifelse(new.data$FIXCALL == new.data$ref & new.data$FIX > Qth,1,0)) #all reference calls
# sum(ifelse(new.data$Category == "Reference Call",1,0))
# sum(ifelse(new.data$ref != new.data$FIXCALL & new.data$FIX > Qth,1,0)) #all non-ref calls
# sum(ifelse(new.data$Category == "Non-reference Call",1,0))
# sum(ifelse(new.data$ref != new.data$FIXCALL & new.data$FIX > Qth,1,0)) - dim(TV_JSE)[1] #all non-ref calls
# sum(ifelse(new.data$Category == "Non-reference Call" & new.data$TV != 'TV',1,0))
# dim(TV_JSE)[1] #TV from list
# sum(ifelse(new.data$Category == "Non-reference Call" & new.data$TV == 'TV',1,0))

o <- order(new.data$location) #order by location
new.data <- new.data[o,]   #replaced all calls where CALL12 = FIXCALL and reference, Q12 > FIX > Q12 

rm(o)


#need to determine which Var from current 3 to use 

TV_JSE <- new.data[new.data$location %in% vars$one,]
if(length(args)>=2){
  TV_JSE <- TV_JSE[TV_JSE$location %in% TV,]
  temp <- vars[vars$WY %in% TV_JSE$location,]
  TV_JSE_check <- cbind(TV_JSE, temp$alt1, temp$alt2)
  TV_JSE_check$TV = "TV"
  TV_JSE_check2 <- TV_JSE_check[(TV_JSE_check$FIXCALL == TV_JSE_check$`temp$alt1` | TV_JSE_check$FIXCALL == TV_JSE_check$`temp$alt2`  ),] #seems to verify and work correctly so far (9 in 24)
  TV_JSE_check <- TV_JSE_check[,c(1:10,13)]
  new.data$TV <- 0  #make new column for TV
  
  new.data <- new.data[!(new.data$location %in% TV_JSE$location),]
  new.data <- rbind(new.data, TV_JSE_check)
  # new.data[new.data$TV == "TV",] #long-form list of all TVs
  # new.data[new.data$Category == "Non-reference Call" & new.data$FIX > Qth,] #long-form list of all variant calls
  
  cat("\nCore0 TVs",
      "\nShort-read confirmed TVs\t",sum(ifelse(new.data$TV == "TV",1,0)),
      "\nTVs called\t",sum(ifelse(new.data$TV == "TV" & new.data$core0 != new.data$ref,1,0)),
      "\nTVs not called\t", sum(ifelse(new.data$TV == "TV" & new.data$core0 == new.data$ref,1,0)))
  
  cat("\nCore12 TVs",
      "\nShort-read confirmed TVs\t",sum(ifelse(new.data$TV == "TV",1,0)),
      "\nTVs called\t",sum(ifelse(new.data$TV == "TV" & new.data$FIXCALL != new.data$ref,1,0)),
      "\nTVs not called\t", sum(ifelse(new.data$TV == "TV" & new.data$FIXCALL == new.data$ref,1,0)))
  
} else{

  TV_JSE <- TV_JSE[TV_JSE$core0 == TV_JSE$CALL12 & TV_JSE$Category != "Reference Call",]
  cat("\nNo list of TVs were provided for", file.name, "\nEstimated True Variants through genotyping probeset:", dim(TV_JSE)[1])
}

o <- order(new.data$location) #order by location
new.data <- new.data[o,]

## output comment as necessary
# cat("\nJSE-fix data + replacement + TV, table")
# table(new.data$Category)
# 
# table(new.data$TV)



####MA-filter generation
ref.match <- new.data[which(new.data$location > 25 & new.data$FIX >0),]#filter for low Q and post >25 bc we only have pos13+

#MA.0.pos <- ref.match$location #correct calls but low Q position

#MA.matrix <- vector()
#i = 1
#j = 1
#pos = MA.0.pos[i]
#cat("\nCalculating MA now... This may take a few minutes\n")
#while(i<=length(MA.0.pos)-12){ #adjust here to make the last 12 bases cutoff, i<= length(MA.0.pos)-12
#  #can add something to look at the 12 bases at the end and beginning to see the errors/correct calls 
#  j = MA.0.pos[i]-12
#  z = new.data[which(new.data$location == j),]
#  old <- ref.match[ref.match$location == pos,]
#  MA.stdev <-vector()
#  MA.original.Q <- z$Q
#  while(j <= MA.0.pos[i]+12){
#    z = new.data[which(new.data$location == j),]
#    MA.stdev <- c(MA.stdev,z$Q)
#    
#    j = j+1
#  }
#  MA.stdev <- as.vector(MA.stdev)
#  MA.matrix.row <- c(pos, mean(MA.stdev),sd(MA.stdev),old$FIX) #row vector with position, moving average, MA-sdev
#  MA.matrix <- rbind(MA.matrix,MA.matrix.row) #merge row vector from above into larger matrix/table
#  
#  i <- i+1
#  pos = MA.0.pos[i]
#}
#cat("Calculations complete\n")
#
#MA.data.frame = as.data.frame(MA.matrix) #convert to formal data frame
#names(MA.data.frame)[c(1,2,3,4)] = c("location", "MA-Q-score", "stdev", "FIX")
# factor <- MA.data.frame$`MA-Q-score`/MA.data.frame$FIX  #factor no longer used as a metric
# MA.data.frame <- cbind(MA.data.frame, factor=factor) #factor looking at how much larger MA-Qscore is compared to fixQ
# cat("\nNew MA-Q cutoff factor\t", sum(ifelse(MA.data.frame$factor > mean(MA.data.frame$factor)+sd(MA.data.frame$factor*2),1,0))) #15, maybe do a t-test or alpha to look at significance?
#
#MA.data.frame <- MA.data.frame[!is.na(MA.data.frame$location),] #get rid of any error NAs
#loc <- MA.data.frame$location


MALoc <- fix$location[26:length(fix$FIX)-12] #JSE Added

#adding MA data, MA-Q-score, stdev,FIX(ie FIXQ) to original

MA <- new.data[which(new.data$location %in% MALoc),]  #cuts off first 12 bases but it should not be a big deal

# MA.data <- new.data[-which(new.data$location %in% MA$location),] #remove Fix.list locations from  new.data

MAQ <- movingAverage(fix$FIX[26:length(fix$FIX)-12],25,TRUE) #JSE Added
MA <- cbind(MA, MAQ, MALoc) #JSE Added

#MA <- cbind(MA, MA.data.frame$`MA-Q-score`,MA.data.frame$stdev)  #very interesting data! high MA lower std, lower MA high std etc, no longer adding factor (,MA.data.frame$factor)

o <- order(MA$location) #order by location
MA <- MA[o,]
if(length(args) >= 2){
  names(MA)[c(12,13)]=c("MAQ","MAstdev")  #MAs all added factor no longer used,"Factor"
} else{
  names(MA)[c(11,12)]=c("MAQ","MAstdev")  #MAs all added factor no longer used,"Factor"
}
  
MA.original <- MA #backup MA into new data.frame incase we need to reference back to it
MA <- MA.original
# rm(MA.data.frame, MA.matrix, old, z,  i, j, loc, MA.0.pos, MA.matrix.row, MA.original.Q, MA.stdev, pos)# factor no longer used (,factor)

#######################################################################
#######################################################################


#####MA filter
#############################MA FILTER
#30_24 has 10 non-reference variants after removing TV
#from MA.calls.nonref.woTV, check for FIX >20 but <30
#there are 8 in 30_24
#add in 'in filter' tag 
MA$MAremove <-0
MA$filter <- 0
# sum(ifelse(MA$Category == "Non-reference Call" & MA$FIX > 20 & MA$FIX < 30,1,0))
temp <- MA[MA$Category == "Non-reference Call" & MA$FIX > 20 & MA$FIX < 30,]
temp$filter <- "filtered"

MA <- MA[!(MA$location %in% temp$location),]
MA <- rbind(MA, temp)
# unique(MA$filter)

MA.filter.20_30 <- MA[MA$FIX > 20 & MA$FIX < 30 & MA$Category == "Non-reference Call" ,] #from 19 non-ref calls > 8
MA.filter.20_30$Category <- 'in filter'
MA.filter.20_30$Qscore <- MA.filter.20_30$FIX
MA.filter.20_30$filter <- "filtered"

#from MA.filter.20_30, if MAQ is less than mean(MAQ)-2*sd(MAQ) then put into MA.filter.remove
#add 'MAQ filter removed, Qs' and FIX value tag for ggplot

temp <- temp[temp$MAQ < mean(MA$MAQ) - 2 * sd(MA$MAQ),]
if(dim(temp)[1] > 0){
  temp$MAremove <-1 
  MA <- MA[!(MA$location %in% temp$location),]
  MA <- rbind(MA, temp)
  # unique(MA$filter)
} 

#make filtered list for combining for graphics
MA.filter.remove <- MA.filter.20_30[MA.filter.20_30$MAQ < mean(MA$MAQ) - 2 * sd(MA$MAQ),] #filter, 20<Q<30 and MAQ less than mean(MAQ)-2*sd(MAQ)

if(dim(MA.filter.remove)[1] > 0){
  MA.filter.remove$Category <- 'MAQ filter removed, Qs'
  Venn_filtered2 <- dim(MA.filter.remove)[1]
  MA.filter.remove$Qscore <- MA.filter.remove$FIX
  MA.filter.remove$filter <- "filtered"
  MA.filter.remove$MAremove <- 1
} else{
  MA.filter.remove <- MA.filter.20_30
  MA.filter.remove$Category <- 'MAQ filter removed, Qs'
  Venn_filtered2 <- dim(MA.filter.remove)[1]
  MA.filter.remove$Qscore <- 0
}

#make exact copy of MA.filter.remove but with MAQ displayed instead of FIX
MA.filter.remove.MAQ <- MA.filter.remove
MA.filter.remove.MAQ$Category <- 'MAQ filter removed, MAQs'
MA.filter.remove.MAQ$Qscore <- MA.filter.remove$MAQ

temp <- MA[MA$Category == "Non-reference Call" & MA$FIX > 20 & MA$FIX < 30 & MA$MAQ > mean(MA$MAQ) - 2 * sd(MA$MAQ),]
if(dim(temp)[1] > 0){
  temp$MAremove <-0 
  temp$filter <-1 
  MA <- MA[!(MA$location %in% temp$location),]
  MA <- rbind(MA, temp)
  # unique(MA$filter)
} 


#do opposite and filter the non-removed non-variant where MAQ is > mean(MA$MAQ) - 2 * sd(MA$MAQ)
MA.filter.nonremove.MAQ <- MA.filter.20_30[MA.filter.20_30$MAQ > mean(MA$MAQ) - 2 * sd(MA$MAQ),]

if(dim(MA.filter.nonremove.MAQ)[1] > 0){
  MA.filter.nonremove.MAQ$Category <- 'high MAQ filter, not removed, Qs'
  MA.filter.nonremove.MAQ$Qscore <- MA.filter.nonremove.MAQ$MAQ
} else{
  MA.filter.nonremove.MAQ <- MA.filter.20_30
  MA.filter.nonremove.MAQ$Category <- 'MAQ filter removed, Qs'
  MA.filter.nonremove.MAQ$Qscore <- 0
}

MA.filter.nonremove.Q <- MA.filter.20_30[MA.filter.20_30$MAQ > mean(MA$MAQ) - 2 * sd(MA$MAQ),]

if(dim(MA.filter.nonremove.Q)[1] > 0){
  MA.filter.nonremove.Q$Category <- 'high MAQ filter, not removed, Qs'
  MA.filter.nonremove.Q$Qscore <- MA.filter.nonremove.MAQ$FIX
  MA.filter.nonremove.MAQ <- MA.filter.nonremove.Q
  MA.filter.nonremove.MAQ$Category <- 'high MAQ filter, not removed, MAQs'
  MA.filter.nonremove.Q$Qscore <- MA.filter.nonremove.MAQ$MAQ
} else{
  MA.filter.nonremove.Q <- MA.filter.20_30
  MA.filter.nonremove.Q$Category <- 'MAQ filter removed, Qs'
  MA.filter.nonremove.Q$Qscore <- 0
}

# sum(ifelse(MA$Category == "Reference Call",1,0)) + sum(ifelse(MA$Category == "Non-reference Call",1,0))
MA.calls <- MA[MA$FIX > Qth,]  #use for accuracy calculations 29650

MA.noncalls <- MA[MA$FIX <Qth,] #these are removed from further use, maybe interesting to see 170
MA.noncalls$Qscore <- MA.noncalls$Q
temp.MA <- MA.noncalls[MA.noncalls$CALL12 == MA.noncalls$ref & MA.noncalls$Q12 > Qth,]

# sum(ifelse(MA$Category == "noncall",1,0))

#from MA.calls, take those that are refernece, 
#add 'ref.calls2' tag and 'FIX' values in 2 new columns for graph
#2 necessary for ordering in for ggplot
# sum(ifelse(MA$Category == "Reference Call",1,0))
MA.calls.ref <- MA.calls[MA.calls$FIXCALL == MA.calls$ref,] #use for accuracy calc and final cbind 29631
MA.calls.ref$Category <- "Reference Calls, Q" #All refernece calls
MA.calls.ref$Qscore <- MA.calls.ref$FIX

#from MA.calls, take those that are not-refernece, 
#from MA.calls.nonref, check for TV from TV_JSE by location
#add 'ref.calls1.tv' tag and 'FIX' values in 2 new columns for graph
#2 necessary for ordering in for ggplot
# sum(ifelse(MA$Category == "Non-reference Call",1,0))
MA.calls.nonref <- MA.calls[!(MA.calls$FIXCALL == MA.calls$ref),] #19
MA.calls.nonref$Category <- 'Non-Reference Calls'
Venn_nonref <- dim(MA.calls.nonref)[1]
MA.calls.nonref$Qscore <- MA.calls.nonref$FIX

# unique(MA$TV)
# sum(ifelse(MA$Category == "Non-reference Call" & MA$TV == "TV",1,0))
if(length(args)==1){
  cat("\nUsing potential TVs, no list provided by user")
}
MA.calls.tv <- MA.calls.nonref[MA.calls.nonref$location %in% TV_JSE$location,] 
MA.calls.tv$Category <- "True Variants" 
Venn_tv <- dim(MA.calls.tv)[1]
MA.calls.tv$Qscore <- MA.calls.tv$FIX


# ifelse(MA.filter.checker >0, MA.filter.nonremove.MAQ$Category <- 'high MAQ filter, not removed, MAQs', MA.filter.nonremove.MAQ <- MA.filter.remove.MAQ[1,])
# ifelse(MA.filter.checker >0, MA.filter.nonremove.MAQ$Qscore <- MA.filter.nonremove.MAQ$MAQ, MA.filter.nonremove.MAQ <- MA.filter.remove.MAQ[2,])
# 
# ifelse(MA.filter.checker == 0, MA.filter.nonremove.MAQ$Category <- 'high MAQ filter, not removed, MAQs','')
# ifelse(MA.filter.checker == 0, MA.filter.nonremove.MAQ$Qscore <- 0,'')
# ifelse(MA.filter.checker == 0, MA.filter.nonremove.MAQ[1] <- "",'')

# colnames(MA.filter.nonremove.MAQ) <- c('location','ref','core0','CALL12','Q','Q12','FIXCALL','FIX','MAQ','MAstdev','Factor','V12','V13')
# if(dim(MA.filter.remove.MAQ)[1]>0 ){
#   cat("\nreads removed by MA filter\t",dim(MA.filter.remove)[1],"\n")
# }else{
#   cat("\nNo reads removed by MA filter\t ")
# }

if(sum(ifelse(MA$filter =="filtered" & MA$MAremove ==1 ,1,0))>0 ){
  cat("\nreads removed by MA filter\t",sum(ifelse(MA$filter =="filtered" & MA$MAremove ==1 ,1,0)),"\n")
}else{
  cat("\nNo reads removed by MA filter\t ")
}

#check to see if any removed are in TV_JSE
if(sum(ifelse(MA$MAremove ==1 & MA$TV=="TV",1,0))>0){
  cat("\nConfirmed TV filtered by MA-filter\t", sum(ifelse(MA$MAremove ==1 & MA$TV=="TV",1,0))>0)
  mis.TV <- MA$location[MA$MAremove ==1 & MA$TV=="TV",]
  cat("\nPosition of removed TV(s) are\t", noquote(paste(MA$location[MA$MAremove ==1 & MA$TV!="TV"],sep="\t")))
}else{
  cat("\nNo TVs were wrongfully removed via MA filter\t")
}
MA.check.TV <- MA.filter.remove[MA.filter.remove$location %in% TV_JSE$location,]
if(dim(MA.check.TV)[1]>0){cat("\nConfirmed TV filtered by MA-filter\t", dim(MA.check.TV))
  Venn_TV_filtered = dim(MA.check.TV)[1]
} else{
  cat("\nNo TVs were wrongfully removed via MA filter\t")
  Venn_TV_filtered = 0
}

#remove is going from calls to non-calls, regardless of reference or non-refence calls
#remove those in MA.filter.remove from MA.calls.nonref.woTV
#MA.calls.nonref.woTV has 10 before filter
#MA.calls.nonref.woTV has 3 remaining non-reference calls after filter


############# MA filter data check
# sum(ifelse(MA$MAremove == 1,1,0))
# sum(ifelse(MA$Category == "Non-reference Call"  & MA$MAremove == 0,1,0)) #non-reference variant calls
# sum(ifelse(MA$Category == "Non-reference Call" & MA$MAremove == 1,1,0)) #non-ref calls removed by MA
# sum(ifelse(MA$Category == "Non-reference Call" & MA$MAremove == 0 & MA$TV == "TV",1,0)) #'correctly' called TVs
# sum(ifelse(MA$Category == "Non-reference Call" & MA$MAremove == 1 & MA$TV == "TV",1,0)) #'wrongly' removed TV calls
# sum(ifelse(MA$Category == "Non-reference Call" & MA$MAremove == 0 & MA$TV != "TV",1,0)) #non-ref calls actual errors (or new variant)
# sum(ifelse(MA$Category == "Non-reference Call" & MA$MAremove == 1 & MA$TV != "TV",1,0)) #non-ref calls, not TV, removed 'correctly'


###Categorical breakdowns for MA filter pulls
# sum(ifelse(MA$Category == "Non-reference Call" & MA$MAremove == 0,1,0)) #25
# MA[MA$Category == "Non-reference Call" & MA$MAremove == 0  ,]
# sum(ifelse(MA$Category == "Non-reference Call" & MA$filter =="filtered",1,0)) #15
# MA[MA$Category == "Non-reference Call" & MA$filter =="filtered"  ,]
# sum(ifelse(MA$filter =="filtered",1,0)) #15
# MA[ MA$filter =="filtered"  ,]
# sum(ifelse(MA$Category == "Non-reference Call" & MA$MAremove == 1,1,0)) #9
# MA[MA$Category == "Non-reference Call" & MA$MAremove == 1  ,] #9
# sum(ifelse(MA$Category == "Non-reference Call" & MA$MAremove == 0 & MA$TV == "TV",1,0)) #7
# MA[MA$Category == "Non-reference Call" & MA$MAremove == 0 & MA$TV == "TV",]

MA.calls.nonref.postMAF <- MA.calls.nonref[!(MA.calls.nonref$location %in% MA.filter.remove$location),] #use for cbind final list
# MA.calls.nonref.postMAF <- MA.calls.nonref[MA.calls.nonref$location %in% MA.filter.nonremove.MAQ$location | !(MA.calls.nonref$location %in% MA.filter.remove$location),]

## output comment as necessary
# cat("\nVariant calls\t",dim(MA.calls.nonref.postMAF)[1])

###################################################### everything done except verification of TV
#verify variant calls with TV list
# MA.calls.nonref.TVcheck <- MA.calls.nonref.postMAF[MA.calls.nonref.postMAF$location %in% TV_JSE$location,]
# MA.calls.nonref.TVcheck$Category <- "True Variants"
# 
# MA.calls.nonref.postTVcheck <- MA.calls.nonref.postMAF[!(MA.calls.nonref.postMAF$location %in% MA.calls.nonref.TVcheck$location),]


#combine reference calls w/ TV, non-reference calls w/o TV, non-calls, and one of filter.removes (this accounts for all reads in MA)
#also add MA.filter.remove.MAQ to add MAQ values for the reads removed by MAQ filter
#add non-removed but in filter with MAQ MA.filter.nonremove.MAQ
# MA.calls.ref.wTV.MAQ <- MA.calls.ref.wTV
# MA.calls.ref.wTV.MAQ$Category <- 'MAQs'
# MA.calls.ref.wTV.MAQ$Qscore <- MA.calls.ref.wTV.MAQ$MAQ

MA.calls.ref.MAQ <- MA.calls.ref
MA.calls.ref.MAQ$Category <- "Reference calls, MAQs"
MA.calls.ref.MAQ$Qscore <- MA.calls.ref.MAQ$MAQ

#need to test ggplot here by addint 1 at a time in correct order

MA.final <- rbind(MA.calls.ref, MA.calls.nonref.postMAF,MA.noncalls, MA.filter.remove)
o <- order(MA.final$location ) #order by location
MA.final <- MA.final[o,]

MA.final.no.probe <- MA.final[MA.final$CALL12 == "N" & MA.final$Category == "Non-Reference Calls",]
MA.final.no.probe$Category <- "No Genotyping Probe"
## output comment as necessary
# cat("\nMA-filtered data, table")
# table(MA.final$Category)
# table(MA.final$TV)

######
######generate plot data,


MA.final.ggplot <- rbind(MA.calls.ref, MA.calls.ref.MAQ, MA.noncalls, MA.calls.nonref, MA.filter.remove, MA.calls.tv,MA.final.no.probe)
#as backup, this works for P
# MA.final.ggplot <- rbind(MA.calls.ref, MA.calls.ref.wTV.MAQ, MA.calls.nonref.woTV,MA.noncalls, MA.filter.remove,MA.filter.remove.MAQ, MA.filter.nonremove.MAQ,MA.calls.tv)

#generate ggplot based on criteria , comment out for figures
if(length(args)>=2){
g <- ggplot(MA.final.ggplot, aes(x = location, y = Qscore)) +
  geom_point(aes(colour = factor(Category),shape = factor(Category),size = factor(Category)), show.legend=FALSE) +
  scale_color_manual(values=c('#56B4E9','#FF0000','#E69F00','#CCCCCC', '#999999','#00FF00'))+ #blue, red, brown, light grey, dark grey, green
  scale_shape_manual(values=c( 16, 16,16,16,16,16)) +
  scale_size_manual(values=c(2,4,1,1,1,2))
g + scale_fill_discrete(name = "Read type and filters", labels = c("High MAQ filter, not removed, MAQs",
                                                                   "MAQ-filter removed, MAQs",
                                                                   "MAQ-filter removed, Qs",
                                                                   "Non-refernce calls",
                                                                   "True Variants",
                                                                   "MAQs",
                                                                   "Qs"))   +
  # theme(legend.position = c(.1,.15), legend.title = element_blank()) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "transparent", colour = NA), #use the following lines to remove background from plots
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey"),
        axis.line = element_line(colour = "black")

  )+
  scale_x_continuous(name = 'Position', breaks = c(0,5000,10000,15000,20000,25000,30000),limits = c(0,30000))
}else{
### do we want to approximate true variants using genotyping probes? First else approximates, second is if we dont use TVs
#####ggplot takes in input category in alphabetical order. colors have to match it 
  cat("\nUsing suspected TVs, TV positions not provided")
  g <- ggplot(MA.final.ggplot, aes(x = location, y = Qscore)) +
    geom_point(aes(colour = factor(Category),shape = factor(Category),size = factor(Category)), show.legend=FALSE) +
    scale_color_manual(values=c('#56B4E9','#000000','#FF0000','#E69F00','#CCCCCC', '#999999','#00FF00' ))+ #blue, black, red, brown, light grey, dark grey, green
    scale_shape_manual(values=c( 16, 16,16,16,16,16,16)) +
    scale_size_manual(values=c(2,2,4,1,1,1,2))
  g + scale_fill_discrete(name = "Read type and filters", labels = c("High MAQ filter, not removed, MAQs",
                                                                     "MAQ-filter removed, MAQs",
                                                                     "MAQ-filter removed, Qs",
                                                                     "Non-refernce calls",
                                                                     "True Variants",
                                                                     "MAQs",
                                                                     "Qs",
                                                                     "No Genotyping Probe"))   +
    # theme(legend.position = c(.1,.15), legend.title = element_blank()) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "transparent", colour = NA), #use the following lines to remove background from plots
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "grey"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey"),
          axis.line = element_line(colour = "black")
          
    )+
    scale_x_continuous(name = 'Position', breaks = c(0,5000,10000,15000,20000,25000,30000),limits = c(0,30000))
}

# else{
#   g <- ggplot(MA.final.ggplot, aes(x = location, y = Qscore)) +
#     geom_point(aes(colour = factor(Category),shape = factor(Category),size = factor(Category)), show.legend=FALSE) +
#     scale_color_manual(values=c('#56B4E9','#FF0000','#E69F00','#CCCCCC', '#999999'))+ #blue, red, brown, light grey, dark grey, green
#     scale_shape_manual(values=c( 16, 16,16,16,16)) +
#     scale_size_manual(values=c(2,4,1,1,1))
#   g + scale_fill_discrete(name = "Read type and filters", labels = c("High MAQ filter, not removed, MAQs",
#                                                                      "MAQ-filter removed, MAQs",
#                                                                      "MAQ-filter removed, Qs",
#                                                                      "Non-refernce calls",
#                                                                      "True Variants",
#                                                                      "MAQs",
#                                                                      "Qs"))   +
#     # theme(legend.position = c(.1,.15), legend.title = element_blank()) +
#     theme(legend.position = "none",
#           panel.background = element_rect(fill = "transparent", colour = NA), #use the following lines to remove background from plots
#           plot.background = element_rect(fill = "transparent", colour = NA),
#           panel.grid.major = element_line(size = 0.5, linetype = 'solid',
#                                           colour = "grey"),
#           panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "grey"),
#           axis.line = element_line(colour = "black")
#           
#     )+
#     scale_x_continuous(name = 'Position', breaks = c(0,5000,10000,15000,20000,25000,30000),limits = c(0,30000))
# }

exfilename <- strsplit(file.name, "/")[[1]]
png.output.name <- paste("R_Results/Figure_",exfilename[length(exfilename)],"_",date,".png", sep="")
ggsave(png.output.name,
       width = 9.86,#19.72,
       height = 5.57,#11.14,
       units = "in",
       dpi = "retina")

# #Venn Diagram Areas
if(length(args)>=2){
Venn_TV <- sum(ifelse(MA.final$TV == "TV", 1,0))
Venn_Nonref <- sum(ifelse(MA.final$Category == "Non-Reference Calls",1,0)) + sum(ifelse(MA.final$Category == "MAQ filter removed, Qs",1,0))
Venn_filtered <- sum(ifelse(MA.final$Category == "MAQ filter removed, Qs",1,0))

#Venn Diagram combinatiosn for tri.venn
Venn_TV_nonref <- sum(ifelse(MA.final$TV == "TV" & MA.final$FIX > Qth & MA.final$FIXCALL != MA.final$ref, 1,0))
Venn_TV_filtered <- sum(ifelse(MA.final$Category == "MAQ filter removed, Qs" & MA.final$TV == "TV", 1,0))
Venn_nonref_filtered <- sum(ifelse(MA.final$Category == "MAQ filter removed, Qs", 1,0))

Venn_tri_inter <- sum(ifelse(MA.final$TV == "TV" & MA.final$Category == "Non-Reference Calls" & MA.final$MAremove == 1, 1,0))

title <- paste(file.name, "Variant Calls breakdown and filter")

tiff.output.name <- paste("test_tri_",file.name, "_venn_",date,".png", sep="")
tiff(filename = tiff.output.name);
venn.plot <- draw.triple.venn(
  filename = tiff.output.name,
  # output = TRUE,
  resolution = 300,
  height = 3000,
  width = 3000,
  units = "px",
  area1 = Venn_TV, #TVs
  area2 = Venn_Nonref, #variant calls
  area3 = Venn_filtered, #variants removed by MA filter
  n12 = Venn_TV_nonref,
  n23 = Venn_nonref_filtered,
  n13 = Venn_TV_filtered,
  n123 = Venn_tri_inter,
  lwd = c(5,5,5),
  category = c("True Variants", "Variant Calls", "Removed by MA filter"),
  cat.cex = c(1.5,1.5,1.5),
  cat.dist = c(0.015,0.05,.03),
  cat.pos = c(180,145,180),
  cex = c(3,3,3,3,3,3,3),
  col = c('#00FF00','#FF0000','#56B4E9' ),   # added color to coincide with graph, needs testing 2/17/21
  imagetype = "png"
  # main = title
);
grid.draw(venn.plot);
dev.off();
}


################################################################################
if(length(args)>=2){
# (sum(ifelse(MA$Category == "Reference Call",1,0)) + sum(ifelse(MA$TV == "TV",1,0)) )/sum(ifelse(MA$Category != "noncall" & MA$MAremove != 1,1,0))
accuracy.MA.filter <- (sum(ifelse(MA.final$Category == "Reference Calls, Q",1,0)) + sum(ifelse(MA.final$TV == "TV",1,0))) / (dim(MA.final)[1]-dim(MA.noncalls)[1]-dim(MA.filter.remove)[1])
# (sum(ifelse(MA$Category == "Reference Call" | MA$Category == "Non-reference CALL",1,0)))/(sum(ifelse(MA.final$location>0,1,0)))
coverage.MA.filter <- (sum(ifelse(MA.final$Category == "Reference Calls, Q",1,0))+ sum(ifelse(MA.final$Category == "Non-Reference Calls",1,0)))/dim(MA.final)[1]
}else{
  accuracy.MA.filter <- (sum(ifelse(MA.final$Category == "Reference Calls, Q",1,0))) / (dim(MA.final)[1]-dim(MA.noncalls)[1]-dim(MA.filter.remove)[1])
  # (sum(ifelse(MA$Category == "Reference Call" | MA$Category == "Non-reference CALL",1,0)))/(sum(ifelse(MA.final$location>0,1,0)))
  coverage.MA.filter <- (sum(ifelse(MA.final$Category == "Reference Calls, Q",1,0))+ sum(ifelse(MA.final$Category == "Non-Reference Calls",1,0)))/dim(MA.final)[1]
}

######## Variant analysis of WY, USGT5 and US  #comment out unless we want to know specific variants

# #for core0, jse, and fix
# base.WY <- raw[raw$location %in% vars$WY,]  #isolate WY positions from raw
# base.USGT5 <- raw[which(raw$location %in% vars$USA_GA5 ),]
# base.US <- raw[which(raw$location %in% vars$USA_VariantSites ),]
# 
# 
# # core0
# TV.count.core0.wy <- sum(ifelse(base.WY$location %in% TV &
#                                   base.WY$core0 != base.WY$ref & #needs to be non-ref bc these are TVs
#                                   base.WY$Q > Qth,1,0))
# TV.count.core0.USGT5 <- sum(ifelse(base.USGT5$location %in% TV &
#                                      base.USGT5$core0 != base.USGT5$ref & #needs to be non-ref bc these are TVs
#                                      base.USGT5$Q > Qth,1,0))
# TV.count.core0.US <- sum(ifelse(base.US$location %in% TV&
#                                   base.US$core0 != base.US$ref & #needs to be non-ref bc these are TVs
#                                   base.US$Q > Qth,1,0))
# #JSE-fix
# TV.count.JSE.wy <- sum(ifelse(        base.WY$location %in% TV &
#                                         base.WY$FIXCALL != base.WY$ref & #needs to be non-ref bc these are TVs
#                                         base.WY$FIX > Qth,1,0))
# TV.count.JSE.USGT5 <- sum(ifelse(     base.USGT5$location %in% TV &
#                                         base.USGT5$FIXCALL != base.USGT5$ref & #needs to be non-ref bc these are TVs
#                                         base.USGT5$FIX > Qth,1,0))
# TV.count.JSE.US <- sum(ifelse(        base.US$location %in% TV&
#                                         base.US$FIXCALL != base.US$ref & #needs to be non-ref bc these are TVs
#                                         base.US$FIX > Qth,1,0))
# #for new.data
# new.data.WY <- new.data[which(new.data$location %in% vars$WY ),]
# new.data.USGT5 <- new.data[which(new.data$location %in% vars$USA_GA5 ),]
# new.data.US <- new.data[which(new.data$location %in% vars$USA_VariantSites ),]
# 
# 
# #in PCR59_WY59, 23453C is supposed to be a TV, although core0/core12 makes a correct reference call
# TV.count.new.data.wy <- sum(ifelse(   new.data.WY$location %in% TV &
#                                         new.data.WY$FIXCALL != new.data.WY$ref &
#                                         new.data.WY$FIX > Qth,1,0))
# TV.count.new.data.USGT5 <- sum(ifelse(new.data.USGT5$location %in% TV &
#                                         new.data.USGT5$FIXCALL != new.data.USGT5$ref &
#                                         new.data.USGT5$FIX > Qth   ,1,0))
# TV.count.new.data.US <- sum(ifelse(   new.data.US$location %in% TV &
#                                         new.data.US$FIXCALL != new.data.US$ref &
#                                         new.data.US$FIX > Qth,1,0))
# 
# #for MA
# MA.WY <- MA.final[which(MA.final$location %in% vars$WY ),]
# MA.USGT5 <- MA.final[which(MA.final$location %in% vars$USA_GA5 ),]
# MA.US <- MA.final[which(MA.final$location %in% vars$USA_VariantSites ),]
# 
# TV.count.MA.wy <- sum(ifelse(   MA.WY$location %in% TV &
#                                   MA.WY$FIXCALL != MA.WY$ref &
#                                   MA.WY$FIX > Qth,1,0))
# TV.count.MA.USGT5 <- sum(ifelse(MA.USGT5$location %in% TV &
#                                   MA.USGT5$FIXCALL != MA.USGT5$ref &
#                                   MA.USGT5$FIX > Qth   ,1,0))
# TV.count.MA.US <- sum(ifelse(   MA.US$location %in% TV &
#                                   MA.US$FIXCALL != MA.US$ref &
#                                   MA.US$FIX > Qth,1,0))
# 
# #calculate accuarcy for core0, JSE (ie core1/2 fix), and new.data (ie core1/2fix + replacement full)
# acc.core0.Wy <- (sum(ifelse(base.WY$ref == base.WY$core0,1,0))+TV.count.core0.wy)/dim(base.WY)[1]
# acc.core0.USGA5 <- (sum(ifelse(base.USGT5$ref == base.USGT5$core0,1,0))+TV.count.core0.USGT5)/dim(base.USGT5)[1]
# acc.core0.US <- (sum(ifelse(base.US$ref == base.US$core0,1,0))+TV.count.core0.US)/dim(base.US)[1]
# 
# acc.JSE.WY <- (sum(ifelse(base.WY$ref == base.WY$FIXCALL,1,0))+TV.count.JSE.wy)/dim(base.WY)[1]
# acc.JSE.USGA5 <- (sum(ifelse(base.USGT5$ref == base.USGT5$FIXCALL,1,0))+TV.count.JSE.USGT5)/dim(base.USGT5)[1]
# acc.JSE.US <- (sum(ifelse(base.US$ref == base.US$FIXCALL,1,0))+TV.count.JSE.US)/dim(base.US)[1]
# 
# acc.new.data.Wy <- (sum(ifelse(new.data.WY$ref == new.data.WY$FIXCALL,1,0)) + TV.count.new.data.wy)/dim(new.data.WY)[1]
# acc.new.data.USGA5 <- (sum(ifelse(new.data.USGT5$ref == new.data.USGT5$FIXCALL,1,0)) + TV.count.new.data.USGT5)/dim(new.data.USGT5)[1]
# acc.new.data.US <- (sum(ifelse(new.data.US$ref == new.data.US$FIXCALL,1,0)) + TV.count.new.data.US)/dim(new.data.US)[1]
# 
# acc.MA.Wy <- (sum(ifelse(MA.WY$ref == MA.WY$FIXCALL,1,0)) + TV.count.MA.wy)/dim(MA.WY)[1]
# acc.MA.USGA5 <- (sum(ifelse(MA.USGT5$ref == MA.USGT5$FIXCALL,1,0)) + TV.count.MA.USGT5)/dim(MA.USGT5)[1]
# acc.MA.US <- (sum(ifelse(MA.US$ref == MA.US$FIXCALL,1,0)) + TV.count.MA.US)/dim(MA.US)[1]

# rm(base.WY, base.USGT5, base.US,new.data.WY, new.data.USGT5, new.data.US,MA.WY, MA.USGT5, MA.US)

##NEW Variant analysis of B117-1P1,one-percent,tenth-percent  #comment out unless we want to know specific variants

#for core0, jse, and fix
# base.WY <- raw[raw$location %in% vars$WY,]  #isolate WY positions from raw
# base.USGT5 <- raw[which(raw$location %in% vars$USA_GA5 ),]
# base.US <- raw[which(raw$location %in% vars$USA_VariantSites ),]
# 
# 
# # core0
# TV.count.core0.wy <- sum(ifelse(base.WY$location %in% TV &
#                                   base.WY$core0 != base.WY$ref & #needs to be non-ref bc these are TVs
#                                   base.WY$Q > Qth,1,0))
# TV.count.core0.USGT5 <- sum(ifelse(base.USGT5$location %in% TV &
#                                      base.USGT5$core0 != base.USGT5$ref & #needs to be non-ref bc these are TVs
#                                      base.USGT5$Q > Qth,1,0))
# TV.count.core0.US <- sum(ifelse(base.US$location %in% TV&
#                                   base.US$core0 != base.US$ref & #needs to be non-ref bc these are TVs
#                                   base.US$Q > Qth,1,0))
# #JSE-fix
# TV.count.JSE.wy <- sum(ifelse(        base.WY$location %in% TV &
#                                         base.WY$FIXCALL != base.WY$ref & #needs to be non-ref bc these are TVs
#                                         base.WY$FIX > Qth,1,0))
# TV.count.JSE.USGT5 <- sum(ifelse(     base.USGT5$location %in% TV &
#                                         base.USGT5$FIXCALL != base.USGT5$ref & #needs to be non-ref bc these are TVs
#                                         base.USGT5$FIX > Qth,1,0))
# TV.count.JSE.US <- sum(ifelse(        base.US$location %in% TV&
#                                         base.US$FIXCALL != base.US$ref & #needs to be non-ref bc these are TVs
#                                         base.US$FIX > Qth,1,0))
# #for new.data
# new.data.WY <- new.data[which(new.data$location %in% vars$WY ),]
# new.data.USGT5 <- new.data[which(new.data$location %in% vars$USA_GA5 ),]
# new.data.US <- new.data[which(new.data$location %in% vars$USA_VariantSites ),]
# 
# 
# #in PCR59_WY59, 23453C is supposed to be a TV, although core0/core12 makes a correct reference call
# TV.count.new.data.wy <- sum(ifelse(   new.data.WY$location %in% TV &
#                                         new.data.WY$FIXCALL != new.data.WY$ref &
#                                         new.data.WY$FIX > Qth,1,0))
# TV.count.new.data.USGT5 <- sum(ifelse(new.data.USGT5$location %in% TV &
#                                         new.data.USGT5$FIXCALL != new.data.USGT5$ref &
#                                         new.data.USGT5$FIX > Qth   ,1,0))
# TV.count.new.data.US <- sum(ifelse(   new.data.US$location %in% TV &
#                                         new.data.US$FIXCALL != new.data.US$ref &
#                                         new.data.US$FIX > Qth,1,0))

#for MA
# table(MA.strainVar$Category)
# exfilename <- strsplit(file.name, "/")[[1]]
# output.file <- paste(exfilename[length(exfilename)])

MA.strainVar <- MA.final[MA.final$location %in% vars$all,]
MA.one <- MA.final[MA.final$location %in% vars$one ,]
MA.tenth <- MA.final[MA.final$location %in% vars$tenth ,]

VarCount.B117 <- sum(ifelse(  MA.strainVar == "Non-Reference Calls",1,0))
VarCount.one <- sum(ifelse(MA.one == "Non-Reference Calls",1,0))
VarCount.tenth <- sum(ifelse(MA.tenth == "Non-Reference Calls",1,0))

MA.strainVar <- MA.strainVar[MA.strainVar$Category == "Non-Reference Calls" ,]
MA.one <- MA.one[MA.one$Category == "Non-Reference Calls" ,]
MA.tenth <- MA.tenth[MA.tenth$Category == "Non-Reference Calls" ,]

#calculate accuarcy for core0, JSE (ie core1/2 fix), and new.data (ie core1/2fix + replacement full)
# acc.core0.Wy <- (sum(ifelse(base.WY$ref == base.WY$core0,1,0))+TV.count.core0.wy)/dim(base.WY)[1]
# acc.core0.USGA5 <- (sum(ifelse(base.USGT5$ref == base.USGT5$core0,1,0))+TV.count.core0.USGT5)/dim(base.USGT5)[1]
# acc.core0.US <- (sum(ifelse(base.US$ref == base.US$core0,1,0))+TV.count.core0.US)/dim(base.US)[1]
# 
# acc.JSE.WY <- (sum(ifelse(base.WY$ref == base.WY$FIXCALL,1,0))+TV.count.JSE.wy)/dim(base.WY)[1]
# acc.JSE.USGA5 <- (sum(ifelse(base.USGT5$ref == base.USGT5$FIXCALL,1,0))+TV.count.JSE.USGT5)/dim(base.USGT5)[1]
# acc.JSE.US <- (sum(ifelse(base.US$ref == base.US$FIXCALL,1,0))+TV.count.JSE.US)/dim(base.US)[1]
# 
# acc.new.data.Wy <- (sum(ifelse(new.data.WY$ref == new.data.WY$FIXCALL,1,0)) + TV.count.new.data.wy)/dim(new.data.WY)[1]
# acc.new.data.USGA5 <- (sum(ifelse(new.data.USGT5$ref == new.data.USGT5$FIXCALL,1,0)) + TV.count.new.data.USGT5)/dim(new.data.USGT5)[1]
# acc.new.data.US <- (sum(ifelse(new.data.US$ref == new.data.US$FIXCALL,1,0)) + TV.count.new.data.US)/dim(new.data.US)[1]
# 
 # acc.MA.Wy <- (sum(ifelse(MA.WY$ref == MA.WY$FIXCALL,1,0)) + TV.count.MA.wy)/dim(MA.WY)[1]
# acc.MA.USGA5 <- (sum(ifelse(MA.USGT5$ref == MA.USGT5$FIXCALL,1,0)) + TV.count.MA.USGT5)/dim(MA.USGT5)[1]
# acc.MA.US <- (sum(ifelse(MA.US$ref == MA.US$FIXCALL,1,0)) + TV.count.MA.US)/dim(MA.US)[1]
# 
# rm(base.WY, base.USGT5, base.US,new.data.WY, new.data.USGT5, new.data.US,MA.WY, MA.USGT5, MA.US)

#######################output text#######################
#comment in out what is wanted
# MA.final <- rbind(MA.calls.ref, MA.calls.nonref.postMAF,MA.noncalls, MA.filter.remove)

# sum(ifelse(MA.final$FIX > Qth,1,0)) #calls before filter
# sum(ifelse(MA.final$FIX > Qth & MA.final$Category != "MAQ-filter removed, Qs",1,0)) #calls after filter
# sum(ifelse(MA.final$FIX < Qth,1,0)) #non-calls, before filter
# sum(ifelse(MA.final$FIX < Qth | MA.final$Category == "MAQ-filter removed, Qs",1,0)) #non-calls, before filter
# sum(ifelse(MA.final$FIXCALL != MA.final$ref & MA.final$FIX > Qth & MA.final$FIX < 30,1,0)) #into filter
# sum(ifelse(MA.final$FIXCALL != MA.final$ref & MA.final$FIX > Qth & MA.final$FIX < 30 & #into filter and filtered out
#              MA.final$MAQ < mean(MA.final$MAQ)-2*sd(MA.final$MAQ),1,0)) 
# sum(ifelse(MA.final$FIXCALL != MA.final$ref & MA.final$FIX > Qth & MA.final$FIX < 30 & #into filter and not filtered
#              MA.final$MAQ > mean(MA.final$MAQ)-2*sd(MA.final$MAQ),1,0)) 
# 
# sum(ifelse(MA.final$Category == "Reference Calls, Q",1,0))
# sum(ifelse(MA.final$Category == "Non-Reference Calls",1,0))
# sum(ifelse(MA.final$Category == "Noncalls",1,0))
# sum(ifelse(MA.final$Category == "MAQ-filter removed, Qs",1,0))

# unique(MA.final$Category)
# 
# cat("\nMA post filter Accuracy\t",accuracy.MA.filter)
# cat("\nMA post filter reads\t",dim(MA.final)[1])
# cat("\nMA post filter total calls\t",dim(MA.final)[1]-dim(MA.filter.remove)[1]-dim(MA.noncalls)[1]-dim(MA.filter.remove.MAQ)[1]-dim(MA.filter.nonremove.MAQ)[1])
# cat("\nMA post filter reference calls\t",dim(MA.calls.ref)[1])
# cat("\nNumber variant calls removed by MA filter\t", dim(MA.filter.remove.MAQ)[1])
# cat("\nNumber variant calls remaining\t", dim(MA.calls.nonref.postMAF)[1])
# # 
# cat("\nNumber of TV found in WY positions, core0\t", TV.count.core0.wy)
# cat("\nNumber of TV found in USGTA5 positions, core0\t", TV.count.core0.USGT5)
# cat("\nNumber of TV found in US positions, core0\t", TV.count.core0.US)
# 
# cat("\nNumber of TV found in WY positions, JSEfix\t", TV.count.JSE.wy)
# cat("\nNumber of TV found in USGTA5 positions, JSEfix\t", TV.count.JSE.USGT5)
# cat("\nNumber of TV found in US positions, JSEfix\t", TV.count.JSE.US)
# 
# cat("\nNumber of TV found in WY positions, fix\t", TV.count.new.data.wy)
# cat("\nNumber of TV found in USGTA5 positions, fix\t", TV.count.new.data.USGT5)
# cat("\nNumber of TV found in US positions, fix\t", TV.count.new.data.US)
# 
# cat("\nNumber of TV found in WY positions, MA\t", TV.count.MA.wy)
# cat("\nNumber of TV found in USGTA5 positions, MA\t", TV.count.MA.USGT5)
# cat("\nNumber of TV found in US positions, MA\t", TV.count.MA.US)



exfilename <- strsplit(file.name, "/")[[1]]
output.file <- paste(exfilename[length(exfilename)])
result.name <- paste("R_Results/",output.file,"_RESULTS.txt",sep="")
cat("RESULTS:\n", "Sample Name:",file.name, file=result.name, append = FALSE)

accuracy.raw <- sum(ifelse(raw$Category == "Reference Call",1,0)) / sum(ifelse(raw$Q > Qth, 1, 0))
accuracy.jse <- sum(ifelse(raw$Category2 == "Reference Call",1,0)) / sum(ifelse(raw$FIX > Qth, 1, 0))

coverage.raw <- sum(ifelse(raw$Q>Qth,1,0))/sum(ifelse(raw$location>1,1,0))
coverage.raw.jse <- sum(ifelse(raw$FIX>Qth,1,0))/sum(ifelse(raw$location>1,1,0))
coverage.new.data <- sum(ifelse(new.data$FIX>Qth,1,0))/sum(ifelse(new.data$location>1,1,0))

if(length(args)>=2){
  accuracy.jse.fix.tv <- (sum(ifelse(new.data$Category == "Reference Call",1,0)) + sum(ifelse(new.data$TV == "TV" & new.data$FIX > Qth, 1,0))) / sum(ifelse(new.data$FIX > Qth, 1, 0))
}else{
  accuracy.jse.fix.tv <- (sum(ifelse(new.data$Category == "Reference Call",1,0)) ) / sum(ifelse(new.data$FIX > Qth, 1, 0))
}

## output comment as necessary
# cat("\nCore0 categorical breakdown")
# table(raw$Category)
# 
# cat("\nJSE filter categorical breakdown")
# table(raw$Category2)

# sum(ifelse(raw$FIXCALL == raw$ref & raw$FIX > Qth,1,0))


# cat("\nJSE filter + additional fix + TV categorical breakdown")
# table(new.data$Category)
# table(new.data$TV)

# sum(ifelse(new.data$FIXCALL == new.data$ref & new.data$FIX > Qth,1,0))

# cat("\nMA filter + TV categorical breakdown")
# table(MA.final$Category)
# table(MA.final$TV)

# as.data.table(temp)
# print(temp, colnames = FALSE)

# cat("\n Number of TVs removed by MA-filter\t", sum(ifelse(MA.final$Category == "MAQ filter removed, Qs" & MA.final$TV == "TV",1,0)))

# cat("\n\tNon-calls\tCalls\tVariant Calls\tReference Calls\tFiltered Calls\tMA filter removed\tTrue Variants\tAccuracy\tcoverage")
# cat("\ncore0 \t", raw.noncall.reads,
#     "\t",raw.call.reads,
#     "\t",raw.nonref.calls,
#     "\t",raw.ref.call,
#     "\t",0,
#     "\t",0,
#     "\t",0,
#     "\t",accuracy.raw*100,
#     "\t",coverage.raw*100,
#     "\t")
# 
# cat("\nJSE-fix\t", sum(ifelse(raw$Category2 == "noncall",1,0)),
#     "\t",raw.jse.call.reads,
#     "\t",raw.jse.potential.variant.calls,
#     "\t",raw.jse.ref.calls,
#     "\t",0,
#     "\t",0,
#     "\t",0,
#     "\t",accuracy.jse*100,
#     "\t",coverage.raw.jse*100,
#     "\t")
# 
# cat("\n2nd replacement with True variants from list\t", sum(ifelse(new.data$Category == "noncall",1,0)),
#     "\t",sum(ifelse(new.data$FIX>Qth,1,0)),
#     "\t",sum(ifelse(new.data$Category == "Non-reference Call",1,0)),
#     "\t",sum(ifelse(new.data$Category == "Reference Call",1,0)),
#     "\t",0,
#     "\t",0,  
#     "\t",sum(ifelse(new.data$TV == "TV",1,0)),
#     "\t",accuracy.jse.fix.tv*100,
#     "\t",coverage.new.data*100,
#     "\t")
# 
# cat("\nJSE-fix with MA filter and true variants\t",sum(ifelse(MA.final$Category == "noncall" & MA.final$MAremove == "0" ,1,0)) + sum(ifelse(MA.final$Category == "MAQ filter removed, Qs",1,0)),
#     "\t",sum(ifelse(MA.final$FIX > Qth,1,0)),
#     "\t",sum(ifelse(MA.final$FIX > Qth & MA.final$FIXCALL != MA.final$ref & MA.final$filter == 0,1,0)),
#     "\t",sum(ifelse(MA.final$Category == "Reference Calls, Q",1,0)),
#     "\t",sum(ifelse(MA.final$filter == "filtered",1,0)),
#     "\t",sum(ifelse(MA.final$MAremove == 1,1,0)),
#     "\t",sum(ifelse(MA.final$TV == "TV",1,0)),
#     "\t",accuracy.MA.filter*100,
#     "\t",coverage.MA.filter*100,
#     "\t")

# cat("\n Core0 WY/USGT5/US accuracy\t",
#     acc.core0.Wy ,"\t",
#     acc.core0.USGA5 ,"\t",
#     acc.core0.US 
# )
# 
# cat("\nReplacement 1 WY/USGT5/US accuracy\t",
#     acc.JSE.WY ,"\t",
#     acc.JSE.USGA5 ,"\t",
#     acc.JSE.US 
# )
# 
# cat("\nReplacement 2 WY/USGT5/US accuracy\t",
#     acc.new.data.Wy ,"\t",
#     acc.new.data.USGA5 ,"\t",
#     acc.new.data.US 
# )
# 
# cat("\nMA-Filter WY/USGT5/US accuracy\t",
#     acc.MA.Wy ,"\t",
#     acc.MA.USGA5 ,"\t",
#     acc.MA.US
# )

cat("\n\tNon-calls\tCalls\tVariant Calls\tReference Calls\tTrue Variants\tAccuracy\tcoverage",
    file=result.name, 
    append = TRUE)
cat("\ncore0 \t", raw.noncall.reads,
    "\t",raw.call.reads,
    "\t",raw.nonref.calls,
    "\t",raw.ref.call,
    "\t",0,
    "\t",accuracy.raw*100,
    "\t",coverage.raw*100,
    "\t",
    file=result.name, 
    append = TRUE)

# cat("\nJSE-fix\t", sum(ifelse(raw$Category2 == "noncall",1,0)),
#     "\t",raw.jse.call.reads,
#     "\t",raw.jse.potential.variant.calls,
#     "\t",raw.jse.ref.calls,
#     "\t",0,
#     "\t",accuracy.jse*100,
#     "\t",coverage.raw.jse*100,
#     "\t")
if(length(args)>1){
cat("\n2nd replacement with True variants from list\t", sum(ifelse(new.data$Category == "noncall",1,0)),
    "\t",sum(ifelse(new.data$FIX>Qth,1,0)),
    "\t",sum(ifelse(new.data$Category == "Non-reference Call",1,0)),
    "\t",sum(ifelse(new.data$Category == "Reference Call",1,0)),
    "\t",sum(ifelse(new.data$TV == "TV",1,0)),
    "\t",accuracy.jse.fix.tv*100,
    "\t",coverage.new.data*100,
    "\t",
    file=result.name, 
    append = TRUE)
}else{
  cat("\nUsing suspected TVs, TV positions not provided")
  cat("\n2nd replacement with True variants from list\t", sum(ifelse(new.data$Category == "noncall",1,0)),
      "\t",sum(ifelse(new.data$FIX>Qth,1,0)),
      "\t",sum(ifelse(new.data$Category == "Non-reference Call",1,0)),
      "\t",sum(ifelse(new.data$Category == "Reference Call",1,0)),
      "\t",dim(TV_JSE)[1],
      "\t",accuracy.jse.fix.tv*100,
      "\t",coverage.new.data*100,
      "\t",
      file=result.name, 
      append = TRUE)
  }
if(length(args)>1){
cat("\nJSE-fix with MA filter and true variants\t",sum(ifelse(MA.final$Category == "noncall" & MA.final$MAremove == "0" ,1,0)) + sum(ifelse(MA.final$Category == "MAQ filter removed, Qs",1,0)),
    "\t",sum(ifelse(MA.final$FIX > Qth,1,0)),
    "\t",sum(ifelse(MA.final$FIX > Qth & MA.final$FIXCALL != MA.final$ref & MA.final$MAremove != 1,1,0)),
    "\t",sum(ifelse(MA.final$Category == "Reference Calls, Q",1,0)),
    "\t",sum(ifelse(MA.final$TV == "TV",1,0)),
    "\t",accuracy.MA.filter*100,
    "\t",coverage.MA.filter*100,
    "\t",
    file=result.name, 
    append = TRUE)
}else{
  cat("\nUsing suspected TVs, TV positions not provided")
  cat("\nJSE-fix with MA filter and true variants\t",sum(ifelse(MA.final$Category == "noncall" & MA.final$MAremove == "0" ,1,0)) + sum(ifelse(MA.final$Category == "MAQ filter removed, Qs",1,0)),
      "\t",sum(ifelse(MA.final$FIX > Qth,1,0)),
      "\t",sum(ifelse(MA.final$FIX > Qth & MA.final$FIXCALL != MA.final$ref & MA.final$MAremove != 1,1,0)),
      "\t",sum(ifelse(MA.final$Category == "Reference Calls, Q",1,0)),
      "\t",dim(TV_JSE)[1],
      "\t",accuracy.MA.filter*100,
      "\t",coverage.MA.filter*100,
      "\t",
      file=result.name, 
      append = TRUE)
}
if(length(args)>=2){
  cat("\nTVs from premade list \n")
  print(TV_JSETV_JSE[,1:9], quote = FALSE, row.names = FALSE)
}
cat("\nPotential genotyping probe-set variant calls \n",
    file=result.name, 
    append = TRUE)
sink(file = result.name,append = TRUE)
print(TV_JSE[,1:9], quote = FALSE, row.names = FALSE, 
      file=result.name, 
      append = TRUE)
cat("\nNumber of variant calls at B.1.1.7.p1 locations:",VarCount.B117 )
cat("\nNumber of variant calls at >1% frequency locations:",VarCount.one )
cat("\nNumber of variant calls at >0.1%frequency locations:",VarCount.tenth )



VarCount.B117 <- sum(ifelse(  MA.strainVar == "Non-Reference Calls",1,0))
VarCount.one <- sum(ifelse(MA.one == "Non-Reference Calls",1,0))
VarCount.tenth <- sum(ifelse(MA.tenth == "Non-Reference Calls",1,0))

temp <- vars[vars$all %in% MA.strainVar$location,]
cat("\nStrain specific mutation breakdown")
table(temp$strain)

cat("\nAll Variant Calls \n")
print(MA.final[MA.final$Category == "Non-Reference Calls" ,c(1:7,10)],row.names = FALSE)

cat("\nVariant Calls, core0=core12 \n")
print(MA.final[MA.final$Category == "Non-Reference Calls" & MA.final$core0 == MA.final$CALL12,c(1:7,10)],row.names = FALSE)
#isolate non-refernce calls, use these calls only for fasta
temp <- MA.final[MA.final$Category == "Non-Reference Calls" & MA.final$core0 == MA.final$CALL12,]
temp$fasta <- temp$FIXCALL
temp.all <- MA.final[!(MA.final$location %in% temp$location),]
temp.all$fasta <- temp.all$ref
temp <- rbind(temp, temp.all)
#convert call noncalls and filter removed to Ns
temp1 <- temp[temp$Category == "noncall" |temp$Category == "MAQ filter removed, Qs",]
temp1$fasta <- "N"
temp <- temp[!(temp$location %in% temp1$location),]
temp <- rbind(temp, temp1)
#order by position after adding above two changes in $fasta column
o <- order(temp$location) #order by location
temp <- temp[o,] 
temp <- temp$fasta
temp.name <- paste("R_Results/",output.file, "_core0=core12", collapse="",sep="")
export.fasta(temp, temp.name)
# print(TV_JSE2, quote = FALSE, row.names = FALSE)

cat("\nCalls removed by MA-filter\n")
print(MA.final[MA.final$MAremove == 1,c(1:7,10)],row.names = FALSE)
#isolate Calls removed by MA-filters, use these calls only for fasta
temp <- MA.final[MA.final$MAremove == 1,]
temp$fasta <- temp$FIXCALL
temp.all <- MA.final[!(MA.final$location %in% temp$location),]
temp.all$fasta <- temp.all$ref
temp <- rbind(temp, temp.all)
#convert call noncalls and filter removed to Ns
temp1 <- temp[temp$Category == "noncall" |temp$Category == "MAQ filter removed, Qs",]
temp1$fasta <- "N"
temp <- temp[!(temp$location %in% temp1$location),]
temp <- rbind(temp, temp1)
#order by position after adding above two changes in $fasta column
o <- order(temp$location) #order by location
temp <- temp[o,] 
temp <- temp$fasta
temp.name <- paste("R_Results/",output.file, "_removed_by_MA-filter", collapse="",sep="")
export.fasta(temp, temp.name)

# print(TV_JSE2, quote = FALSE, row.names = FALSE)

cat("\n Strain specific Variants from list\n")
print(MA.strainVar[,c(1:8,11,13)],row.names = FALSE)
#isolate  Strain specific Variants, use these calls only for fasta
temp <- MA.strainVar
temp$fasta <- temp$FIXCALL
temp.all <- MA.final[!(MA.final$location %in% temp$location),]
temp.all$fasta <- temp.all$ref
temp <- rbind(temp, temp.all)
#convert call noncalls and filter removed to Ns
temp1 <- temp[temp$Category == "noncall" |temp$Category == "MAQ filter removed, Qs",]
temp1$fasta <- "N"
temp <- temp[!(temp$location %in% temp1$location),]
temp <- rbind(temp, temp1)
#order by position after adding above two changes in $fasta column
o <- order(temp$location) #order by location
temp <- temp[o,] 
temp <- temp$fasta
temp.name <- paste("R_Results/",output.file, "_Strain_specific", collapse="",sep="")
export.fasta(temp, temp.name)

cat("\n > 1% mutation frequency from list\n")
print(MA.one[,c(1:8,11,13)],row.names = FALSE)
#isolate 1%, use these calls only for fasta
temp <- MA.one
temp$fasta <- temp$FIXCALL
temp.all <- MA.final[!(MA.final$location %in% temp$location),]
temp.all$fasta <- temp.all$ref
temp <- rbind(temp, temp.all)
#convert call noncalls and filter removed to Ns
temp1 <- temp[temp$Category == "noncall" |temp$Category == "MAQ filter removed, Qs",]
temp1$fasta <- "N"
temp <- temp[!(temp$location %in% temp1$location),]
temp <- rbind(temp, temp1)
#order by position after adding above two changes in $fasta column
o <- order(temp$location) #order by location
temp <- temp[o,] 
temp <- temp$fasta
temp.name <- paste("R_Results/",output.file, "_oneP", collapse="",sep="")
export.fasta(temp, temp.name)

cat("\n > 0.1% mutation frequency from list\n")
print(MA.tenth[,c(1:8,11,13)],row.names = FALSE)
#isolate 0.1% , use these calls only for fasta
temp <- MA.tenth
temp$fasta <- temp$FIXCALL
temp.all <- MA.final[!(MA.final$location %in% temp$location),]
temp.all$fasta <- temp.all$ref
temp <- rbind(temp, temp.all)
#convert call noncalls and filter removed to Ns
temp1 <- temp[temp$Category == "noncall" |temp$Category == "MAQ filter removed, Qs",]
temp1$fasta <- "N"
temp <- temp[!(temp$location %in% temp1$location),]
temp <- rbind(temp, temp1)
#order by position after adding above two changes in $fasta column
o <- order(temp$location) #order by location
temp <- temp[o,] 
temp <- temp$fasta
temp.name <- paste("R_Results/",output.file, "_tenthP", collapse="",sep="")
export.fasta(temp, temp.name)


sink(file=NULL)



exfilename <- strsplit(file.name, "/")[[1]]
output.file <- paste("R_Results/",exfilename[length(exfilename)], sep="")


export.temp <- MA.final[,c(1,2,7,8)]
export.data(export.temp, file.name.handle = output.file, type = "ALL", notes = "All Calls") #testing

export.temp <- MA.final[MA.final$Category == "Non-Reference Calls",c(1,2,7,8)]
export.data(export.temp, file.name.handle = output.file, type = "VarCall", notes = "Variant Calls") #testing

export.temp <- MA.final[MA.final$Category == "Non-Reference Calls" & MA.final$core0 == MA.final$CALL12,c(1,2,7,8)]
export.data(export.temp, file.name.handle = output.file, type = "HQVarCall", notes = "HQVariant Calls, core1=core12") #testing

export.temp <- MA.final[MA.final$Category == "Non-Reference Calls" & MA.final$core0 != MA.final$CALL12,c(1:8)]
export.data(export.temp, file.name.handle = output.file, type = "LQVarCall", notes = "LQVariant Calls, core1!=core12") #testing

export.temp <- MA.final[MA.final$Category == "noncall" | MA.final$Category == "MAQ filter removed, Qs",]
export.temp$FIXCALL <- "N"

# high quality variants where core 0 1/2 agree 

export.temp1 <- MA.final[!(MA.final$location %in% export.temp$location),]
export.temp1 <- rbind(export.temp1, export.temp)
o <- order(export.temp1$location) #order by location
export.temp1 <- export.temp1[o,] 
export.temp1 <- export.temp1$FIXCALL

export.fasta(export.temp1,output.file)
}


# input.argument <- "PCR64_WY64.csv" #by default
input.argument <- "all_calls_11-Mar-2021_012_R_19_Clade.csv" #by default
# input.argument <- "file.names.txt" #second test 
input.argument <- commandArgs(trailingOnly = TRUE)


if(str_sub(input.argument, -4) == ".txt"){
  
  i=1
  list <- read.table(input.argument)
  
  while(i< dim(list)[1]-3){
    sequence.sample(list[i,])
    i<-i+1
    cat("\nCompleted sample", list[i,],"\n\n")
  }
  cat("\nFinished generating results\n\n")
  unlink("file.names.txt")
}else{
  sequence.sample(input.argument)
}

# # table(MA.final$Category)
# 
# # var_table <- read.table(file = "variant_table.csv", sep=",", header = TRUE)
# # var_table <- var_table[21:62,]
# var_table <- var_table[1:59,c(7,9)]
# var_table$variant.1 <- str_sub(temp$V9, start = 2L)
# 
# temp <- var_table[c(1:12,16:59),]
# temp$variant.1 <- str_sub(temp$variant.1, start = 2L,end = str_length(temp$variant.1)-1)
# 
# var_table <- var_table[13:15,]
# var_table <- rbind(var_table, temp)
# # temp$variant.1 <- str_pad(temp$variant.1, width = 6, pad = "A", side ="left")
# tempV10<- str_c(temp$V9, sep=" ")
# # var_table <- var_table[c(1:29, 3:123),]
# # var_table <- var_table[c(1:3, 7:62),]
# 
# export.data(var_table, file.name = "var_list")

