###############################
######## Assignment 7 #########
###############################

*------------------------------------
# Measuring pairwise dissimilarities
*------------------------------------

library (TraMineRextras)

#1. Create a state sequence object with explicit short and long labels 
    #from the biofam data.
data(biofam)
mycol<-brewer.pal(,"RdBu")
biofam.lab<-c("Parent", "Left", "Married","Left+Marr","Child","Left+Child","Left+Marr+Child","Divorced")
biofam.shortlab<-c("P", "L", "M", "LM", "C","LC","LMC", "D")
biofamseq <- seqdef(biofam[,10:25],states=biofam.shortlab,labels=biofam.lab)
summary (biofamseq)

#2. Compute the matrix of pairwise HAM distances between sequences and 
    #display the results for the first 5 sequences.
distHAM <- seqdist(biofamseq, method = "HAM")
head(distHam)

#3. Plot the first 2 sequences and check that the HAM distance is the 
    #number of non matching positions between them.
plot(biofam.seq[1:2,])

#4. Check on the biofam data that the LCS distance provides the same 
    #non-normalized distances as OM with indel=1 and a constant substitution 
    #cost of 2 which we denote as OM(1,2).
distLCS<- seqdist(biofamseq, method = "LCS")
distOM <- seqdist(biofamseq, method = "HAM", sm="CONSTANT")
summary (distLCS!=distOM)

#5. Compute the normalized LCS and OM distances by setting norm=TRUE and verify 
    #that the normalized OM(1,2) distance is twice the normalized LCS. How do you 
    #explain that?
distLCSnorm<- seqdist(biofamseq, method = "LCS", norm=TRUE)
distOMnorm <- seqdist(biofamseq, method = "HAM", sm="CONSTANT", norm=TRUE)
max(distLCSnorm)
max(distOMnorm)
sum((2*distLCSnorm)!=distOMnorm)
#6. Check that normalizing OM(1,2) distances with norm="gmean" you get the 
    #normal- ized LCS distance. Conversely check that normalizing LCS with 
    #norm="maxlength" you get the normalized OM(1,2) distance.

distLCSnorm2<- seqdist(biofamseq, method = "LCS", norm="maxlength" )
distOMnorm2<- seqdist(biofamseq, method = "HAM", sm="CONSTANT", norm="gmean")
max(distLCSnorm2)
max(distOMnorm2)
sum(distLCSnorm!=distOMnorm2)
sum(distLCSnorm2!=distOMnorm)