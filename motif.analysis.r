

### PWM ####
#require(BiocManager)
#BiocManager::install("seqLogo")
require(seqLogo)

mFile <- system.file("extdata/pwm1", package="seqLogo")
m <- read.table(mFile)

p <- makePWM(m)

seqLogo(p)



### motif discovery ###

#require(BiocManager)
#BiocManager::install("rGADEM")
#BiocManager::install("rtracklayer")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
require(rGADEM)


path <- system.file("extdata/Test_100.fasta",package="rGADEM")
FastaFile<-paste("",path,sep="")
Sequences <- readDNAStringSet(FastaFile, "fasta")
gadem<-GADEM(Sequences,verbose=1, genome=Hsapiens) ### motif discovery from given sequences

consensus(gadem) ### consensus sequence
pwm = getPWM(gadem) ### PWM matrix

seqLogo(pwm$GTGTTTACATG) ### visualization of PWM matrix


matchPWM(as.matrix(pwm$GTGTTTACATG),Sequences[[1]]) ### motif search
matchPWM(as.matrix(pwm$GTGTTTACATG),Sequences[[2]]) ### motif search
matchPWM(as.matrix(pwm$GTGTTTACATG),Sequences[[3]]) ### motif search
