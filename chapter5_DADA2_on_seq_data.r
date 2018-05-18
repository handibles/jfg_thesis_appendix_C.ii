## ZSH
# cd /media/cart/cart_sandbox/MV/mv_reads
# cp /media/cart/cart_sandbox/z.raw_sequences/2017_02_15_JamieFitzGerald_MVOELKLEIN_un-trimmed_data.zip /media/cart/cart_sandbox/MV/
# gzip -drf *.zip
# rm 2017_02_15_JamieFitzGerald_MVOELKLEIN_un-trimmed_data.zip

#setwd("~/jf/mv_r/")

## Q UE R I E S 
#
# why not merge and then do error estimates? Twice the work for incomplete infos?
# sequence variant sequences?
# what data do we want of the sequence processing?
#  where and when data is lost..
#  


#                       < ! >                   #
##   H E A V Y   C O P Y P A S T A   Z O N E   ##
#                     do not shame              #


##   P R I M E 
# from howto parallel
library('parallel')
no_cores <- detectCores() - 1 
cl <- makeCluster(no_cores)
# stopCluster(cl)
library('dada2')


##   P R E P
# keep path in same partition
path <-'~/jf/mv_reads'
list.files(path)

## make an ordered vector of the filenames, instaed of calling the list up repeatedly.
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1.fastq"))
fnRs <- sort(list.files(path, pattern="_R2.fastq"))

## clever
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)


# # plot all qscores, only run once
# pdf()
# # would be better to do all the analyses in one go
# # and then index into them. 
# plotQualityProfile(fnFs[1:9])
# plotQualityProfile(fnFs[10:18])
# plotQualityProfile(fnFs[19:27])
# plotQualityProfile(fnFs[28:33])
# plotQualityProfile(fnRs[1:9])
# plotQualityProfile(fnRs[10:18])
# plotQualityProfile(fnRs[19:27])
# plotQualityProfile(fnRs[28:33])
# dev.off()


## Q U A L I T Y   F I L T E R I N G   &   T R I M M I N G
## decide what sort of trimming you can do! 
# Product should be 410 bp, need a 430 bp overlap
# problematic, but try 241-201 - gives overlap and Q>25.
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory with original reads
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# FOOL! Need to trim off the primers at fore (19bp) and aft (20bp) also. 
# added trimLeft=c(19,20).
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(241,201), trimLeft=c(19,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)


##  E R R O R - M O D E L   E S T I M A T I O N
# estimate the error model for the data from 1e6 reads
# through progressive leave 'n' out approach. jfg added randomize.
# multithread looks a good call
errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE)
save.image("~/jf/mv_r/.RData",compress=FALSE)

# check convergence levels differ by orders of magnitude 
dada2:::checkConvergence(errF)
dada2:::checkConvergence(errR)
    # plot error substitution rates
    pdf('mv_errorF-R.pdf')
    plotErrors(errF, nominalQ=TRUE)
    plotErrors(errR, nominalQ=TRUE)
    dev.off()
    # can give an error message: "Self-consistency loop terminated before convergence."
    # see https://benjjneb.github.io/dada2/faq.html ; tl:dr if plots look reasonable, not a problem. 
    # Can allow larger estiamtes, but MAX_CONSIST seems to be a retired option
    # errF.20 <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE, MAX_CONSIST=20)
    # errR.20 <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE, MAX_CONSIST=20)
    


##  D E R E P L I C A T I O N   (seq.var picking)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# save on the system partition: otherwise hangs in a timely fashion
save.image("~/jf/mv_r/.RData",compress=FALSE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##  ' S A M P L E   I N F E R E N C E '  -  combine EE model and SV's
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
# save on the system partition: otherwise hangs in a timely fashion
save.image("~/jf/mv_r/.RData",compress=FALSE)
head(dadaFs[[1]])


##   M E R G E   P A I R S
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


##   S E Q T A B L E 
seqtab <- makeSequenceTable(mergers)
## The sequences being tabled vary in length.
dim(seqtab)
# Inspect *distribution* of sequence lengths
table(nchar(getSequences(seqtab)))


##   N O   C H I M
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# [1]   33 4167
sum(seqtab.nochim)/sum(seqtab)
# [1] 0.7733314
## 23% of reads removed from the run - much displeasure at your announce.
## Tut loses: 4%. Advise "upstream processing may need to be revisited" 



##   R E A D S   I N  /  R E A D S   O U T 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
View(track)
write.table(track,'mv_seqprogression.txt',sep='\t')

##   T A X O N O M Y
# add  tryRC=TRUE :: retry tax assignments so can iron out euk seqs and close issue: https://github.com/benjjneb/dada2/issues/366
taxa_rc <- assignTaxonomy(seqtab.nochim, "silva_nr_v128_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
# unname(head(taxa_rc))
taxa_rc_sp <- addSpecies(taxa_rc, "silva_species_assignment_v128.fa.gz")

##   A C C U R A C Y   (workflow contains mock community, possibly not appropriate, but could add to our set!)
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
## "DADA2 inferred 20 sample sequences present in the Mock community."
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")


## D U P L I C A T E   T H E   O U T P U T
taxa.orig <- taxa_rc
seqtab.nochim.orig <- seqtab.nochim



## S A V I N G
# saveRDS(file='.RData.RDS',compress=FALSE)

save.image("~/jf/mv_r/.RData",compress=FALSE)


