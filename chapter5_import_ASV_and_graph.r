#! /usr/bin/Rscript


library('heatmaply')
library('ape')
library('phyloseq')
library('ggplot2')
library('DESeq2')
library('scales')


# =====================================================================================

#				f r o m   D A D A 2 . . . 

# =====================================================================================





##   R E S O L V E   N A M E S :
##   See also separate .R file (mv_renames.r) for more explanation / validation
identical(rownames(taxa.orig), colnames(seqtab.nochim.orig))   # check TRUE
uniqueseqs_rc <- as.data.frame(as.character(colnames(seqtab.nochim)))  ; colnames(uniqueseqs_rc) <-c('UniqueSeq')
svnames <- paste0('mvsv',formatC( 1:ncol(seqtab.nochim), width=4,format='d',flag='0')) # make generic 'mvsv###' ids for uniqueseqs
rownames(uniqueseqs_rc) <- svnames
colnames(seqtab.nochim) <- svnames
rownames(taxa_rc_sp) <- svnames          # note 'sp' designation
taxa_rc_sp <- cbind(as(taxa_rc_sp,'matrix'), as(uniqueseqs_rc,'data.frame'))    # input unique seqs as final row/col in taxa_rc_sp object
taxa_rc_sp <- as.matrix(taxa_rc_sp)                                             # back to matrix


##   C O N S T R U C T   F R O M   D A D A 2
# sample.names <- rownames(seqtab.nochim)       # remake as was in DADA2 flow
# Config <- sapply(strsplit(sample.names,'-'),'[',2)
#timepoint <- sapply(strsplit(sample.names,'-'),'[',3)
# mv_env2$timepoint <- c('one','two','three','four','five','six',
#             'one','two','three','four','five','six','seven','eight','nine',
#             'one','two','three','four','five','six','seven','eight','nine',
#             'one','two','three','four','five','six','seven','eight','nine')
# mv_env2$timepoint <- c(1:6,1:9,1:9,1:9)


## original MV-ENV fabricant
  # mv_env <- data.frame(Config=Config, timepoint=timepoint)
  # mv_env$situ <- "Ex-Situ"
  # mv_env$situ[mv_env$Config=="ISA"] <- "In-Situ"
  # mv_env$situ[mv_env$Config=="ISB"] <- "In-Situ"
  # mv_env$Reactor <- NULL
  # mv_env$Reactor[mv_env$Config=="ISB"] <- "2.ISB"
  # mv_env$Reactor[mv_env$Config=="ISA"] <- "1.ISA"
  # mv_env$Reactor[mv_env$Config=="BES"] <- "3.BES"
  # mv_env$Reactor[mv_env$Config=="CES"] <- "4.CES"
  # rownames(mv_env) <- samples.out
## hereon use expanded data, which is a little awkward and largely guessed at (ESPECIALLY H2 Rates: ISB 7-9)
mv_env2 <- read.table('mv_env_exp.txt',sep='\t',header=TRUE,row.names = 1)
mv_env2$H2[is.na(mv_env2$H2)] <- c(74)     ## guessed at H2 Rates: ISB 7-9, i.e. double medium rate
mv_env2$CO2[is.na(mv_env2$CO2)] <- c(0)     ## Totally invented CO2 Rates: na = 0?
mv_env <- sample_data(mv_env2)
# or jsut schuck direct into MV
# sample_data(mv) <- sample_data(mv_env2)

# Construct phyloseq object (straightforward from dada2 outputs), note Reverse-Complimented taxa w. SPec
mv_orig <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
#             read_tree('mv_raxml-ng.tre'),   # if available
              sample_data(mv_env),
              tax_table(taxa_rc_sp)
               )  
mv_orig
##  H O O R A Y
tax_table(mv_orig)[,1:8]
tax_table(mv_orig) <- tax_table(mv_orig)[,-c(8)]  #not sure why this is happening, but.

    # # another one for RAXMLng inclusion
    # mv_rax <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
    #                     read_tree('mv_raxml-ng.tre'),   # if available
    #                     sample_data(mv_env),
    #                     tax_table(taxa_rc_sp)
    # )  ; mv_rax   # note, is 2444 not 2300
  
write.table( cbind( t(otu_table(mv_orig)) , as.data.frame(tax_table(mv_orig)) ) ,'mv_abundance_taxonomy.txt',sep='\t')


## check for those lost Eukaryota
sum(sample_sums(subset_taxa(mv_orig,Kingdom=='Eukaryota')))
#> 2780  ## OR ##  #> 2902 (???) #> ORRR 2780 again


## R E M O V E   S V 's   B E L O W   T H R E S H O L D
z.txlen <- ncol(tax_table(mv_orig))  # allow index to end of tax_table, i.e. UniqueSeqs
## now check for and cut out SVs with SV length below e.g. 350 (based on histogram of SV lengths)
# SVs below (e.g.) 350
z.offcuts <- otu_table (mv_orig)[ , which (nchar (tax_table (mv_orig)[ , z.txlen ] ) <350), ]
sum( taxa_sums(z.offcuts))
#> 2752
# only keep taxa NOT found in z.offcuts - note inversion of logical! 
# previous approach was to re-invoke cutoff, prune (keep) those GREATER than 350 
mv <- prune_taxa(!(colnames(otu_table(mv_orig)) %in% colnames(z.offcuts)),mv_orig)
mv
sum(sample_sums(subset_taxa(mv,Kingdom=='Eukaryota')))
#> 528 # good enough?

    ## R E P E A T   for MV_RAX - consistency/compatability
#    mv_rax <- prune_taxa(!(colnames(otu_table(mv_rax)) %in% colnames(z.offcuts)),mv_rax) ; sum(sample_sums(subset_taxa(mv_rax,Kingdom=='Eukaryota')))
    ## should be 528

    
#  T A X   F I X E S: (for plotting purposes)
# see below for all the MBAO3 SV:
tax_table(subset_taxa(mv,Order=='MBA03'))[,1:6]
# fix largest MBAO3 SV in two steps: 
z.tx <- rownames(tax_table(subset_taxa(mv,Order=='MBA03')))
# need to index to them directly rather than apply/subsetting
tax_table(mv)[z.tx,5] <-c('MBA03 F.')
tax_table(mv)[z.tx,6] <-c('MBA03 G.')
# trying to set in the same call gets weird.


## largest Lentimicrobiaceae SV
z.tx <- rownames(tax_table(subset_taxa(mv,Family=='Lentimicrobiaceae'))[,1:7])
tax_table(mv)[z.tx,1:6]
# fix largest Lenti SV: 
tax_table(mv)[z.tx,6] <- c('Lentimicrobiaceae G.')


## Anaerolineaceae group in CES
tax_table(subset_taxa(mv,Family=='Anaerolineaceae'))[,1:7]
z.tx <- rownames(tax_table(subset_taxa(mv,Family=='Anaerolineaceae'))[is.na(tax_table(subset_taxa(mv,Family=='Anaerolineaceae'))[,6]),6]) 
tax_table(mv)[z.tx,6] <-c('Anaerolineaceae G.')   #fancier methods seem to fail


## Sequencing depth
sample_data(mv)$'SeqDepth' <- sample_sums(mv)
## Fixed timepoint
sample_data(mv)$timepoint <- as.factor(sample_data(mv)$timepoint)






## T R A N S F O R M   -   R E L . B U N
mv_rb <- transform_sample_counts(mv, function(x)x/sum(x))
## T R A N S F O R M   -   H E L L I N G E R
## hellinger of abund_value_i is the SQRT( abund_value_i / sum_all_abunds_site_j)
mv_h<-transform_sample_counts(mv, function(x) sqrt(x/sum(x)))   # hellinger   ## NO IT'S NOT   ### but is it tho  ## NO   ### but ...is it?  # it is, but is it the right way around? ### assume so as its a phylofunction
    # ## indications that row/column order of matrix affects transformation - do it in decostand:
    # otu_table(mv_h)[1:10,1:10]
    # dim(otu_table(mv)) # taxa are columns, samples are rows
    # z_h <- decostand(otu_table(mv), MARGIN=1, method='hellinger') #Margin=1 is the defaultm just being explicit
    # z_h[1:10,1:10]
    # # they're identical! even if identical() says otherwise.
    # summary( z_h[1:10,1:10] - otu_table(mv_h)[1:10,1:10] )   # see, all zeroes



    # ## T R A N S F O R M   -   R - L O G
    # library(DESeq2) ; mv_d2 <-phyloseq_to_deseq2(mv, ~Reactor)
    # # DESeq2 is amazing and carries you when you feel low.
    # # observe the effects of the following normalisations.
    # mv_norm <- normTransform(mv_d2)
    # mv_rlog <- rlogTransformation(mv_d2)
    # mv_vst <- varianceStabilizingTransformation(mv_d2)
    # library(vsn) # use meanSdPlot() to pleasantly illustrate how well your counts are standardised
    # pdf('mv_norm_rlog_vst.pdf')
    # meanSdPlot(assay(mv_norm))
    # meanSdPlot(assay(mv_rlog))
    # meanSdPlot(assay(mv_vst))
    # dev.off()    # no doubt, all the craic in the world.


## S U B S E T  # why 0.23?..
prunA = genefilter_sample(mv, filterfun_sample(function(x) x >=200), A=0.23*nsamples(mv)) # different from subset used in heatmapping # (10, 0.23) misses taxa picked up by LEfSe
mv_10 = prune_taxa(prunA, mv)
mv_un10 = prune_taxa(!prunA, mv)

prunB = genefilter_sample(mv_rb, filterfun_sample(function(x) x >=0.01), A=3) # arbitrary rel.abundance cutoff: >2%  #,A=0.23*nsamples(mv_rb))
mv_01 = prune_taxa(prunB, mv_rb)
mv_un01 = prune_taxa(!prunB, mv_rb)    # N O T E: negation of logical vector with prepending '!' !!

prunC = genefilter_sample(mv_un01, filterfun_sample(function(x) x >=0.001), A=3) # arbitrary rel.abundance cutoff: >2%  #,A=0.23*nsamples(mv_rb))
mv_001 = prune_taxa(prunC, mv_un01)
mv_un001 = prune_taxa(!prunC, mv_un01)    # N O T E: negation of logical vector with prepending '!' !!



##   F I N I S H
## record the taxa and SVs
# write.table(tax_table(mv),'mv_tax.txt',sep='\t')
# write.table(t(otu_table(mv)),'mv_sv_transp.txt',sep='\t')
# write.table(mv_env,'mv_env.txt',sep='\t')

## K E P T O M A N I A
# rm everything EXCEPT your 'kept' material (including the 'kept' list!)
#kept <- c('mv','mv_orig','mv_env2','seqtab.nochim','seqtab.nochim.orig','taxa_rc','taxa_rc_sp','taxa.orig','uniqueseqs_rc')  # keep input
kept <- c('mv','mv_h','mv_orig','mv_rax','mv_env2','mv_10','mv_un10','mv_02','mv_un02','mv_01','mv_un01','mv_001','mv_un001','mv_rb')             # keep output
rm(list=setdiff(ls(),kept))



# =====================================================================================

#				. . . t o   p h y l o s e q

# =====================================================================================




##   P R I M E 
# load.image('phyloforeward.RData')
# hallowed, hollow ground
#install.packages('tidyr')
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library('tidyr')
library('DESeq2')
library('vegan')
library('scales')
library('heatmaply')

##   S E Q U E N C E   L E N G T H S - infer quality of identity
pdf('mv_distribution of sequences.pdf')
hist(as.vector(apply( tax_table(mv)[,8], 1, function(x) nchar(x) )),main='histogram of all SV lengths')
hist(as.vector(apply( tax_table(subset_taxa(mv,Kingdom=='Bacteria'))[,8], 1, function(x) nchar(x) )),main='histogram of Bact SV lengths')
hist(as.vector(apply( tax_table(subset_taxa(mv,Kingdom=='Archaea'))[,8], 1, function(x) nchar(x) )),main='histogram of Arch SV lengths')
hist(as.vector(apply( tax_table(subset_taxa(mv,Kingdom=='Eukaryota'))[,8], 1, function(x) nchar(x) )),main='histogram of "Euk" SV lengths')
## Eukaryota seqs are alll over the place! That one Archaeal seq ~230 (M'bacterium, total abundance of 5) doesn't look great either.
taxa_sums(subset_taxa(mv,Kingdom=='Eukaryota'))   # 2902 in total, but first 15 comprise 2039.
sum(taxa_sums(subset_taxa(mv,Kingdom=='Eukaryota'))[1:30])
tax_table(subset_taxa(mv,Kingdom=='Eukaryota'))[1:30,1:6]          # empty taxa you dope
## see euk specific plot later
z.hieuk <- names(taxa_sums(subset_taxa(mv,Kingdom=='Eukaryota')))[1:30]
z.hieuk <- prune_taxa(z.hieuk, mv) ; sum(sample_sums(z.hieuk)) ; plot_bar(z.hieuk,fill='UniqueSeq') +
  theme(legend.position="none") +
  facet_grid('~Reactor',scales='free_x',space='free_x')
rownames(tax_table(z.hieuk))
dev.off()


#dev.off()





## S U M M A R Y   S T A T S  ##
# phylum breakdown across whole study, per-reactor, and maybe even per sample...
# not the same as the average % per sample, i.e. irrespective of library size.
#
##   M A K E   M E   I N T O   A   F U N C T I O N
### take phylo object, take sample variable: cut phylo by variable, by rank
## A L S O   A   F U N C T I O N   F O R   O T U   I D s

mv_phyglom <- tax_glom(mv,taxrank='Phylum')  # one glom instance preserves order in later sub-sets

# sum_phy inherentl ordered by SV abundance
sum_phy <- as.data.frame( round( (taxa_sums(mv_phyglom)/sum(taxa_sums(mv) ) *100), digits=4 ),  # grab rounded percentages
                          row.names = as.character(tax_table(mv_phyglom)[,2]) )                 # w. phyla as row names , cant set col.names?

# can illustrate taxa tables are intact based on sort(taxa_sums(z.isa/b/bes/ces)) etc.
z.isa <- subset_samples(mv_phyglom,Reactor=='1.ISA' )
sum_phy[,2] <- round( (taxa_sums(z.isa)/sum(taxa_sums(z.isa) ) *100), digits=4 )  #keep NArm to keep comparable order between reactors
z.isb <- subset_samples(mv_phyglom,Reactor=='2.ISB' )
sum_phy[,3] <- round( (taxa_sums(z.isb)/sum(taxa_sums(z.isb) ) *100), digits=4 )
z.bes <- subset_samples(mv_phyglom,Reactor=='3.BES' )
sum_phy[,4] <- round( (taxa_sums(z.bes)/sum(taxa_sums(z.bes) ) *100), digits=4 )
z.ces <- subset_samples(mv_phyglom,Reactor=='4.CES' )
sum_phy[,5] <- round( (taxa_sums(z.ces)/sum(taxa_sums(z.ces) ) *100), digits=4 )

colnames(sum_phy) <- c('TOT%','ISA%','ISB%','BES%','CES%') # cant seem to sort (order) properly either

## add number of genera per phylum (overall) #OR just the 3-5 phyla you'll actually talk about..
#length(get_taxa_unique(subset_taxa(mv,Phylum=='Firmicutes'),'Genus')) # this will overlook NAs etc
dim(tax_table(subset_taxa(mv,Phylum=='Firmicutes')))
dim(tax_table(subset_taxa(mv,Phylum=='Euryarchaeota')))
dim(tax_table(subset_taxa(mv,Phylum=='Bacteroidetes')))

View(sum_phy)

## Stats per mv_01 / mv_un01
# total ASV_orig, total sequences_orig
mv_orig ; sum(sample_sums(mv_orig)) ; mean(sample_sums(mv_orig))
# total ASV, total sequences
mv ; sum(sample_sums(mv)) ; mean(sample_sums(mv))
# reads / sample
sum_phy 
# total un_ASV, total un_sequences
mv_01 ; sum(sample_sums(otu_table(mv)[,rownames(tax_table(mv_01))])) ; mean(sample_sums(otu_table(mv)[,rownames(tax_table(mv_01))]))
# total un_ASV, total un_sequences
mv_un01 ; sum(sample_sums(otu_table(mv)[,rownames(tax_table(mv_un01))])) ; mean(sample_sums(otu_table(mv)[,rownames(tax_table(mv_un01))]))
# Firm / Bact / Eury % range
sum_phy # see above
# #firm / #bact / #eury
dim(tax_table(subset_taxa(mv,Phylum=='Firmicutes')))
dim(tax_table(subset_taxa(mv,Phylum=='Euryarchaeota')))
dim(tax_table(subset_taxa(mv,Phylum=='Bacteroidetes')))

# NON- Firm / Bact / Eury % range
c(100, 100, 100, 100, 100) - colSums(sum_phy[1:3,])
get_taxa_unique(mv,taxonomic.rank = 'Phylum') ; length(get_taxa_unique(mv,taxonomic.rank = 'Phylum'))
#subtract this from 2300:
(dim(tax_table(subset_taxa(mv,Phylum=='Firmicutes'))) ) + ( dim(tax_table(subset_taxa(mv,Phylum=='Euryarchaeota'))) ) + ( dim(tax_table(subset_taxa(mv,Phylum=='Bacteroidetes'))) )
2300 - 1662
# 




pdf('mv_dada_phylo_rc.pdf',width=10, height=7.5)


###   P L O T T I N G

# C O L O U R S   -   B Y   S E Q V A R
# more illustrative to assign a randomised colour vector of the correct length to catch all SVs
library(scales,RColorBrewer)
z.col <- hue_pal()( nrow(tax_table(mv_01)) )  ;  z.col1 <- sample(z.col1, nrow(tax_table(mv_01)) )
z.col1 <- hue_pal(l=55, c=80)( nrow(tax_table(mv_01)) )  ;  z.col1 <- sample(z.col1, nrow(tax_table(mv_01)) )
z.col2 <- hue_pal(l=50, c=85)( nrow(tax_table(mv_un01)) )  ;  z.col2 <- sample(z.col2, nrow(tax_table(mv_un01)) )
show_col((z.col)[1:37])
# diversity metrics: 
hist(as.vector(apply( tax_table(mv)[,8], 1, function(x) nchar(x) )),main='histogram of all SV lengths (<350 trimmed)')
z.rm <- plot_richness(mv, x = "Reactor",
               color = "SeqDepth",
               shape = 'Reactor',
               measures=c('Observed',
                          'Shannon',
                          'InvSimpson')
               ) +
  geom_point(size=3) +
  scale_shape_manual(values = c(17, 15, 16, 18)) +
  scale_fill_distiller(palette='Spectral')
  
#+ geom_boxplot()


## B A R   C H A R T S
# dont need grid lines now?
# bw theme?
# plot_bar(mv,'Reactor',fill='Phylum',title='MV Absolute Read Abundance per Reactor')
plot_bar(mv,'timepoint',fill='Phylum',title='MV Absolute Read Abundance: 6,114,729 of 7,941,206 reads total') +
  facet_grid(~Reactor,scales='free_x',space='free_x') +
  scale_fill_manual(values = z.col42)

plot_bar(mv_rb,'timepoint',fill='Phylum',title='MV Relative Read Abundance - Phylum Breakdown') +
  facet_grid(~Reactor,scales='free_x',space='free_x')
plot_bar(mv_rb,'timepoint',fill='UniqueSeq',title='MV Relative Read Abundance - Decorative Purposes') +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +     # rotate
  facet_grid(~Reactor,scales='free_x',space='free_x') +
  scale_fill_manual(values=z.col) +
  theme(legend.position="none")

  # ## PHYLOGENETICS
  # # tree composed in RAXML-ng
  # # does subsetting frack the tree?
  # # pdf('mv_trees.pdf')
  # plot_tree(mv_orig, color = 'Phylum',title = 'RAXMLng Tree of MV (2444 SVs)', ladderize = TRUE) + theme(legend.position="none")
  # plot_tree(mv, color = 'Phylum', title='RAXMLng Tree of MV length filtered (>350bp, 2300 SVs)')
  # plot_tree(mv_10, color = 'Phylum', title = 'RAXMLng Tree of MV >10reads in 50% prevalence (197 SVs)')
  # # dev.off()

###   M E L T I E S
# + geom_bar(aes(color=CATEG, fill=CATEG), stat ="identity", position="dodge")
## Melt, unify, and send to one dataframe object
z.rm3 <- psmelt(mv_01)
z.rm4 <- psmelt(mv_un01)
# sort mv_un02
z.rm4 <- z.rm4[with(z.rm4, order(Kingdom, Phylum, Class, Order, Family, Genus)), ]
#plot the agglomerated SVs for reference of phylum abundance
## pool all non->2% to dummy phylum, append to bottom of melted frame 
z.rm4$Kingdom <- c('D: Misc')
z.rm4$Phylum <- c('P: Misc')
z.rm4$Class <- c('C: Misc')
z.rm4$Order <- c('O: Misc')
z.rm4$Family <- c('F: Misc')
z.rm4$Genus <- c('G: Misc')
z.rm4$Species <- c('Sp: Misc')
#z.rm4$UniqueSeq <- c('SV: Misc')
z.rm5 <- rbind(z.rm3, z.rm4)


## C O M B O   b y   S e q V a r  -  more space - favoured plot
ggplot(z.rm5, aes(x=timepoint,y=Abundance,fill=UniqueSeq)) +
  #  facet_grid(Phylum+Class+Order+Family+Genus~Reactor,scales='free_x',space='free') +                                  # scale/size control
  facet_grid(Phylum+Genus~Reactor,scales='free',space='free') +                                  # scale/size control
  geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack",  colour="black") +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 9)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +  # rotate
  scale_fill_manual(values = (c(sort(z.col1), (z.col2)))) +  #improve SV differentiation #c(z.col42[1:37]) , 
  theme(panel.background = element_rect(fill='white', colour='grey70')) +
  theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
  theme(panel.grid.minor.y = element_line(colour='grey85', size = 0.1)) +                               # horiz lines 2
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
  ggtitle("Major ASVs (Rel.Ab >1%, note variable scale)") + theme(legend.position="none")

ggplot(z.rm5, aes(x=timepoint,y=Abundance,fill=Genus)) + 
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") + 
  facet_grid(Class+Order+Family~Reactor,scales='free',space='free') + 
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) + 
  ggtitle("DW SLV OTUs with Rel.Abundances above 2% by Genus/Order/Class") +
  theme(legend.position="none") 


## Ladies and Gentlemen; THE ARCHAEA!
# re-scale colours to length of SVs in Archaea
z.colA <- sample(z.col, nrow(tax_table(subset_taxa(mv_rb,Kingdom=='Archaea'))) )
z.rmA <- psmelt(subset_taxa(mv_rb,Kingdom=='Archaea'))
ggplot(z.rmA, aes(x=timepoint,y=Abundance,fill=z.col)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Order+Genus~Reactor,scales='free',space='free_x') +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 11)) +      # biggerness
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +     # rotate
  ggtitle("MV Archaeal SVs by Class:Genus")  + theme(legend.position="none")
## ...   *clap*,     *clap*.

      # 
      # ## C O M B O M I N O R S 
      # z.rm4 <- psmelt(mv_un02)   # reset z.rm4
      # ggplot(z.rm4, aes(x=timepoint, y=Abundance, fill=UniqueSeq)) +
      #   geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack",  colour="black") +
      #   facet_grid(Kingdom~Reactor,scales='free',space='free_x') +                                     # scale/size control
      #   theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +      # bigger
      #   theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +     # rotate
      #   theme(legend.position="none") + ggtitle("MV Overview of 'Minor' SVs Rel.Abundances (<2%)")
      # 
      # ## C O M B O M I N O R B A C T E R I A 
      # z.rm4 <- psmelt(subset_taxa(mv_un02,Kingdom=='Bacteria'))   # small bacterial taxa
      # ggplot(z.rm4, aes(x=timepoint, y=Abundance, fill=UniqueSeq)) +
      #   geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack",  colour="black") +
      #   facet_grid(Phylum~Reactor,scales='free',space='free_x') +                                     # scale/size control
      #   theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +      # bigger
      #   theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +     # rotate
      #   theme(legend.position="none") + ggtitle("MV Overview of 'Minor' Bacterial SVs Rel.Abundances (<2%)")
      # 
      # ## C O M B O M I N O R P: F I R M I E S 
      # z.rm4 <- psmelt(subset_taxa(mv_un02,Phylum=='Firmicutes'))   # small bacterial taxa
      # ggplot(z.rm4, aes(x=timepoint, y=Abundance, fill=UniqueSeq)) +
      #   geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack",  colour="black") +
      #   facet_grid(Class~Reactor,scales='free',space='free_x') +                                     # scale/size control
      #   theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +      # bigger
      #   theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +     # rotate
      #   theme(legend.position="none") + ggtitle("MV Overview of 'Minor' P:Firmicutes  SVs Rel.Abundances (<2%)")
      # 
      # ## C O M B O M I N O R C: C L O S T 
      # z.rm4 <- psmelt(subset_taxa(mv_un02,Class=='Clostridia'))   # small bacterial taxa
      # ggplot(z.rm4, aes(x=timepoint, y=Abundance, fill=UniqueSeq)) +
      #   geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack",  colour="black") +
      #   facet_grid(Order~Reactor,scales='free',space='free_x') +                                     # scale/size control
      #   theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +      # bigger
      #   theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +     # rotate
      #   theme(legend.position="none") + ggtitle("MV Overview of 'Minor' C:Clostridial SVs Rel.Abundances (<2%)")
      # 
      # ## C O M B O M I N O R O: C L O S T R I D I A L E S  - Caldicop, FamilyXI, Rumino, Syntropho, Christen, 
      # z.rm4 <- psmelt(subset_taxa(mv_un02,Order=='Clostridiales'))   # small bacterial taxa
      # ggplot(z.rm4, aes(x=timepoint, y=Abundance, fill=UniqueSeq)) +
      #   geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack") + # remove colour="black" to de-outline
      #   facet_grid(Family~Reactor,scales='free_x',space='free_x') +                                     # scale/size control
      #   scale_fill_manual(values = (z.col)) +     # manual FILL; z.col invoked above
      #   theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +      # bigger
      #   theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +     # rotate
      #   theme(legend.position="none") + ggtitle("MV Overview of 'Minor' O:Clotridiales SVs Rel.Abundances (<2%)")

## see mv_minorplots.R for more of this nonsense


# ## " E U K " - not necessary, informative etc.
# plot_bar(subset_taxa(mv,Kingdom=='Eukaryota'),fill='UniqueSeq',title='MV "Eukaryote" SV plot') +
#   theme(legend.position="none") +
#   facet_grid('~Reactor',scales='free_x',space='free_x')



### O R D I N A T E
plot_ordination(mv,ordinate(mv,'NMDS',distance='bray'),
                color='Reactor', shape='situ', title='NMDS(B-C) of MV - Samples',
                label='null', type='samples')
plot_ordination(mv,ordinate(mv,'CCA'), #,formula= (~'Acetic.acid')
                color='Phylum', shape='Reactor', title='CCA of MV - Separating Samples by Community Structure', type='split') #+ 

plot_ordination(mv, ordinate(mv, 'CCA', formula=(~H2 )), #+ CO2
                color='nominal.rate', shape='Reactor', title='CCA of MV - H2', type='samples') #+ 

plot_ordination(mv_h,ordinate(mv,'DCA'),#,formula= (~'Acetic.acid')),
                color='SeqDepth', shape='Reactor', title='CCA of MV - Samples by SeqDepth - No obvious clustering by coverage', type='samples') +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  geom_point(size=5)


dev.off()
# rm(list=ls(pattern='z.'))
