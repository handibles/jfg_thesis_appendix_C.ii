
# # ##                        ## # #
# # ##   T  E  S  T  S    3   ## # #  
# # ##                        ## # #

# see tests_1 for non-correlation tests, tests_2 for mental correlation matrix roundaboutage

## DCA over PCA (indirect gradient analysis)
## H-transform to account for arch (rather than horseshoe) effect

## CCA - model building in step / AIC / Akaike Information Content
# from the vegantutor.pdf
# knock all the NAs out of the sampledata:
z.mv_env3 <- mv_env2              # duplicate
z.mv_env3[is.na(z.mv_env3)] <- 0  # kill NAs 
View(z.mv_env3)                    # View(as.data.frame(as.matrix(sample_data(mv))))

# kill all the dross too
names(z.mv_env3)
z.mv_env3 <- z.mv_env3[,-c(1, 2, 12, 13, 24)]

z.mod1 <- cca(otu_table(mv) ~., z.mv_env3)  # model with everything...
z.mod0 <- cca(otu_table(mv) ~1 , z.mv_env3) # model with nothing...
z.mod <- step(z.mod0, scope = formula(z.mod1), test = "perm")

sample_data(mv_h)$'SeqDepth' <- sample_sums(mv)

# SeqDepth CCA
# plot_ordination(mv_h,ordinate(mv_h, 'CCA', formula=(~SeqDepth)),
#                 color='SeqDepth', shape='Reactor', title='CCA of MV - Samples by SeqDepth - No obvious clustering by coverage',
#                 type='Samples') +
#   scale_shape_manual(values = z.col42) +
#   geom_point(size=5)

# CCA setup (i.e. IS/B/C or IS/ES versus H2 concentration - how much variation accounted for?)
sample_data(mv_cca)$timepoint <- as.factor(sample_data(mv_h)$timepoint)

#c(15, 16, 17, 18, 1)




## fancy ordinations
# 

tax_table(mv_un10)[,1:7] <- c('Misc <10 reads')
mv_cca <- merge_phyloseq(mv_un10, mv_10)
mv_ccah <- transform_sample_counts(mv_cca, function(x) sqrt(x/sum(x)))   # hellinger

# cant figure out how to cast trasnparent or less offensive species dots. smaller points?
z.col42[16] <- c('#6CEAEA') # Order
z.col42[12] <- c('#98C4F0') # Family
#z.deed <- rgb(108, 234, 234, max = 255, alpha = 50, names = "goopblue")
show_col(z.col42)
z.col42[12]
z.col42[12] <- c('grey85')

plot_ordination(mv_cca,ordinate(mv_ccah,'DCA'),  #, formula = (~H2)
                color='SeqDepth',
                shape='Reactor',
                title='DCA of Samples Emphasising Library Size (SeqDepth)',
                type='samples'
) + 
#  scale_color_viridis_d() + #_colour_manual(values = (z.col42)) +
  geom_point(size=5) +  #, alpha=1/10
  scale_shape_manual(values = c(15, 20, 17, 18)) +
  theme(panel.background = element_rect(fill='white', colour='grey20')) +
  theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +         # horiz lines 1
  theme(panel.grid.major.x = element_line(colour='grey75', size = 0.2))          # horiz lines 2
  
# geom_point(size=5)





## ANCOM
# plotting results, not doing the test (test is in continued support)
z.anc <- c('mvsv0020','mvsv0024','mvsv0044','mvsv0051','mvsv0060',
           'mvsv0065','mvsv0085','mvsv0090','mvsv0113','mvsv0119',
           'mvsv0138','mvsv0159','mvsv0186','mvsv0225','mvsv0227',
           'mvsv0247','mvsv0221','mvsv0037','mvsv0064','mvsv0097',
           'mvsv0114','mvsv0127','mvsv0187','mvsv0218','mvsv0222',
           'mvsv0380','mvsv0010')
z.anc <- prune_taxa(z.anc,mv_rb)
plot_bar(z.anc,'timepoint',fill='UniqueSeq',title='MV IS/ES ANCOM.i') +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +     # rotate
  facet_grid(Class+Genus~Reactor,scales='free',space='free_x') +
  scale_fill_manual(values=z.col) +
  theme(legend.position="none")
#
# STRIKINGLY LACKING
#


## L E F S E 
# copy lefse names to gedit, remove non mvsv stuff with Ctrl H: \S*mv, save as table
as.vector(z.sv) <- read.table('lefse_sv.txt',header=FALSE)
rownames(z.sv) <- z.sv[,1] # there MUST be a better way of doing this, but..
mv_lef <- prune_taxa(colnames(otu_table(mv)) %in% rownames(z.sv),mv_rb)

ggplot(psmelt(mv_lef), aes(x=timepoint,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Order+Genus~Reactor,scales='free',space='free') +                                  # scale/size control
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 10)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +  # rotate
  ggtitle("LEfSe SV's >3 by Genus (RB)") + theme(legend.position="none")



## D E S e q 2   mk.iii
# 
pdf(width=8, height=6,'mv1000__I-E__B-C_0.01.pdf')
# can play with lfcShrink here if intersted in ranking and vis, as default shrinkage can
# be overly strong for some datasets
## results() does auto indpendent filtering
## alpha is the adjusted p value (padj) cutoff
  ## L R T
  #https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#likelihood-ratio-test
  # LRT is different rom usual DESeq2 test (Wald) in that it comapres two models;
  # with (full)  and without (reduced) your term of interest t odetermine effect of
  # the removed/duced term.
  
  # useful when comparing a factor with >2 levels (i.e. Reactor), would use:
  #mv_in_r<-DESeq(mv_in, test= 'LRT', reduced=~1, fitType='local', minReplicatesForReplace =5)
  # where reduced= ~1 is the 'basic' instance of model.
  # Don't need reduced model if comapring in situ, as only two cases (in/ex) so using Wald. 


mv_in <- phyloseq_to_deseq2(mv, ~situ)
# prefiltering to remove outgrageous LFC from near absent taxa.
keep <- rowSums(counts(mv_in)) >= 1000
mv_in <- mv_in[keep,]

# use relevel to specify the base condition, otherwise works alphabetically (i.e. 'ex-situ' = baseline)
## W A L D
mv_in$situ <-relevel(mv_in$situ,ref='In-Situ')
mv_in_s<-DESeq(mv_in) #fitType='local',minReplicatesForReplace =5)

mv_in_r1_res <- results(mv_in_s, cooksCutoff=TRUE, alpha=0.01, lfcThreshold=1) #, contrast=c('setup','CES','IS'))   #
mv_in_r1_sig <- mv_in_r1_res[which(mv_in_r1_res$padj<0.01),]
mv_in_r1_sig <- cbind(as(mv_in_r1_sig,'data.frame'),as(tax_table(mv)[rownames(mv_in_r1_sig),],'matrix'))
ggplot(mv_in_r1_sig, aes(x=Genus, y=log2FoldChange, color= Phylum)) +
  geom_hline(yintercept=0, color='grey75') +
  theme(axis.text = element_text(size=12)) +
  labs(title='Differential Abundance: in situ (L) v. ex situ (R)   (p/padj=0.01/0.01)') +
  coord_flip() + scale_size_continuous(name='root sum abundance') +
  scale_colour_manual(values = (z.col57[seq(0,58,4)])) +  #custom fill
  geom_point(alpha=0.6, aes(size=( sqrt(33*baseMean) ) ), width=0.1, height=0.1) +  # see also geom_point
  geom_point(color='grey35',size=0.1) #+
#  guides(col=guide_legend(ncol=2))
# lfcThreshold really curtails results.. especially Methanothermobacter!

## export DES2 data
write.table( cbind(as(mv_in_r1_res,'data.frame'),as(tax_table(mv)[rownames(mv_in_r1_res),],'matrix')), 'mv_in_r1_res.txt', sep='\t')



# catch NAs
## Catch and explore the NAs in another plot: not really anything to see here! :)
z.rm <- mv_in_r1_sig[which(is.na(mv_in_r1_sig$Genus)),]
g.d2_isesNA <- ggplot(z.rm, aes(x=Phylum, y=log2FoldChange, color= Family)) +
  geom_hline(yintercept=0, color='grey75') +
  theme(axis.text = element_text(size=12)) +
  labs(title='DA: in situ (L) v. ex situ (R) NA\'s   (p/padj=0.01/0.01)') +
  coord_flip() + scale_size_continuous(name=' total abundance') +
  scale_colour_manual(values = (sample(z.col,22)), na.value='red') +  #custom fill, custom NAs
  geom_jitter(alpha=0.6, aes(size=((15*baseMean)) ), width=0.1, height=0.1) +  # see also geom_point
  geom_point(color='grey35',size=0.8) #+


    # shows a dichotomy in results, likely taxa only in one of the conditions. 
#View(mv_in_r1_sig) # note -0-9, then 20-25 (-/+)
    # T O T A L L Y   R E L A T I V E 
      # cut out middlemen with lfcThreshold.
      z.rm <- c(order(mv_in_r1_sig$log2FoldChange)[191:233], order(mv_in_r1_sig$log2FoldChange)[1:11])
      mv_new <- mv_in_r1_sig[z.rm,]
      ggplot(mv_new, aes(x=Family, y=log2FoldChange, color= Order)) +
        geom_hline(yintercept=0, color='grey75') +
        theme(axis.text = element_text(size=12)) +
        labs(title='MV-Reactor - g.TE v. silage (LRT, 0.05, local)') +
        coord_flip() + scale_size_continuous(name='root sum abundance') +
        geom_jitter(alpha=0.6, aes(size=( sqrt(33*baseMean) ) ), width=0.1, height=0.1) +  # see also geom_point
        geom_point(color='grey35',size=0.1) +
        scale_fill_manual(values = (z.col)) +
        guides(col=guide_legend(ncol=2))
      
      mv_old <- mv_in_r1_sig[order(mv_in_r1_sig$log2FoldChange)[12:190],] 
      ggplot(mv_old, aes(x=Family, y=log2FoldChange, color= Order)) +
        geom_hline(yintercept=0, color='grey75') +
        theme(axis.text = element_text(size=12)) +
        labs(title='MV-Reactor - g.TE v. silage (LRT, 0.05, local)') +
        coord_flip() + scale_size_continuous(name='root sum abundance') +
        geom_jitter(alpha=0.6, aes(size=( sqrt(33*baseMean) ) ), width=0.1, height=0.1) +  # see also geom_point
        geom_point(color='grey35',size=0.1) +
        scale_fill_manual(values = (z.col)) +
        guides(col=guide_legend(ncol=2))


##  SUBSETTING:  E X   S I T U
##
# consider culling the majority of SVs: 100 top SVs are 90% of sequences: 10% ~=300,000 reads
  #sum(taxa_sums(otu_table(mv)[,1:100]))/ sum(taxa_sums(mv))
  # however, not appropriate for ES only as many taxa abundant in IS not ES
  # filter z.rm for > 1% etc. 

z.rm <- subset_samples(mv, situ=='Ex-Situ')
prunB <- genefilter_sample(z.rm, filterfun_sample(function(x) x >0)) ; z.rm <- prune_taxa(prunB, z.rm)
z.rb <- transform_sample_counts(z.rm, function(x)x/sum(x))

prunB = genefilter_sample(z.rb, filterfun_sample(function(x) x >=0.01), A=2) # arbitrary rel.abundance cutoff: >2%  #,A=0.23*nsamples(mv_rb))
z.02 = prune_taxa(prunB, z.rb)
z.un02 = prune_taxa(!prunB, z.rb)    # N O T E: negation of logical vector with prepending '!' !!

z.rm3 <- psmelt(z.02) ; z.rm4 <- psmelt(z.un02)
z.rm4[,] <- c('Misc <1%')
z.rm5 <- rbind(z.rm3, z.rm4)

#pdf(width=8, height=6,'mv1000_B-C_0.01.pdf')

z.col <- hue_pal(l=50)( nrow(tax_table(z.rb)) )  ;  z.col <- sample(z.col, nrow(tax_table(z.rb)) )
ggplot(z.rm3, aes(x=timepoint,y=Abundance,fill=UniqueSeq)) +
  facet_grid(Phylum+Family+Genus~Reactor,scales='free',space='free') +                                  # scale/size control
  geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack",  colour="black") +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 9)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +  # rotate
  scale_fill_manual(values = (z.col)) +  #improve SV differentiation
  ggtitle("Major ex situ ASvs (Rel.Ab >1%, note variable scale)") + theme(legend.position="none")
  # add BG
  # add lines


##  D E S E Q 2 :   C E S   v.   B E S 
##
mv_es <- phyloseq_to_deseq2(z.rm, ~Reactor)
# plot((rowSums(counts(mv_es)))) ; summary((rowSums(counts(mv_es)))) 
keep <- rowSums(counts(mv_es)) >= 1000
mv_es <- mv_es[keep,]
# use relevel to specify the base condition, otherwise works alphabetically (i.e. 'ex-situ' = baseline)
mv_es$Reactor <-relevel(mv_es$Reactor,ref='3.BES')
mv_es_s<-DESeq(mv_es) #fitType='local',minReplicatesForReplace =5)
## Wald test with a=0.05 on BES/CES
mv_es_r1_res <- results(mv_es_s, lfcThreshold=1, alpha=0.01)  
                        # altHypothesis= c('greater'),
                        # contrast=c('setup','CES','BES'),
                        # cooksCutoff=meaningless-for-2-groups, ) 
mv_es_r1_sig <- mv_es_r1_res[which(mv_es_r1_res$padj<0.01),]
mv_es_r1_sig <- cbind(as(mv_es_r1_sig,'data.frame'),as(tax_table(mv)[rownames(mv_es_r1_sig),],'matrix'))


## Export
write.table( cbind(as(mv_es_r1_res,'data.frame'),as(tax_table(mv)[rownames(mv_es_r1_res),],'matrix')), 'mv_es_r1_res.txt', sep='\t')

## MV BES v. CES 
  ggplot(mv_es_r1_sig, aes(x=Genus, y=log2FoldChange, color= Phylum)) +
    geom_hline(yintercept=0, color='grey75') +
    theme(axis.text = element_text(size=12)) +
    labs(title='MV Setups:BES v. CES   0.01/0.01') +
    coord_flip() + scale_size_continuous(name=' total abundance') +
    scale_colour_manual(values = (sample(z.col42,42)), na.value='red') +  #custom fill
    geom_jitter(alpha=0.6, aes(size=((15*baseMean)) ), width=0.1, height=0.1) +  # , shape=Kingdom
    geom_point(color='grey35',size=0.5) #+


## Catch and explore the NAs in another plot: not really anything to see here! :)
  z.rm <- mv_es_r1_sig[which(is.na(mv_es_r1_sig$Genus)),]
  ggplot(z.rm, aes(x=Phylum, y=log2FoldChange, color= Order)) +
    geom_hline(yintercept=0, color='grey75') +
    theme(axis.text = element_text(size=12)) +
    labs(title='MV Setups:BES v. CES NA\'s   0.01/0.01') +
    coord_flip() + scale_size_continuous(name=' total abundance') +
    scale_colour_manual(values = (sample(z.col,22)), na.value='red') +  #custom fill, custom NAs
    geom_jitter(alpha=0.6, aes(size=((15*baseMean)) ), width=0.1, height=0.1) +  # see also geom_point
    geom_point(color='grey35',size=0.8) #+


  z.rm6 <- psmelt(prune_taxa(rownames(mv_es_r1_sig),z.rb))
  ggplot(z.rm6, aes(x=timepoint,y=Abundance,fill=UniqueSeq)) +
    facet_grid(Phylum+Family~Reactor,scales='free_x',space='free') +                                  # scale/size control
    geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack") + #,  colour="black") +
    theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 9)) +   # bigger
    theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +  # rotate
    scale_fill_manual(values = (z.col)) +
    ggtitle("MV DA in BES - CES (Rel.Ab >2%, note variable scale)") + theme(legend.position="none")
  
dev.off()

## catch those hydrogen metabolisers, trashypastying done in excel
z.h2 <- c('mvsv0065','mvsv0068','mvsv0072','mvsv0085',
          'mvsv0157','mvsv0189','mvsv0243','mvsv0257',
          'mvsv0031','mvsv0040','mvsv0098','mvsv0099',
          'mvsv0150','mvsv0155','mvsv0158','mvsv0181',
          'mvsv0183','mvsv0204','mvsv0219','mvsv0230',
          'mvsv0265','mvsv0273','mvsv0074','mvsv0133','mvsv0173')
rowSums(otu_table(mv_rb)[,z.h2])
#check content
tax_table(mv_rb)[z.h2,1:7]
  ## = = = = = = = = =
  ## = = = = = = = = =
  

## D E S E Q 2 :   H 2    (H2 shown to be fairly ineffective in CA/step mdoel building - better variable?)
##
mv_env3 <- mv_env2
names(mv_env3)
mv_env3$H2.rate.L.day
mv_env3[16:18,6]<-c(74,74,74)    # dont quote, will un-numeric it 
names(sample_data(mv_10))[[6]] <- c('H2')

mv_inh <- phyloseq_to_deseq2(mv_10, ~H2)
keep <- rowSums(counts(mv_inh)) >= 100 ; mv_inh <- mv_in[keep,]

#mv_inh$H2 <-relevel(mv_inh$H2,ref=0) # relevel is only for factors
mv_inh_s<-DESeq(mv_inh) #fitType='local',minReplicatesForReplace =5)
## W A L D
mv_inh_r1_res <- results(mv_inh_s, cooksCutoff=TRUE,alpha=0.05) #, contrast=c('Reactor','4.CES','2.ISB'),
mv_inh_r1_sig <- mv_inh_r1_res[which(mv_inh_r1_res$padj<0.05),]
mv_inh_r1_sig <- cbind(as(mv_inh_r1_sig,'data.frame'),as(tax_table(mv)[rownames(mv_inh_r1_sig),],'matrix'))
mv_inh_r1_sig <- mv_inh_r1_sig[-c(46),]
ggplot(mv_inh_r1_sig, aes(x=Family, y=log2FoldChange, color= Order)) +
  geom_hline(yintercept=0, color='grey75') +
  theme(axis.text = element_text(size=12)) +
  labs(title='MV-Reactor - g.TE v. silage (LRT, 0.05, local)') +
  coord_flip() + scale_size_continuous(name='Avg. Ab in situ') +
  geom_jitter(alpha=0.6, aes(size=(abs(baseMean*2^log2FoldChange)) ), width=0.1, height=0.1) +
  geom_point(color='grey35',size=0.1) +
  scale_fill_manual(values = (z.col)) +
  guides(col=guide_legend(ncol=2))





## D E S e q 2    ( for  P L O T T I N G   a n d   L E f S e  mk.ii )
# DESEq2 advises rlog is just for exploration, and to use DESeq2's main transformation when testing.
# what doe sthat look like?
mv_in_sD2 <- counts(mv_in_s,normalized=TRUE)s
mv_D2 <- mv
otu_table(mv_D2) <- otu_table(mv_in_sD2, taxa_are_rows =TRUE)
prunA = genefilter_sample(mv_D2, filterfun_sample(function(x) x >=500), A=0.2*nsamples(mv)) # different from subset used in heatmapping # (10, 0.23) misses taxa picked up by LEfSe
mv_10 = prune_taxa(prunA, mv) ; mv_10
mv_un10 = prune_taxa(!prunA, mv) ; mv_un10

z.rm3 <- psmelt(mv_10)
z.rm4 <- psmelt(mv_un10)
#plot the agglomerated SVs for reference of phylum abundance
z.rm4$Kingdom <- c('D: Misc')
z.rm4$Phylum <- c('P: Misc')
z.rm4$Class <- c('C: Misc')
z.rm4$Order <- c('O: Misc')
z.rm4$Family <- c('F: Misc')
z.rm4$Genus <- c('G: Misc')
z.rm4$Species <- c('Sp: Misc')
z.rm5 <- rbind(z.rm3, z.rm4)

## DESeq2 transformation plot - not what that transformation is for, but hey.
ggplot(z.rm5, aes(x=timepoint,y=Abundance,fill=UniqueSeq)) +
  facet_grid(Phylum+Genus~Reactor,scales='free',space='free') +                                  # scale/size control
  geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack") + #,  colour="black") +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 9)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +  # rotate
  scale_fill_manual(values = (z.col)) +  #improve SV differentiation
  ggtitle("MV Top SV's by UniqueSeq (Rel.Ab >2%, note variable scale)") + theme(legend.position="none")

# not much different, because youre wasting your time getting distracted. 

plot_bar(mv_D2,'timepoint',fill='UniqueSeq',title='MV DESeq2 Transformed') +
  facet_grid(Phylum~Reactor,scales='free_x',space='free_x') +
  scale_fill_manual(values=z.col) +
  theme(legend.position="none")

####




