###
### mv_images 
###

# 
# DO:

# IMAGES
  VFA profiles, FOSTAC, Biogas?
# #Diversity - 	change colour - green?
  Rarefaction
  CCA ... of what though?
# ?plot_network
# 
# 
# EQUATION
# Hellinger
# 
# 
# APPENDICES:
# A: DNA extraction, recipes etc.
# B: Sequence Provenance
# C: ASV tables?

# 
# ==========
# 
# intro
#   biofuel: RED, transport,fuel etc. 3rd generation weightings.
#   AD: include flann o brien. 4 step rpogam. 
#   mechanisms: hydrolysis, fermentation, acidogenesis

## C O L O U R S
z.col <- hue_pal()( nrow(tax_table(mv)) )  ;  z.col <- sample(z.col, nrow(tax_table(mv)) )
  z.col1 <- hue_pal(l=55, c=80)( nrow(tax_table(mv_01)) )  ;  z.col1 <- sample(z.col1, nrow(tax_table(mv_01)) )
  z.col2 <- hue_pal(l=50, c=85)( nrow(tax_table(mv_un01)) )  ;  z.col2 <- sample(z.col2, nrow(tax_table(mv_un01)) )
    z.col42 = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3",
                "#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0",
                "#117878","#18A5A5","#3FE4E4","#3A5FCD","#4F94CD",             # #6CEAEA, #98F0F0
                "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4",
                "#787811","#A5A518","#D2D21E","#EAAB6C","#EAEA6C","#F0F098","#F7F7C5",
                "#784511","#A55E18","#D2781E","#E4913F","#E4E43F","#F0C498",
                "#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7")
    show_col(z.col42)
      z.col35 <- z.col42[-c(6,12,17,23,30,36,42)]
      z.col35 <-c(z.col35,colours()[584:588],colours()[593:598])
      show_col(z.col35)
        z.col57 <-c(z.col42,colours()[584:598])    # show_col(colours()[584:598])
        z.col57 <- sample(z.col57,57)
  


# pdf('mv_images.pdf',width=10, height=6)

## V F A   +   G A S
#overly precious method
z.env <- as.data.frame(sample_data(mv)[,c(2,5,7,8,10,11,15:23)], row.names = rownames(sample_data(mv)) )
z.vfa <- melt(z.env[,-c(5,6)], id.vars=c(1,2) )   #exclude non-target data, assign 1+2 as variables
z.vfa[1:66,5] <- c('Gas') ; z.vfa[67:363,5] <- c('Acid')   #ugly but works
colnames(z.vfa)[4] <- c('mg.L')
g.vfa <- ggplot(z.vfa, aes(timepoint, mg.L, fill=variable)) +
  facet_grid(V5~Reactor,scales='free',space='free') +                                  # scale/size control
  geom_bar(aes(fill=variable), stat ="identity", position="stack",  colour="black") +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 0)) +  # rotate
  scale_fill_manual(values = c( rev(PiYG(11)[c(2,10)]), rev(YlGnBu(12)[2:10] ) ))  +
  theme(panel.background = element_rect(fill='white', colour='grey70')) +  # 
  theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
  theme(panel.grid.minor.y = element_line(colour='grey85', size = 0.1)) +                               # horiz lines 2
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
  ggtitle("Acid and gas profiles during operation")


### P R O C E S S   V A R I A B L E S
# similarly: (see also:  "plot.ts(z.env[,3:6], na.omit())" )
z.var <- melt(z.env[,c(1,2,5,6)], id.vars=c(1,2) )   #exclude non-target data, assign 1+2 as variables
g.var <- ggplot(z.var, aes(timepoint, value, group=variable)) +
  facet_grid(variable~Reactor,scales='free',space='free') +                                  # scale/size control
  geom_line(aes(color=variable), data = na.omit(z.var)) +
  geom_point(aes(color=variable) ) +
  scale_colour_brewer(palette = "Set1") +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 0)) +  # rotate
  theme(panel.background = element_rect(fill='white', colour='grey70')) +
  theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
  theme(panel.grid.minor.y = element_line(colour='grey85', size = 0.1)) +                               # horiz lines 2
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
  ggtitle("Process variables during operation")

      
## L E N G T H 
g.hist <- ( hist(as.vector(apply( tax_table(mv)[,8], 1, function(x) nchar(x) )),main='histogram of all SV lengths (<350 trimmed)')  )

## R I C H N E $ $
g.div <- plot_richness(mv, x = "Reactor",
                      color = "SeqDepth",
                      shape = 'Reactor',
                      measures=c('Observed',
                                 'Shannon',
                                 'InvSimpson')
) +
  geom_point(size=3) +
  scale_shape_manual(values = c(17, 15, 16, 18)) +
  scale_fill_distiller(palette='Spectral')            # cant change colour

## T O T A L   R E A D S
g.reads <- plot_bar(mv,'timepoint',fill='Phylum',title='Total Processed Reads: 6,111,977') +
  facet_grid(~Reactor,scales='free_x',space='free_x') +
  theme(panel.background = element_rect(fill='white', colour='grey70')) +
  theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                   # horiz lines 1
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
  scale_fill_manual(values = z.col42)

## D C A   -  S e q d e p t h 
g.dca <- plot_ordination(mv,ordinate(mv,'DCA'),  #, formula = (~H2)
                         color='SeqDepth',
                         shape='Reactor',
                         title='DCA of Samples in Hellinger Distance, Emphasising Library Size (SeqDepth)',
                         type='samples') + 
          theme(panel.background = element_rect(fill='white', colour='grey20')) +
          theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +         # horiz lines 1
          theme(panel.grid.major.x = element_line(colour='grey75', size = 0.2)) +          # horiz lines 2
          # text option for label, offset for lable
          geom_text(label = rownames(sample_data(mv)), size=3, nudge_x = 0.2, check_overlap = TRUE, color='black') +
          scale_color_viridis() +
          scale_fill_viridis() +
          scale_shape_manual(values = c(17, 15, 16, 18)) +
          geom_point(size=6) 



##
## R E L A T I V E   A B U N D A N C E S
##
z.rm3 <- psmelt(mv_01)

z.rm4 <- psmelt(mv_001)
z.rm5 <- psmelt(mv_un001) ; z.rm5[,29:35] <- c('Misc. <0.1%')
z.rm6 <- rbind(z.rm4, z.rm5)

###    -   M A J O R
g.rama <- ggplot(z.rm3, aes(x=timepoint,y=Abundance,fill=UniqueSeq)) +
  #  facet_grid(Phylum+Class+Order+Family+Genus~Reactor,scales='free_x',space='free') +                                  # scale/size control
  facet_grid(Phylum+Genus~Reactor,scales='free',space='free') +                                  # scale/size control
  geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack",  colour="black") +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 9)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +  # rotate
  scale_fill_manual(values = (c(sort(z.col1), (z.col2)))) +  #improve SV differentiation #c(z.col42[1:37]) , 
  theme(panel.background = element_rect(fill='white', colour='grey70')) +
  theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
#  theme(panel.grid.minor.y = element_line(colour='grey85', size = 0.1)) +                               # horiz lines 2
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
  ggtitle("Major ASVs (Rel.Ab >1%, note variable scale)") + theme(legend.position="none")

###    -   M I N O R
g.rami <- ggplot(z.rm6, aes(x=timepoint,y=Abundance,fill=UniqueSeq)) +
  #  facet_grid(Phylum+Class+Order+Family+Genus~Reactor,scales='free_x',space='free') +                                  # scale/size control
  facet_grid(Phylum~Reactor,scales='free',space='free') +                                  # scale/size control
  geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack") +                      #,  colour="black") +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 9)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +  # rotate
  scale_fill_manual(values = (c(sort(z.col1), (z.col2)))) +  #improve SV differentiation #c(z.col42[1:37]) , 
  theme(panel.background = element_rect(fill='white', colour='grey70')) +
  theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
  ggtitle("Minor ASVs (Rel.Ab <1%, note variable scale)") + theme(legend.position="none")

###  A R C H A E A
z.colA <- sample(z.col57, nrow(tax_table(subset_taxa(mv_rb,Kingdom=='Archaea'))) )
z.rmA <- psmelt(subset_taxa(mv_rb,Kingdom=='Archaea'))
g.arc <- ggplot(z.rmA, aes(x=timepoint,y=Abundance,fill=z.colA)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Order+Genus~Reactor,scales='free',space='free_x') +
  scale_fill_manual(values = (sample(z.col57, length(z.colA)))) +  #improve SV differentiation #c(z.col42[1:37]) , 
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 11)) +      # biggerness
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +     # rotate
  theme(panel.background = element_rect(fill='white', colour='grey70')) +
  theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
  ggtitle("All Archaeal ASVs")  + theme(legend.position="none")
show_col(z.col57)





## 
####    D E S E Q 2  
##

##      D E S E Q 2  -    I S   v.   E s
mv_in <- phyloseq_to_deseq2(mv, ~situ)
# prefiltering to remove outgrageous LFC from near absent taxa.
keep <- rowSums(counts(mv_in)) >= 1000
mv_in <- mv_in[keep,]
# Wald, relevel
mv_in$situ <-relevel(mv_in$situ,ref='In-Situ')
mv_in_s<-DESeq(mv_in) #fitType='local',minReplicatesForReplace =5)
# cut copy preach
mv_in_r1_res <- results(mv_in_s, cooksCutoff=TRUE, alpha=0.01, lfcThreshold=1) #, contrast=c('setup','CES','IS'))   #
mv_in_r1_sig <- mv_in_r1_res[which(mv_in_r1_res$padj<0.01),]
mv_in_r1_sig <- cbind(as(mv_in_r1_sig,'data.frame'),as(tax_table(mv)[rownames(mv_in_r1_sig),],'matrix'))
g.d2ises <-  ggplot(mv_in_r1_sig, aes(x=Genus, y=log2FoldChange, color= Phylum)) +
      geom_hline(yintercept=0, color='grey75') +
      theme(axis.text = element_text(size=12)) +
      labs(title='DA: in situ (L) v. ex situ (R)   (p/padj=0.01/0.01)') +
      coord_flip() + scale_size_continuous(name='root sum abundance') +
      scale_colour_manual(values = (z.col57[seq(2,58,4)])) +  #custom fill
      geom_jitter(alpha=0.6, aes(size=( sqrt(33*baseMean) ) ), width=0.1, height=0.1) +  # see also geom_point
      theme(panel.background = element_rect(fill='white', colour='grey70')) +
      theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
      theme(panel.grid.major.x = element_line(colour='grey75', size = 0.2)) +                                             # rm vert lines
      geom_point(color='grey35',size=0.1) 

# catch NAs
z.rm <- mv_in_r1_sig[which(is.na(mv_in_r1_sig$Genus)),]
g.d2_isesNA <- ggplot(z.rm, aes(x=Phylum, y=log2FoldChange, color= Family)) +
  geom_hline(yintercept=0, color='grey75') +
  theme(axis.text = element_text(size=12)) +
  labs(title='DA: in situ (L) v. ex situ (R) NA\'s   (p/padj=0.01/0.01)') +
  coord_flip() + scale_size_continuous(name=' total abundance') +
  scale_colour_manual(values = (sample(z.col,22)), na.value='red') +  #custom fill, custom NAs
  geom_jitter(alpha=0.6, aes(size=((15*baseMean)) ), width=0.1, height=0.1) +  # see also geom_point
  geom_point(color='grey35',size=0.8) #+
g.d2_isesNA


##   S U B S E T   T O   E X   S I T U
#
z.rm <- subset_samples(mv, situ=='Ex-Situ')
prunB <- genefilter_sample(z.rm, filterfun_sample(function(x) x >0)) ; z.rm <- prune_taxa(prunB, z.rm)
z.rb <- transform_sample_counts(z.rm, function(x)x/sum(x))

prunB = genefilter_sample(z.rb, filterfun_sample(function(x) x >=0.01), A=2)  
z.02 = prune_taxa(prunB, z.rb)
z.un02 = prune_taxa(!prunB, z.rb)    # N O T E: negation of logical vector with prepending '!' !!

z.rm3 <- psmelt(z.02) ; z.rm4 <- psmelt(z.un02)
z.rm4[,29:35] <- c('Misc <1%')  # subset further?
z.rm5 <- rbind(z.rm3, z.rm4)

    ## EX SITU REL ABUND
g.raes  <-  ggplot(z.rm5, aes(x=timepoint,y=Abundance,fill=UniqueSeq)) +
      facet_grid(Phylum+Genus~Reactor,scales='free', space='free') +                      # scale/size control
      geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack",  colour="black") +
      theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 9)) +   # bigger
      theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +  # rotate
      scale_fill_manual(values = (z.col)) +  #improve SV differentiation
      ggtitle("Major ex situ ASVs (Rel.Ab >1%, note variable scale)") + theme(legend.position="none") +
      theme(panel.background = element_rect(fill='white', colour='grey70')) +
      theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +            # horiz lines 1
#      theme(panel.grid.minor.y = element_line(colour='grey85', size = 0.1)) +            # horiz lines 2
      theme(panel.grid.major.x = element_blank())                                        # rm vert lines
      

##  D E S E Q 2 :   C E S   v.   B E S 
##
mv_es <- phyloseq_to_deseq2(z.rm, ~Reactor)
keep <- rowSums(counts(mv_es)) >= 1000  # cut out smallies
mv_es <- mv_es[keep,]
# relevel, Wald test, 0.01/0.01
mv_es$Reactor <-relevel(mv_es$Reactor,ref='3.BES')
mv_es_s<-DESeq(mv_es) #fitType='local',minReplicatesForReplace =5)
mv_es_r1_res <- results(mv_es_s, lfcThreshold=1, alpha=0.01)  
mv_es_r1_sig <- mv_es_r1_res[which(mv_es_r1_res$padj<0.01),]
mv_es_r1_sig <- cbind(as(mv_es_r1_sig,'data.frame'),as(tax_table(mv)[rownames(mv_es_r1_sig),],'matrix'))

    ## MV BES v. CES 
g.d2es  <-  ggplot(mv_es_r1_sig, aes(x=Genus, y=log2FoldChange, color= Phylum)) +
      geom_hline(yintercept=0, color='grey75') +
      theme(axis.text = element_text(size=12)) +
      labs(title='DA ASVs: BES v. CES   (p/padj=0.01/0.01)') +
      coord_flip() + scale_size_continuous(name=' total abundance') +
      scale_colour_manual(values = (z.col57[seq(0,58,4)]), na.value='red') +  #custom fill
      geom_jitter(alpha=0.6, aes(size=((15*baseMean)) ), width=0.1, height=0.1) +  # , shape=Kingdom
      geom_point(color='grey35',size=0.5) +
      theme(panel.background = element_rect(fill='white', colour='grey70')) +
      theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +  # horiz lines 1
      theme(panel.grid.major.x = element_line(colour='grey75', size = 0.2))    # vert lines
  
    ## Catch and explore the NAs in another plot: not really anything to see here! :)
    z.rm <- mv_es_r1_sig[which(is.na(mv_es_r1_sig$Genus)),]
g.d2esNA  <-    ggplot(z.rm, aes(x=Phylum, y=log2FoldChange, color= Order)) +
      geom_hline(yintercept=0, color='grey75') +
      theme(axis.text = element_text(size=12)) +
      labs(title='DA \'NA\' ASVs: BES v. CES') +
      coord_flip() + scale_size_continuous(name=' total abundance') +
      scale_colour_manual(values = (z.col57[seq(0,58,4)]), na.value='red') +  #custom fill, custom NAs
      geom_jitter(alpha=0.6, aes(size=((15*baseMean)) ), width=0.1, height=0.1) +  # see also geom_point
      geom_point(color='grey35',size=0.8) +
      theme(panel.background = element_rect(fill='white', colour='grey70')) +
      theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +  # horiz lines 1
      theme(panel.grid.major.x = element_line(colour='grey75', size = 0.2))    # vert lines

    # ES sig rel abunds
    z.rm6 <- psmelt(prune_taxa(rownames(mv_es_r1_sig), z.rb))   #lose too much : [rownames(mv_es_r1_sig) %in% rownames( tax_table(mv_01))] 
g.d2ra  <-  ggplot(z.rm6, aes(x=timepoint,y=Abundance,fill=UniqueSeq)) +
      facet_grid(Phylum+Order~Reactor,scales='free',space='free') +                                  # scale/size control
      geom_bar(aes(fill=UniqueSeq), stat ="identity", position="stack",  colour="black") +
      theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 9)) +   # bigger
      theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +  # rotate
      scale_fill_manual(values = (z.col)) +
      ggtitle("DA ASVs: BES v. CES (Rel.Ab >1%, note variable scale)") + theme(legend.position="none") +
      theme(panel.background = element_rect(fill='white', colour='grey70')) +
      theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
#      theme(panel.grid.minor.y = element_line(colour='grey85', size = 0.1)) +                               # horiz lines 2
      theme(panel.grid.major.x = element_blank())                                             # rm vert lines
   
 
###
### H E A T M A P L Y
###
    ## IMPORT COR & PP MATRICES
##BLIP
#    mv_cor <- read.table('C:/Users/Blip the Tall/Downloads/mv_sparcc/mv_cor.txt',sep='\t',header=TRUE,row.names=1)
#    mv_pp <- read.table('C:/Users/Blip the Tall/Downloads/mv_sparcc/mv_two_sided_pp.txt',sep='\t',header=TRUE,row.names=1)
##JFG
#    mv_cor <- read.table('/home/jfg/BioApps/sparcc/mv_sparcc/mv_cor.txt',sep='\t',header=TRUE,row.names=1)
#    mv_pp <- read.table('/home/jfg/BioApps/sparcc/mv_sparcc/mv_two_sided_pp.txt',sep='\t',header=TRUE,row.names=1)

    # sort phylo object and then pass to decendants 
    # mv_sort <- merge_phyloseq(otu_table(mv),tax_table(mv),sample_data(mv)) ## if necessary to remove tree
    mv_sort <- mv
    tax_table(mv_sort) <- tax_table(mv_sort)[with(as.data.frame(tax_table(mv_sort)), order(Kingdom,Phylum,Class,Order,Family,Genus)),]
    otu_table(mv_sort) <- otu_table(mv_sort)[,with(as.data.frame(tax_table(mv_sort)), order(Kingdom,Phylum,Class,Order,Family,Genus))]
    tax_table(mv_sort)[45:55,1:5]
    
    # propagate subset of 2300x2300
    mv_sort_id <- rownames(tax_table(mv_sort)) ; head(mv_sort_id)
    mv_cor2 <- mv_cor[mv_sort_id,mv_sort_id] ;   mv_cor2[1:7,1:7]
    mv_pp2 <- mv_pp[mv_sort_id,mv_sort_id] ;     mv_pp2[1:10,1:10]
    
    ## NUKE CORR via PP >0.1  - doing early avoids confounding later due to autocorelation=1
    z.pNA <- apply(mv_pp2,c(1,2),function(x)x>0.01)    # logical matrix of values >0.01
    # consider a more stringent value to reduce size of set
    mv_cor2[z.pNA] <- c(0)                              # index cor_10 via log.mat and apply NA
    
    z.cor1 <- mv_cor2[rownames(tax_table(mv_01)), rownames(mv_in_r1_sig)]
    rownames(z.cor1) <- paste( sub('mvsv', '', rownames(tax_table(mv_01)) ) , tax_table(mv_01)[,6] ) 
    colnames(z.cor1) <- paste( sub('mvsv', '', rownames(mv_in_r1_sig))  , tax_table(mv)[rownames(mv_in_r1_sig),6] ) #, 1:nrow(cor_10) ) 
    # NOW remove non-informative rows, helping x-axis legibility:
    z.1 <- Filter(function(x) all(x < 0.75 & x > -0.75), z.cor1) ; dim(z.1)         # identify the shit ones
    z.2 <- !(colnames(z.cor1) %in% colnames(z.1))
    z.co1a <- z.cor1[  , z.2]
g.heat  <-  heatmaply(z.co1a, col=c("#440154FF","#21908CFF","#21908CFF","#27AD81FF","#21908CFF","#21908CFF","#FDE725FF"), xlab='DA ASVs (IS v. ES)', ylab='Major Abundant ASVs (>1%)') #, file='mv_co1a_stretch.html')
    

    #     
    # ## P L O T   N E T W O R K
    #     plot_net(mv, type='samples', color='Acetic.acid', shape='Reactor', maxdist = 0.3) +
    #       geom_point(shape=c(17, 15, 16, 18))
    

# pdf('mv_graphics.pdf',width=10,height=12)
#ls(pattern='g.')
g.vfa
g.var
g.arc  #!
g.d2es    #57
g.d2esNA  #57
g.d2ises  #57  
g.d2_isesNA  
g.d2ra   #cut?
g.dca
g.div
#g.heat
plot(g.hist)
g.raes  # quash minors
g.rama
g.rami   
g.reads
dev.off()

## E X I T 
#    dev.off()
z.rm <- ls(pattern='g.') ; save(z.rm,file='mv.graphical_calls.RData')
svg((z.rm))

