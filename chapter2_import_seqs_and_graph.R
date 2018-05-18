#install.packages('tidyr')
library('tidyr')
library('phyloseq')
library('ggplot2')
library('vegan')
library('scales')



## PHYLO

# this map differs from other map by way of column names
ea_map <- 'env.txt'
ea_map <- import_qiime_sample_data(ea_map)
ea_map$Reactor <- c('R1','R1','R1','R1','R1','R1','R6','R6','R6','R6','R6')

# pre-moulded to fit our needs
ea_ab <-read.table('RDE_3up_otu.txt',sep='\t',header=TRUE,row.names = 1)
dim(ea_ab)
ea_ab <- otu_table( ea_ab, taxa_are_rows = TRUE)
ea_tax <- read.table('RDE_3up_tax.tsv',sep='\t',header=TRUE,row.names = 1)
ea_tax <- tax_table( as.matrix(ea_tax,rownames.force = TRUE) )
# ea_tax[,7] <- NULL ; ea_tax[,7]<-c(rownames(ea_tax))
# variations on:
sample_names(ea_map) <- c(colnames(ea_ab))
#  c('Domain','Phylum','Class','Order','Family','Genus','OTU')
ea<-merge_phyloseq(ea_ab,ea_tax,ea_map)
sample_data(ea)[,'SeqDepth']<-sample_sums(ea)
sample_data(ea)[,'Timepoint']<-sample_names(ea)
ea
# fix OPB54
tax_table(ea)[rownames(otu_table(subset_taxa(ea,Family=='OPB54'))),6] <-c('OPB54 G.')

# fix g
tax_table(ea)[rownames(otu_table(subset_taxa(ea,Genus=='g'))),6] <-c('uncultured G.')


## T R A N S F O R M   -   R E L . B U N
# rel.bun
ea_rb <- transform_sample_counts(ea, function(x)x/sum(x))
# hellinger
ea_h<-transform_sample_counts(ea,function(x) sqrt(x/sum(x)))

## S U B S E T
prun.ea_s.3.23 = genefilter_sample(ea, filterfun_sample(function(x) x >=10), A=3 )
ea_s.3.23 = prune_taxa(prun.ea_s.3.23, ea)
ea_sun.3.23 = prune_taxa(!prun.ea_s.3.23, ea)
ea_s.3.23_rb <- transform_sample_counts(ea_s.3.23, function(x)x/sum(x))
# arbitrary rel.abundance cutoff: >2%, x0.23 samples
prun.ea_s.01.23 = genefilter_sample(ea_rb, filterfun_sample(function(x) x >=0.01), A=2 )
ea_s.01.23 = prune_taxa(prun.ea_s.01.23, ea_rb)
# N O T E: negation of logical vector with prepending '!' !!
ea_sun.01.23 = prune_taxa(!prun.ea_s.01.23, ea_rb)
    # flatten all family for CCA plot
    #tax_table(ea_sun.01.23)[,'Family'] <- c('F:Misc & NA')
    #ea_s.sun <- merge_phyloseq(ea_s.01.23,ea_sun.01.23)

## record the taxa and OTUs
write.table(tax_table(ea),'ea_out_tax.txt',sep='\t')
write.table(otu_table(ea),'ea_out_otu.txt',sep='\t')




## S U M M A R Y   S T A T S  ##
# phylum breakdown across whole study, per-reactor, and maybe even per sample...
# not the same as the average % per sample, i.e. irrespective of library size.
#
##   M A K E   M E   I N T O   A   F U N C T I O N
### take phylo object, take sample variable: cut phylo by variable, by rank
## A L S O   A   F U N C T I O N   F O R   O T U   I D s

ea_phyglom <- tax_glom(ea,taxrank='Phylum')  # one glom instance preserves order in later sub-sets

# sum_phy inherentl ordered by SV abundance
sum_phy <- as.data.frame( round( (taxa_sums(ea_phyglom)/sum(taxa_sums(ea) ) *100), digits=4 ),  # grab rounded percentages
                          row.names = as.character(tax_table(ea_phyglom)[,2]) )                 # w. phyla as row names , cant set col.names?

# can illustrate taxa tables are intact based on sort(taxa_sums(z.isa/b/bes/ces)) etc.
z.1 <- subset_samples(ea_phyglom,Reactor=='R1' )
sum_phy[,2] <- round( (taxa_sums(z.1)/sum(taxa_sums(z.1) ) *100), digits=4 )  #keep NArm to keep comparable order between reactors
z.6 <- subset_samples(ea_phyglom,Reactor=='R6' )
sum_phy[,3] <- round( (taxa_sums(z.6)/sum(taxa_sums(z.6) ) *100), digits=4 )

colnames(sum_phy) <- c('TOT%','R1%','R6%') # cant seem to sort (order) properly either

## add number of genera per phylum (overall) #OR just the 3-5 phyla you'll actually talk about..
#length(get_taxa_unique(subset_taxa(ea,Phylum=='Firmicutes'),'Genus')) # this will overlook NAs etc
dim(tax_table(subset_taxa(ea,Phylum=='Firmicutes')))
dim(tax_table(subset_taxa(ea,Phylum=='Euryarchaeota')))
dim(tax_table(subset_taxa(ea,Phylum=='Bacteroidetes')))

View(sum_phy)

  # ## Stats per ea_01 / ea_un01
  # # # total ASV_orig, total sequences_orig
  # # ea_orig ; sum(sample_sums(ea_orig)) ; mean(sample_sums(ea_orig))
  # # total ASV, total sequences
  # ea ; sum(sample_sums(ea)) ; mean(sample_sums(ea))
  # # reads / sample
  # sum_phy 
  # # total un_ASV, total un_sequences
  # ea_01 ; sum(sample_sums(otu_table(ea)[,rownames(tax_table(ea_01))])) ; mean(sample_sums(otu_table(ea)[,rownames(tax_table(ea_01))]))
  # # total un_ASV, total un_sequences
  # ea_un01 ; sum(sample_sums(otu_table(ea)[,rownames(tax_table(ea_un01))])) ; mean(sample_sums(otu_table(ea)[,rownames(tax_table(ea_un01))]))
  # # Firm / Bact / Eury % range
  # sum_phy # see above
  # # #firm / #bact / #eury
  # dim(tax_table(subset_taxa(ea,Phylum=='Firmicutes')))
  # dim(tax_table(subset_taxa(ea,Phylum=='Euryarchaeota')))
  # dim(tax_table(subset_taxa(ea,Phylum=='Bacteroidetes')))
  # 
  # # NON- Firm / Bact / Eury % range
  # c(100, 100, 100, 100, 100) - colSums(sum_phy[1:3,])
  # get_taxa_unique(ea,taxonomic.rank = 'Phylum') ; length(get_taxa_unique(ea,taxonomic.rank = 'Phylum'))
  # #subtract this from 2300:
  # (dim(tax_table(subset_taxa(ea,Phylum=='Firmicutes'))) ) + ( dim(tax_table(subset_taxa(ea,Phylum=='Euryarchaeota'))) ) + ( dim(tax_table(subset_taxa(ea,Phylum=='Bacteroidetes'))) )
  # 2300 - 1662
  # # 




## P L O T T I N G
# open colour palette, appropriately sized
#z.col<-hue_pal()(length(get_taxa_unique(ea,'Family')))
#z.samp<-(length(get_taxa_unique(ea_s.01.23,'Genus')))
z.col42 = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3",
            "#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0",
            "#117878","#18A5A5","#3FE4E4","#3A5FCD","#4F94CD",             # #6CEAEA, #98F0F0
            "#117845","#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4",
            "#787811","#A5A518","#D2D21E","#EAAB6C","#EAEA6C","#F0F098","#F7F7C5",
            "#784511","#A55E18","#D2781E","#E4913F","#E4E43F","#F0C498",
            "#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7")
#z.col <- sample(z.col42,z.samp)
#or, by seq
z.seq <- z.col42[ seq(from=1, to=42, by=c(3)) ] ; z.seq2 <- c(z.seq, '#6CEAEA', '#98F0F0', "#D21E2C")
show_col(z.seq2)
z.seq2 <- c( z.col42[ seq(from=1, to=42, by=c(2)) ] )#, z.col42[ seq(from=2, to=42, by=c(2)) ] )
z.seq3 <- c( "#781156","#E43FAD","#114578","#3F91E4","#117878","#18A55E",
             "#3A5FCD","#D21E2C","#A5A518","#6CEAAB","#EAEA6C","#784511",
             "#E4913F","#781122","#E43F5B","#6CEAEA","#98F0F0")




pdf('ea-outversions_i.pdf',width=9, height=9)



## P R O C E S S   P L O T S 

# acids cut

# colours
# show_col(brewer.pal(n=15, name='RdYlGn'))
# "#D73027" , "#1A9850"

z.env1 <- read.table('ea_env_out_r1.csv', header=TRUE, sep='\t') #row.names=1,
z.env6 <- read.table('ea_env_out_r6.csv', header=TRUE, sep='\t') #row.names=1,
# 1
z.var1 <- melt(z.env1, id.vars=c(1) ) ;  z.var1[1:80,'measurement']<-c('Inhibitor') ; z.var1[81:160,'measurement']<-c('Indicator') ; z.var1[,'Reactor'] <- c('R1')
z.var6 <- melt(z.env6, id.vars=c(1) ) ;  z.var6[1:80,'measurement']<-c('Inhibitor') ; z.var6[81:160,'measurement']<-c('Indicator') ; z.var6[,'Reactor'] <- c('R6')
z.var <- rbind( z.var1 , z.var6 )

z.var.hib <- z.var[z.var$'measurement'=='Inhibitor',] ; colnames(z.var.hib)[3] <- c('mg.L')
ggplot(z.var.hib, aes(Week, mg.L, group=Reactor)) +
  facet_grid(~variable,scales='free',space='free') +                                  # scale/size control
  geom_line(aes(color=Reactor), size=0.8, data = na.omit(z.var.hib)) +
  geom_point(aes(color=Reactor, shape=Reactor) ) +
  scale_colour_manual(values= c( "#D73027" , "#1A9850" )) +
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 12)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 0)) +  # rotate
  theme(panel.background = element_rect(fill='white', colour='grey70')) +
  theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
  theme(panel.grid.minor.y = element_line(colour='grey85', size = 0.2)) +                               # horiz lines 2
  theme(panel.grid.major.x = element_line(colour='grey75', size = 0.2)) +                               # vert lines
  theme(panel.grid.minor.x = element_line(colour='grey85', size = 0.2)) #+                               # vert lines
#  geom_hline(aes(yintercept = 7500, colour = 'red'), z.var.hib) +                                      # cant as plots on both reactors 
#ggtitle("Inhibitory variables during operation")


z.var.dic <- z.var[z.var$'measurement'=='Indicator',] ; colnames(z.var.dic)[3] <- c('Ratio')
ggplot(z.var.dic, aes(Week, Ratio, group=Reactor)) +
  facet_grid(~variable,scales='free',space='free') +                                  # scale/size control
  geom_line(aes(color=Reactor), size=0.8, data = na.omit(z.var.dic)) +
  geom_point(aes(color=Reactor, shape=Reactor) ) +
  scale_colour_manual(values= c( "#D73027" , "#1A9850" )) +
  scale_y_continuous(breaks = seq(0, 1.6, by = 0.2)) +
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 12)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 0)) +  # rotate
  theme(panel.background = element_rect(fill='white', colour='grey70')) +
  theme(panel.grid.major.y = element_line(colour='grey75', size = 0.2)) +                               # horiz lines 1
  theme(panel.grid.minor.y = element_line(colour='grey85', size = 0.2)) +                               # horiz lines 2
  theme(panel.grid.major.x = element_line(colour='grey75', size = 0.2)) +                               # vert lines
  theme(panel.grid.minor.x = element_line(colour='grey85', size = 0.2)) #+                               # vert lines
# ggtitle("Indicator variables during operation")




# diversity metrics: 
plot_richness(ea, x = "Reactor", color = "SeqDepth") #+ geom_boxplot()


## A B U N D A N C E ##
plot_bar(ea,'Timepoint',fill='Class',title='ea-SILVA Absolute Read Abundance')
#plot_bar(ea_rb,'Timepoint',fill='Class',title='ea-SILVA Relative Read Abundance') + facet_grid(~State,scales='free_x',space='free_x')
plot_bar(ea_s.3.23_rb ,fill='Class',title='ea-SILVA.3.23 Relative Read Abundance') + facet_grid(~Reactor,scales='free_x',space='free_x')

z.rm3 <- psmelt(ea_s.01.23) #(subset_samples(ea_s.01.23,Reactor=='G'))
z.rm4 <- psmelt(ea_sun.01.23) #(subset_samples(ea_sun.01.23,Reactor=='G'))
z.rm4[,25] <- c('P:Other')
z.rm4[,26] <- c('C:Other')
z.rm4[,27] <- c('O:Other')
z.rm4[,28] <- c('F:Other')
z.rm4[,29] <- c('G:Other')
z.rm5 <- rbind(z.rm3, z.rm4)

#show_col(z.seq)
ggplot(z.rm5, aes(x=Timepoint,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  scale_fill_manual(values = (z.seq3), na.value='grey40') +
  facet_grid(Phylum~Reactor,scales='free',space='free') +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14)) +   # bigger
  theme(panel.grid.minor.y = element_line(colour='grey85')) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.background = element_rect(fill='white')) +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270))
#  scale_color_manual(values = (z.seq), na.value='grey40') +
#  ggtitle(">1% Relative Abundance, 93% taxonomic ID, 100% identity clustering") +
  #  theme(legend.position="none")


  ## visualise for minors
      z.seq2 <- c( z.col42[ seq(from=1, to=42, by=c(2)) ] , z.col42[ seq(from=2, to=42, by=c(2)) ] )

          # new 'Other' grouping 
          z.rm2 <- (ea_sun.01.23)
          prun.clos = genefilter_sample(z.rm2, filterfun_sample(function(x) x >=0.005), A=1 )
          z.rm3 = prune_taxa(prun.clos, z.rm2) ; z.rm3 <- psmelt(z.rm3)
          z.rm4 = prune_taxa(!prun.clos, z.rm2)
          # View(z.rm4)
          z.rm4 <- psmelt(z.rm4) ; z.rm4[,29] <- c('G:Other, <0.5%') ; z.rm4[,28] <- c('F:Other, <0.5%') ; z.rm4[,27] <- c('O:Other, <0.5%') 
          
          z.rm6 <- rbind(z.rm3, z.rm4)
          ggplot(z.rm6, aes(x=Timepoint,y=Abundance,fill=Genus)) + 
            geom_bar(aes(fill=Genus), stat ="identity", position="stack", color='black') +
#            scale_color_manual(values = (z.seq2), na.value='grey40') +
#            scale_fill_manual(values = (z.seq2), na.value='grey40') +
            facet_grid(Order+Family~Reactor,scales='free',space='free') + 
            theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) + 
            theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +   # bigger
            ggtitle("Minor (<0.5%) OTUs") + 
            theme(legend.position="none") 


## Ladies and Gentlemen; THE ARCHAEA!
z.rmA <- psmelt(subset_taxa(ea_rb,Kingdom=='Archaea'))
ggplot(z.rmA, aes(x=Timepoint,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Family~Reactor,scales='free_x',space='free_x') +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +   # bigger
  ggtitle("ea SLV Archaeal OTUs by Family/Orders")
## ...   *clap*,     *clap*.


## O R D I N A T E
plot_ordination(ea,ordinate(ea,'DCA'),color='SeqDepth',shape='Reactor',title='CCA of ea-RB',label='Timepoint',type='samples')



## import   L E F S E   output for plotting
ea_lefse <- as.data.frame(read.table('ea_lefse_edit.csv',sep='\t',header=TRUE,row.names = 1))
#colnames(ea_lefse) <- c('taxa','KW.RST','Setup','LDA','p.val') # already in file, ported with header/row.names above

## table com, binned with LefSe output in excel to show the relevant stuff: many higher clades are merely
## representations of their OTU, but at a higher confidence due to bootstrapping: these do not contribute
## but do clutter - smooth out!
ea_lefse <- ea_lefse[with(ea_lefse, order(Reactor, -KW.RST)), ]    # sort by A, then by B!  *** 
ea_lefse$Taxa <- factor(ea_lefse$Taxa, levels = ea_lefse$Taxa)     # 'fix' order to factors *** 

ggplot(ea_lefse, aes(x=Taxa, y=KW.RST, fill= Reactor, shape=Reactor)) +
  theme(axis.text = element_text(size=13), legend.text=element_text(size=13), axis.title = element_text(size=13)) +
  theme(strip.text.x = element_text(size = 15)) + 
  facet_grid(~Reactor, space = 'free', drop = ) +
  geom_point(aes(size=LDA)) +
  coord_flip() +
  guides(fill=guide_legend(title.theme=element_text(size=21, face= 'bold', colour='green'))) +
#  labs(title='Taxa Associations with Reactor Setup') +
  labs(y='Magnitude of Difference (H score)') +
  scale_size_continuous(name='LDA Effect Size (log10)') +
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 12)) +   # bigger
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 0)) +  # rotate
  theme(panel.background = element_rect(fill=NA, colour=NA)) +                              # rm backgrnd
  theme(panel.grid.major.x = element_line(colour='grey80')) +                               # horiz lines 1
  theme(panel.grid.major.y = element_line(colour='grey80')) +                               # horiz lines 1
  scale_shape_manual(values = c(22, 21, 24)) +
  scale_fill_manual(values = c(brewer.pal(9,'YlOrRd')[8],brewer.pal(3,'Greens')[3],brewer.pal(3,'YlGnBu')[3])) #+



dev.off()

