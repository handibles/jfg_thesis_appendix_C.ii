## dw_mk.ii images
## S T A G E   D R E S S I N G
# for either Windows or BioUbuntu etc., change dir to contain the ####supplementary files###XYZ#######:
#setwd("/home/user/dw")
#setwd("C:/Users/user/dw")
## packages as necessary
#install.packages('tidyr','phyloseq','ggplot2','DESeq2','ancom.R','vegan','RColorBrewer','scales')

# # wave at the cavalry
# source('https://bioconductor.org/biocLite.R')
# biocLite()

## queue lights
library('tidyr')
library('phyloseq')
library('ggplot2')
library('DESeq2')
library('vegan')
library('RColorBrewer')
library('scales')
## 3, 2, 1:




##                                                       ##
#                                                         #
#                                                         #
##    Import, shape, analyse and finally export data:    ##
#                                                         #
#                                                         #
##                                                       ##


## M E T A D A T A
# mapfile contains the metadata: here take oportunity to add some more 
dw_slv_map <- '../dw_map.txt'
dw_slv_map <- import_qiime_sample_data(dw_slv_map)
## additional data
dw_slv_map$Treat <- c('silage','silage','silage','silage','silage','silage','silage','silage','silage w.TE','silage w.TE','slurry-silage','slurry-silage','slurry-silage','slurry-silage')
dw_slv_map$TE <- c('low-TE','low-TE','low-TE','low-TE','low-TE','low-TE','low-TE','low-TE','hi-TE','hi-TE','hi-TE','hi-TE','hi-TE','hi-TE')
# Beff values calculcated from BMP values given in Wall et al., 2014
dw_slv_map$Beff <- dw_slv_map$SMY/c(400,400,400,400,400,400,400,400,400,400,348,348,348,348)
##RE-MAKE ENTIRE MAP FILE HERE?

## I M P O R T   S I L V A   O T U

### REMAKE OTU FILE? SPLIT/ APPEND ACC#?
## import from LEfSe output, and convert to phyloseq object
#dw_slv_tax<-(separate(data = dw_slv_input, col = 'taxonomy', into = c(Domain","Phylum","Class","Order","Family","Genus","OTU"), sep = ";"))[,16:21]
# file tidied before addition to R
dw_slv_input<-read.table('dw_splibbed---ssu---fingerprint----Total---sim_95---tax_silva---td_20.csv',sep='\t',header=TRUE,row.names = 1)
dw_slv_ab <- dw_slv_input[,1:14]
# cut out NO RELATIVES (row 345), as cannot contribute to analysis
dw_slv_ab <- dw_slv_ab[1:344,]
dw_slv_ab <- otu_table(dw_slv_ab,taxa_are_rows = TRUE)
# cut out NO RELATIVES (row 345), as cannot contribute to analysis
dw_slv_tax<-(dw_slv_input[1:344,15:21])
dw_slv_tax<-(dw_slv_input[,15:21])
dw_slv_tax <- as.matrix(dw_slv_tax,rownames.force = TRUE)
dw_slv_tax<-tax_table(dw_slv_tax)
# set row and column names
row.names(dw_slv_tax) <- row.names(dw_slv_input)
colnames(dw_slv_tax) <- c('Domain','Phylum','Class','Order','Family','Genus','OTU')
dw_slv<-merge_phyloseq(dw_slv_ab,dw_slv_tax,dw_slv_map)
# add sequencing depth metadata via phyloseq
sample_data(dw_slv)[,'SeqDepth']<-sample_sums(dw_slv)

# view phyloseq object
dw_slv


## T I D Y   U P   O T U s
# tidy MBA03 OTU
tax_table(dw_slv)['silva4_253',] <-c('Bacteria','Firmicutes','Clostridia','MBA03','MBA03 F.','MBA03 G.','silva4_253')
# tidy vadinBC27 OTU
tax_table(dw_slv)['silva4_051',] <-c('Bacteria','Bacteroidetes','Bacteroidia','Bacteroidales','Rikenellaceae','vadinBC27 group','silva4_051')

# Abbreviate accession from Ruminococcaceae UCG-012 as clunks display
tax_table(dw_slv)['silva4_236',] <-c('Bacteria','Firmicutes','Clostridia','Clostridiales','Ruminococcaceae','Rum.\'aceae UCG-012','silva4_236')

# Abbreviate accession from Christensenellaceae R-7 as clunks display
tax_table(dw_slv)['silva4_113',] <-c('Bacteria','Firmicutes','Clostridia','Clostridiales','Christensenellaceae','Chri.\'aceae R-7 group','silva4_236')



## T R A N S F O R M   -   R E L A T I V E    A B U N D A N C E
# rel.bun
dw_slv_rb <- transform_sample_counts(dw_slv, function(x)x/sum(x))

## S U B S E T
## subset OTUs for display purposes:
## agglomerate OTUs with which have less than 2% relative abundance as they cannot be displayed

# arbitrary rel.abundance cutoff: >2%, in 23% of samples
prun.dw_s.02.23 = genefilter_sample(dw_slv_rb, filterfun_sample(function(x) x >=0.02), A=2) #0.23*nsamples(dw_slv_rb))
dw_s.02.23 = prune_taxa(prun.dw_s.02.23, dw_slv_rb)
dw_sun.02.23 = prune_taxa(!prun.dw_s.02.23, dw_slv_rb)
# flatten all pruned families for plotting
tax_table(dw_sun.02.23)[,'Family'] <- c('F:Misc & NA')
dw_s.sun <- merge_phyloseq(dw_s.02.23,dw_sun.02.23)



## P L O T T I N G

pdf('dw_trace_element_supplements.pdf',width=12, height=9)

## P R O C E S S - note change in order of samples to Timepoint
z.rm <- sample_data(dw_slv)[with(sample_data(dw_slv),order(Timepoint)),]
z.rm[,17:19]<-NULL
z.rm[,1:5]<-NULL
z.rm[,'TE']<-NULL
z.rm[,'w_teor']<-NULL
z.rm[,'NH3']<-NULL
z.rm[,'CH4_pc']<-NULL
z.rm[,'Timepoint']<-NULL
z.rm[,'Week']<-NULL

plot.ts(z.rm[1:10,1:8],main='Process for G samples')
plot.ts(z.rm[10:14,1:8],main='Process for SG samples')


## S U M M A R Y   S T A T S  ##
# diversity metrics: 
plot_richness(dw_slv, x = "Treat", color = "SeqDepth", shape = 'Reactor')

# corrleation analysis of sequence depth by community structure:
plot_ordination(dw_slv,ordinate(dw_slv,'CCA'),color='SeqDepth',shape='Treat',title='CA of dw showing sequencing depth',label='Timepoint',type='biplot')

# raw read abundances, coloured by Phylum
plot_bar(dw_slv,'Timepoint',fill='Phylum',title='Sample Absolute Read Abundance') + facet_grid(~Treat,scales='free_x',space='free_x')
# relative abundances, coloured by Order
plot_bar(dw_s.sun ,'Timepoint',fill='Family',title='Relative Read Abundance (taxa <2% grouped)') + facet_grid(~Treat,scales='free_x',space='free_x')


## A B U N D A N C E S
f2.a <- psmelt(dw_s.02.23)
f2.b <- psmelt(dw_sun.02.23)
#plot the agglomerated OTUs for reference of phylum abundance
## pool all non->2% to dummy phylum, append to bottom of melted frame 
f2.b[,27] <- c('Phyla <2%')
f2.b[,28] <- c('Classes <2%')
f2.b[,29] <- c('Orders <2%')
f2.b[,30] <- c('Families <2%')
f2.b[,31] <- c('Genera <2%')
f2.c <- rbind(f2.a, f2.b)
# make abundance a percentage rather than a proportion
f2.d <- f2.c ; f2.d$Abundance <- (f2.d$Abundance)*100

# improve colours 
z.col <- colors(distinct=TRUE)[1:135] ; z.col <- c(z.col,colors(distinct=TRUE)[233:501])
z.col <- sample(z.col,40)   # pick 40 random colours
show_col(z.col)

z.col42 = c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3",  # each block of 6 = a hue
            "#114578","#185EA5","#1E78D2","#3F91E4","#6CABEA","#98C4F0",
            "#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0","#117845",
            "#18A55E","#1ED278","#3FE491","#6CEAAB","#98F0C4","#787811",
            "#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5",
            "#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498",
            "#781122","#A5182F","#D21E2C","#E43F5B","#EA6C81","#F098A7")
show_col(z.col42)
z.col12 <- c('#A51876','#E43FAD','#1E78D2','#094582','#117878','#3FE4E4',
             '#18A55E','#3FE491','#b8f9d3','#D2D21E','#D2781E','#E43F5B',
             'grey50')
show_col((z.col12))


f2 <-
  ggplot(f2.c, aes(x=Timepoint,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Family~Treat,scales='free',space='free', switch = 'y') +    # on top of the rest, switch labels for y axis
  scale_fill_manual(values=z.col12) +
  theme(strip.text.y = element_text(angle = 0, size =12), axis.text.y = element_text(angle = 0, size=9)) +
  theme(strip.text.x = element_text(size = 14), axis.text.x = element_text(angle = 270, size=11)) +#, strip.text.y = element_text(size = 12))# +
  #stuff from MIHUA's rmarkdown
  theme(panel.background = element_rect(fill=NA, colour=NA)) +                              # rm backgrnd
  theme(panel.grid.major.y = element_line(colour='grey75')) +                               # horiz lines 1
  theme(panel.grid.minor.y = element_line(colour='grey85')) +                               # horiz lines 2
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
#  theme(legend.position = 'none')
  ggtitle("AD of Grass Silage / Slurry: OTUs >2% grouped by Family")
f2
# abundance value overplotting can be taken care of as SVG externally!

# for reference
ggplot(f2.c, aes(x=Timepoint,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Order~Treat,scales='free_x',space='free_x') +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270, size=11)) +
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14)) + 
  ggtitle("the previous version of Rel Ab testing")



## plot all the agglomerated taxa (<0.02% relative abundance), for comparison. 
f2.b <- psmelt(dw_sun.02.23)   # reset f2.b
te_agglom <- ggplot(f2.b, aes(x=Timepoint,y=Abundance,fill=Class)) +
  geom_bar(aes(fill=Class), stat ="identity", position="dodge",  colour="black") +
  facet_grid(Phylum~Treat,scales='free_x',space='free_x' +
               theme(strip.text.y = element_text(angle = 270))) + theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +
  theme(legend.position="none") +
  ggtitle("Relative Abundance Overview of 'Minor' OTUs (<2%)")
te_agglom

## Ladies and Gentlemen; T H E   A R C H A E A!
# everyone's favourite phylogenetic neighbour
# Archaeal community abundance
te_arch <- psmelt(subset_taxa(dw_slv_rb,Domain=='Archaea'))
te_arch_plot <- ggplot(te_arch, aes(x=Timepoint,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Family~Treat,scales='free_x',space='free_x') +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +
  ggtitle("Archaeal OTUs by Family/Orders") +
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14))
te_arch_plot
## ...   *clap*,     *clap*.


##  D I F F E R E N T I A L   A B U N D A N C E    T E S T S   -    L E f S e

##export for testing
## tried doing this intelligently, doing it reliably instead via export to Excel
#write.table(tax_table(dw_slv),'dw_slv_tax.txt',sep='\t')
#write.table(otu_table(dw_slv),'dw_slv_otu.txt',sep='\t')

## import LEfSe output for plotting
dw_lefse <- as.data.frame(read.table('dw_te_edited.txt',sep='\t',header=TRUE,row.names = 1))
#colnames(dw_lefse) <- c('taxa','KW.RST','Setup','LDA','p.val') # already in file, ported with header/row.names above
dw_lefse <- dw_lefse[with (dw_lefse, order(Taxa)),]



## table com, binned with LefSe output in excel to show the relevant stuff: many higher clade are merely
## representations of their OTU, but at a higher confidence due to bootstrapping: these do not contribute
## but do clutter - smooth out!
#dw_lefse_te<-read.table('dw_noted_otus_tabular.txt',sep='\t',header=TRUE,row.names = 1)
dw_lefse_te<-read.table('dw_te_edited2.txt',sep='\t',header=TRUE,row.names = 1)
##


## kept getting bad orders of numeric data: FIX: sort by column A, then by coulmn B:
dw_lefse_te <- dw_lefse_te[with(dw_lefse_te, order(Condition, -KW.RST)), ]
## still need to stop ggplot reordering: so 'fix' the order by factorising:
dw_lefse_te$Taxa <- factor(dw_lefse_te$Taxa, levels = dw_lefse_te$Taxa)



ggplot(dw_lefse_te, aes(x=Taxa, y=KW.RST, fill= Condition, shape=Condition)) +
  theme(axis.text = element_text(size=13), legend.text=element_text(size=13), axis.title = element_text(size=13)) +
  theme(strip.text.x = element_text(size = 15)) + 
  facet_grid(~Condition) +
  geom_point(aes(size=LDA)) +
  coord_flip() +
  guides(fill=guide_legend(title.theme=element_text(size=21, face= 'bold', colour='green'))) +
  labs(title='Taxa Associations with Reactor Setup') +
  labs(y='Magnitude of Difference (H score)') +
  scale_size_continuous(name='LDA Effect Size (log10)') +
  theme(panel.background = element_rect(fill=NA, colour=NA)) +                              # rm backgrnd
  theme(panel.grid.major.x = element_line(colour='grey80')) +                               # horiz lines 1
  theme(panel.grid.major.y = element_line(colour='grey80')) +                               # horiz lines 1
  scale_shape_manual(values = c(22, 21, 24)) +
  scale_fill_manual(values = c(brewer.pal(9,'YlOrRd')[8],brewer.pal(3,'Greens')[3],brewer.pal(3,'YlGnBu')[3])) #+
#  xlim(35)   #fuuuuck do it in inkscape
#  scale_x_continuous(breaks = c(3,4,5))
#  theme(panel.grid.minor.y = element_line(colour='grey85')) #+                           # horiz lines 2
#  theme(panel.grid.major.x = element_blank())                                             # rm vert lines
  



## V E G A N    C C A   P L O T
# vegan friendly map: DF
cca_map <- as.data.frame(as.matrix(dw_slv_map))
# fix Beff values as numeric values for CCA 
cca_map$Beff<-as.numeric(as.character(cca_map$Beff))
# fix order of rows for map, as sample #04 is timepoint #12 due to barcode use (not important, just a little awkward)
cca_map<-cca_map[order(rownames(cca_map)),]
# Weeks as samplenames for clarity
rownames(cca_map)<-cca_map$Timepoint


##from OTU_table to sorted, coloured plot:
veg.3up=as.data.frame(t(otu_table(dw_slv))) # 14 X 344 df of the OTU abundances from phyloseq
rownames(veg.3up)<-cca_map$Timepoint # match the rownames in cca_map
#veg.5.Fam=c(tax_table(dw_slv)[,5]) # the Fam slot data for dw_slv

# names of OTUs significant in LEfSe are changed to a dummy value for LDA score bins (a-h), so can be plotted by vegan.
# colouring by veg.rassoclev below plots that category as a colour in conjunction with veg.9cols 
# non-significant OTUs simply relabelled 'unassoc. taxa'
veg.rassoc<-read.table('dw_noted_otus.txt',sep='\t',header=FALSE,row.names = 1)
# transpose so maches otu table (OTUs as col), to DF, and use to rename OTU table columns::
veg.rassoc<-as.data.frame(t(veg.rassoc))
veg.3up2 <- as.matrix(veg.3up)
colnames(veg.3up2) <-as.matrix(veg.rassoc[1,1:344])  # couldnt make this work any other way
veg.3up <- as.data.frame(veg.3up2)
# vector to carry association levels as factors:
veg.rassoclev <-c('R-G:         LDA ≦3.25',
                  'R-G-TE:   LDA ≦3.25',
                  'R-G-TE:   LDA 3.26-3.5',
                  'R-G-TE:   LDA 3.6-4.0',
                  'R-G-TE:   LDA >4.0',
                  'R-SG:      LDA ≦3.25',
                  'R-SG:      LDA 3.26-3.5',
                  'R-SG:      LDA 3.6-4.0',
                  'unassoc. taxa')

# setup CCA object
veg.cca3c<-vegan::cca(veg.3up~Beff,cca_map)
# set plotting scale
veg.scl<-c(-2)

# the colour and the shape
# make three sets; green red blue for setups
# 14 colours for 14 timepoints
veg.14cols<-c(brewer.pal(9,'YlOrRd')[2:9],brewer.pal(7,'Greens')[6:7],brewer.pal(9,'YlGnBu')[5:9])
veg.chars <- c(22, 21, 24)
veg.2col <- c('rosybrown1','lightblue')   # for ordiellipse bg

# OTUs
v.blues <- brewer.pal(9,'Blues')
v.greens <- brewer.pal(9,'Greens')
veg.9cols <- c('maroon')
veg.9cols<-c(veg.9cols,v.greens[7:9],v.greens[4])
veg.9cols <- c(veg.9cols,v.blues[7:8],v.blues[5])
veg.9cols <- c(veg.9cols,'gray85')
#show_col(veg.9cols)

# plot that stuff out!
veg.phy <- plot(veg.cca3c, type='n', display='species', scaling=veg.scl, main='  ') # plot empty
with (cca_map, ordiellipse(veg.cca3c, Reactor, col=veg.2col, draw='polygon', kind='sd', label=TRUE, cex=1.2, alpha=0))  # SD ellipse
with (veg.3up, points(veg.cca3c, display='spec', col = 'gray30', bg=veg.9cols[as.factor(colnames(veg.3up))], pch=21, cex=0.9)) # colour by fav families
with (cca_map, points(veg.cca3c, display='sites', col='black', bg=veg.14cols[cca_map$Timepoint], pch=veg.chars[as.factor(cca_map$Treat)],cex=2)) # reactor points and colours
with(cca_map,text(veg.cca3c,display='bp',col='black',cex=1.4))
with(veg.phy, legend("topleft", legend = levels(as.factor(veg.rassoclev)), bty = "n", col = 'gray30', pch = 21, pt.bg = veg.9cols, cex=1.2, pt.cex=1))
with(cca_map, legend("topright", legend = levels(as.factor(dw_slv_map$Timepoint)), bty = "n", col = 'gray30', pch = veg.chars[as.factor(dw_slv_map$Treat)], pt.bg = veg.14cols, cex=1.2, pt.cex=1.4))

# end workflow
dev.off()







#
#
#
## ## F I G U R I N G   I T   O U T 
# #dw refiguring 
# 
# 
# ## doublecheck association
# z <- c('silva4_107','silva4_108','silva4_109','silva4_110','silva4_111','silva4_112','silva4_113','silva4_114','silva4_115','silva4_116','silva4_117','silva4_118','silva4_119','silva4_120','silva4_121','silva4_122','silva4_123','silva4_124','silva4_125','silva4_126','silva4_127','silva4_128','silva4_129','silva4_130','silva4_131','silva4_132','silva4_133','silva4_134','silva4_135','silva4_136','silva4_137','silva4_138','silva4_139','silva4_140','silva4_141','silva4_142','silva4_143','silva4_144','silva4_145','silva4_146','silva4_147','silva4_148','silva4_149','silva4_150','silva4_151','silva4_154','silva4_156','silva4_157','silva4_158','silva4_159','silva4_160','silva4_161','silva4_162','silva4_163','silva4_164','silva4_165','silva4_166','silva4_167','silva4_169','silva4_170','silva4_171','silva4_172','silva4_173','silva4_174','silva4_175','silva4_177','silva4_178','silva4_179','silva4_180','silva4_181','silva4_182','silva4_183','silva4_184','silva4_185','silva4_186','silva4_187','silva4_188','silva4_189','silva4_190','silva4_191','silva4_193','silva4_194','silva4_195','silva4_196','silva4_197','silva4_198','silva4_199','silva4_200','silva4_202','silva4_203','silva4_204','silva4_205','silva4_206','silva4_207','silva4_208','silva4_209','silva4_210','silva4_211','silva4_212','silva4_213','silva4_214','silva4_215','silva4_216','silva4_217','silva4_218','silva4_219','silva4_220','silva4_221','silva4_222','silva4_223','silva4_224','silva4_225','silva4_226','silva4_227','silva4_229','silva4_230','silva4_231','silva4_232','silva4_233','silva4_234','silva4_235','silva4_236','silva4_237','silva4_238','silva4_239','silva4_240','silva4_241','silva4_243','silva4_244','silva4_245','silva4_247','silva4_248','silva4_249','silva4_251','silva4_252')
# View(as.data.frame(tax_table(dw_slv)[z,]))
# # lots of clostridia! explains emphasis in first graph?
#
# # get clostridial totals for R-G
# ## "associated strongly with increased Class Clostridia abundance (LDA: 4.4), 
# ##  with reads increasing by Y% over average of X%% to 90% (Figure 2)." 
# 
# z <- subset_samples(dw_slv, Reactor=='G')
# length(get_taxa_unique(subset_taxa(z, Class=='Clostridia'),taxonomic.rank = 'Genus'))
# zz <- subset_samples(dw_slv, Reactor=='SG')
# length(get_taxa_unique(subset_taxa(zz, Class=='Clostridia'),taxonomic.rank = 'Genus'))
# otu_table(z)[(taxa_sums(z) == 0),]
# sample_sums(subset_taxa(z, Class=='Clostridia')) / sample_sums(z)
# 
# zz <- subset_samples(z, Treat=='silage')
# mean(sample_sums(subset_taxa(dw_slv, Class=='Clostridia')) / sample_sums(dw_slv))
# zzz <- subset_samples(z, Treat=='silage w.TE')
# mean(sample_sums(subset_taxa(zzz, Class=='Clostridia')) / sample_sums(zzz))
# 
# # Clostridia NOT more diverse in R-G despite claims
# plot_richness(dw_slv,x='Reactor',color='WeekNo') #subset_samples(dw_slv, WeekNo != '1')
# sample_data(dw_slv)
# 
# z <- rownames(otu_table(dw_slv)[(taxa_sums(dw_slv) >= 3),])
# z <- prune_taxa(z, dw_slv)
# plot_richness(z, x='Reactor', color='WeekNo')
# 
# #z <- subset_samples(dw_slv, Reactor=='G')
# z <- (subset_taxa(dw_slv, Domain=='Archaea')) #/ sample_sums(z)
# (taxa_sums(z)) / sum(taxa_sums(dw_slv))
# zz <- subset_samples(z, Treat=='silage')
# mean(sample_sums(subset_taxa(zz, Class=='Clostridia')) / sample_sums(zz))
# # [1] 0.8714049
# zzz <- subset_samples(z, Treat=='silage w.TE')
# mean(sample_sums(subset_taxa(zzz, Class=='Clostridia')) / sample_sums(zzz))
# # [1] 0.8964196
# 





##                                                       ##
#                                                         #
#                                                         #
##      Specifically recreate the paper figures:         ##
#                                                         #
#                                                         #
##                                                       ##



# FIGURE 1A, 1B
# Created in excel, not reproduced here. Data supplied as tabular table in supporting info.

# FIGURE 2
f2 <-
  ggplot(f2.c, aes(x=Timepoint,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Family~Treat,scales='free',space='free', switch = 'y') +    # on top of the rest, switch labels for y axis
  theme(strip.text.y = element_text(angle = 0, size =12), axis.text.y = element_text(angle = 0, size=11)) +
  theme(strip.text.x = element_text(size = 14), axis.text.x = element_text(angle = 0, size=11)) +#, strip.text.y = element_text(size = 12))# +
  #stuff from MIHUA's rmarkdown
  theme(panel.background = element_rect(fill=NA, colour=NA)) +                              # rm backgrnd
  theme(panel.grid.major.y = element_line(colour='grey75')) +                               # horiz lines 1
  theme(panel.grid.minor.y = element_line(colour='grey85')) +                               # horiz lines 2
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
  #  theme(legend.position = 'none')
  ggtitle("AD of Grass Silage / Slurry: OTUs with Rel.Abundances above 2% by Genus/Orders")
f2


# FIGURE 3a
ggplot(dw_lefse_te, aes(x=Taxa, y=KW.RST, fill= Condition, shape=Condition)) +
  theme(axis.text = element_text(size=13), legend.text=element_text(size=13), axis.title = element_text(size=13)) +
  labs(title='Taxa Associations with Reactor Setup') +
  scale_size_continuous(name='LDA Effect Size (log10)') +
  geom_point(aes(size=LDA)) +
  coord_flip() +
  guides(fill=guide_legend(title.theme=element_text(size=21, face= 'bold', colour='green'))) +
  theme(strip.text.x = element_text(size = 15)) + facet_grid(~Condition) +
  labs(y='Magnitude of Difference (H score)') +
  theme(panel.background = element_rect(fill=NA, colour=NA)) +                              # rm backgrnd
  theme(panel.grid.major.x = element_line(colour='grey80')) +                               # horiz lines 1
  theme(panel.grid.major.y = element_line(colour='grey80')) +                               # horiz lines 1
  scale_shape_manual(values = c(22, 21, 24)) +
  scale_fill_manual(values = c(brewer.pal(9,'YlOrRd')[8],brewer.pal(3,'Greens')[3],brewer.pal(3,'YlGnBu')[3])) #+

# FIGURE 3B
# ordihull alpha set to zero!
veg.phy <- plot(veg.cca3c, type='n', display='species', scaling=veg.scl, main='  ') # plot empty
with (cca_map, ordiellipse(veg.cca3c, Reactor, col=veg.2col, draw='polygon', kind='sd', label=TRUE, cex=1.2, alpha=0))  # SD ellipse
with (veg.3up, points(veg.cca3c, display='spec', col = 'gray30', bg=veg.9cols[as.factor(colnames(veg.3up))], pch=21, cex=0.9)) # colour by fav families
with (cca_map, points(veg.cca3c, display='sites', col='black', bg=veg.14cols[cca_map$Timepoint], pch=veg.chars[as.factor(cca_map$Treat)],cex=2)) # reactor points and colours
with(cca_map,text(veg.cca3c,display='bp',col='black',cex=1.4))
with(veg.phy, legend("topleft", legend = levels(as.factor(veg.rassoclev)), bty = "n", col = 'gray30', pch = 21, pt.bg = veg.9cols, cex=1.2, pt.cex=1))
with(cca_map, legend("topright", legend = levels(as.factor(dw_slv_map$Timepoint)), bty = "n", col = 'gray30', pch = veg.chars[as.factor(dw_slv_map$Treat)], pt.bg = veg.14cols, cex=1.2, pt.cex=1.4))




##                                                       ##
#                                                         #
#                                                         #
##      Specifically recreate Supporting figures:        ##
#                                                         #
#                                                         #
##                                                       ##


## SF1: Rarefaction curve: produced by silva, not in R.

## SF2: CA with respect to Sequencing Depth
plot_ordination(dw_slv,ordinate(dw_slv,'CCA'),color='SeqDepth',shape='Treat',title='CA of dw showing sequencing depth',label='Timepoint',type='biplot')

## SF3: Archaeal Community Abundance
te_arch <- psmelt(subset_taxa(dw_slv_rb,Domain=='Archaea'))
te_arch_plot <- ggplot(te_arch, aes(x=Timepoint,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Family~Treat,scales='free_x',space='free_x') +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +
  ggtitle("Archaeal OTUs by Family/Orders") +
  theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14))

ggplot(te_arch, aes(x=Timepoint,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Family~Treat,scales='free',space='free', switch = 'y') +    # on top of the rest, switch labels for y axis
  theme(strip.text.y = element_text(angle = 0, size =19), axis.text.y = element_text(angle = 0, size=19)) +
  theme(strip.text.x = element_text(size = 15), axis.text.x = element_text(angle = 270, size=13)) +#, strip.text.y = element_text(size = 12))# +
  theme(panel.background = element_rect(fill=NA, colour=NA)) +                              # rm backgrnd
  theme(panel.grid.major.y = element_line(colour='grey75')) +                               # horiz lines 1
  theme(panel.grid.minor.y = element_line(colour='grey85')) +                               # horiz lines 2
  theme(panel.grid.major.x = element_blank()) +                                             # rm vert lines
  #  theme(legend.position = 'none')
  theme(legend.text=element_text(size=19)) +
  scale_fill_manual(values =brewer.pal(8, name = 'Set1')) #, limits = c("4", "6", "8", "10"))
#scale_fill_manual(values = sample(z.col42, 8)) #, limits = c("4", "6", "8", "10"))

  ggtitle("AD of Grass Silage / Slurry: OTUs with Rel.Abundances above 2% by Genus/Orders")


te_arch_plot

## un-SF, unmentioned: Diversity indices, coloured by Sequencing depth
#plot_richness(dw_slv, x = "Treat", color = "SeqDepth", shape = 'Reactor') #+ geom_boxplot()

