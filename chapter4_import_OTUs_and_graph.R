## displaced suplhur reducer - one-C?
## heat selection> State D?
## loss of A-Proteo/mysterion: inoculum deg?
## methanobacteria happiest wiht otherwise shitty conditions - where is the blockage?State C higly uniform but for Mollicutes/M'bacter split
#

setwd("/media/cart/Dropbox/SilentGeno/R/R_AJ/R_AJ_SILVA")
setwd("C:/Users/Blip the Tall/Documents/Jamie blipping/Dropbox/SilentGeno/R/R_AJ/R_AJ_SILVA")
#install.packages('tidyr')
library('tidyr')
library('phyloseq')
library('ggplot2')
library('DESeq2')
library('ancom.R')
library('vegan')

# this map differs from other map by way of column names
aj_map <- 'aj_map.txt'
aj_map <- import_qiime_sample_data(aj_map)
# aj_map$State <- NULL
# aj_map$State <- c('B','B','B','C','C','C','D','D','D')
# aj_map$Beff <- aj_map$SMY/c(400,400,400,400,400,400,400,400,400,400,348,348,348,348)
# aj_map$facOLR <- c('two','two','three','three','three','three-five','three-five','four','four','four','two','two','three','four')
# aj_map$WeekNo<- as.numeric(as.character(c('01','20','37','40','43','49','54','64','65','68','01','20','37','65')))


# pre-moulded to fit our needs
aj_input<-read.table('aj_silva_edOTU.tsv',sep='\t',header=TRUE,row.names = 1)
aj_ab <- aj_input[,1:9]
dim(aj_ab)
# cut out NO RELATIVES
aj_ab <- aj_ab[1:114,]
aj_ab <- otu_table(aj_ab,taxa_are_rows = TRUE)
#aj_tax<-(separate(data = aj_input, col = 'taxonomy', into = c(Domain","Phylum","Class","Order","Family","Genus","OTU"), sep = ";"))[,16:21]
aj_tax<-(aj_input[1:114,10:16])
aj_tax <- as.matrix(aj_tax,rownames.force = TRUE)
aj_tax<-tax_table(aj_tax)
row.names(aj_tax) <- row.names(aj_input)
# aj_tax[,7]<-c(rownames(aj_tax))
colnames(aj_tax) <- c('Domain','Phylum','Class','Order','Family','Genus','OTU')
aj<-merge_phyloseq(aj_ab,aj_tax,aj_map)
sample_data(aj)[,'SeqDepth']<-sample_sums(aj)
aj
# fix MBA03 - done in exceloid
#tax_table(aj)['silva4_253',] <-c('Bacteria','Firmicutes','Clostridia','MBA03','MBA03 F.','MBA03 G.','silva4_253')
# curtail vadinBC27
#tax_table(aj)['silva4_051',] <-c('Bacteria','Bacteroidetes','Bacteroidia','Bacteroidales','Rikenellaceae','vadinBC27 group','silva4_051')


## T R A N S F O R M   -   R E L . B U N
# rel.bun
aj_rb <- transform_sample_counts(aj, function(x)x/sum(x))
# hellinger
aj_h<-transform_sample_counts(aj,function(x) sqrt(x/sum(x)))

#jfgformation
## identical in barchart (as just a common sum) but different in DA testing?
aj_jfg<-transform_sample_counts(aj, function(x) (x/sum(x))*max(sample_sums(aj)))

## S U B S E T
prun.aj_s.3.23 = genefilter_sample(aj, filterfun_sample(function(x) x >=3),A=0.23*nsamples(aj))
aj_s.3.23 = prune_taxa(prun.aj_s.3.23, aj)
aj_sun.3.23 = prune_taxa(!prun.aj_s.3.23, aj)
aj_s.3.23_rb <- transform_sample_counts(aj_s.3.23, function(x)x/sum(x))
# arbitrary rel.abundance cutoff: >1%, x0.23 samples
prun.aj_s.01.23 = genefilter_sample(aj_rb, filterfun_sample(function(x) x >=0.01),A=0.23*nsamples(aj_rb))
aj_s.01.23 = prune_taxa(prun.aj_s.01.23, aj_rb)
# N O T E: negation of logical vector with prepending '!' !!
aj_sun.01.23 = prune_taxa(!prun.aj_s.01.23, aj_rb)
# # flatten all family for CCA plot
# tax_table(aj_sun.02.23)[,'Family'] <- c('F:Misc & NA')
# aj_s.sun <- merge_phyloseq(aj_s.02.23,aj_sun.02.23)




# just aj .02 x0.23

## record the taxa and OTUs

write.table( c( otu_table(aj), as.data.frame(tax_table(aj))), 'aj_abundance_taxonomy.txt', sep='\t')
# write.table(tax_table(aj),'aj_tax.txt',sep='\t')
# write.table(otu_table(aj),'aj_otu.txt',sep='\t')

## P L O T T I N G

pdf('aj_silva-diversions.pdf',width=10, height=7.5)

## P R O C E S S - note change in order of samples to Replicate
# z.rm <- sample_data(aj)[with(sample_data(aj),order(Replicate)),]
# z.rm[,17:19]<-NULL
# z.rm[,1:5]<-NULL
# z.rm[,'State']<-NULL
# z.rm[,'w_teor']<-NULL
# z.rm[,'NH3']<-NULL
# z.rm[,'CH4_pc']<-NULL
# z.rm[,'Replicate']<-NULL
# z.rm[,'Week']<-NULL

# plot.ts(z.rm[1:10,1:8],main='Process for G samples')
# plot.ts(z.rm[10:14,1:8],main='Process for SG samples')
# #plot.ts(z.rm[,1:8],main='Process for G-SG samples')


## S U M M A R Y   S T A T S  ##
# diversity metrics: 
plot_richness(aj, 'Replicate', x = "State", color = "SeqDepth", shape = 'Replicate') #+ geom_boxplot()
# seq depth
plot_ordination(aj, ordinate(aj,'CCA'), color='SeqDepth', shape='State', title='CA-plot of AJ Samples w./ Sequencing Depth', label='Replicate', type='samples') +
  geom_point(size=5) #+
  # geom_label()
  # geom_text(aes(colu)) +
  # theme(strip.background = element_blank()) #+


plot_bar(aj, 'Replicate', fill='Class',title='aj-SILVA Absolute Read Abundance') + facet_grid(~State,scales='free_x',space='free_x') #
plot_bar(aj_rb, 'Replicate', fill='Family',title='aj-SILVA Relative Read Abundance') + facet_grid(Domain~State,scales='free',space='free') #
# plot_bar(aj_s.01.23 ,'Replicate',fill='Family',title='aj-SILVA.3.23 Relative Read Abundance') +
#   facet_grid(State~Class,scales='free',space='free_y') +
#   theme(axis.text.x = element_text (angle=-90, hjust=0, vjust=0.5, size=10 )) +
#   theme(strip.text.y = element_text(angle=0, size=10))


z.rm3 <- psmelt(aj_s.01.23)
z.rm4 <- psmelt(aj_sun.01.23)
#plot the agglomerated OTUs for reference of phylum abundance
## pool all non->2% to dummy phylum, append to bottom of melted frame 
z.rm4[,14] <- c('P:Other')
z.rm4[,15] <- c('C:Other')
z.rm4[,16] <- c('O:Other')
z.rm4[,17] <- c('F:Misc & NA')
z.rm4[,18] <- c('G:Misc')
z.rm5 <- rbind(z.rm3, z.rm4)

ggplot(z.rm5, aes(x=Replicate,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(State~Class,scales='free_x',space='free_x', switch='y') +
  theme(strip.text.y = element_text(angle = 180, size=14), strip.text.x = element_text(angle = 90, size=11), axis.text.x = element_text(angle = 270, size=10)) +
  ggtitle("AJ OTUs with Rel.Abundances above 1% by Class/Genus") +
  theme(strip.background = element_blank()) #+

## sort colouring of e.g. genus by class
# z.rm5 <- z.rm5[with(z.rm5, order(Class, Genus)), ]
## fix new order by  making a factor
# z.rm5$Genus <- factor(z.rm5$Genus, levels = z.rm5$Genus)

## Ladies and Gentlemen; THE ARCHAEA!
z.rmA <- psmelt(subset_taxa(aj_rb,Domain=='Archaea'))
ggplot(z.rmA, aes(x=Replicate,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Genus~State,scales='free',space='free_x') +
  theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 270)) +
  ggtitle("AJ Archaeal OTUs by Genus/Order (note y-axis scales)")
## ...   *clap*,     *clap*.


# everything else
z.rm4 <- psmelt(subset_taxa(aj_sun.01.23,Phylum!='Firmicutes'))   # reset z.rm4
ggplot(z.rm4, aes(x=Replicate,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Order~State,scales='free_x',space='free') + #, switch='y'
  theme(strip.text.y = element_text(angle = 0, size=10), strip.text.x = element_text(size=14), axis.text.x = element_text(angle = 270, size=14)) +
  theme(legend.position="none") +
  ggtitle("AJ Overview of 'Minor' OTUs Rel.Abundances (<1%, except Firmicutes)")

# try subsetting to minor firmicutes only
z.rm4 <- psmelt(subset_taxa(aj_sun.01.23,Phylum=='Firmicutes'))   # reset z.rm4
ggplot(z.rm4, aes(x=Replicate,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Order~State,scales='free_x',space='free') + #, switch='y'
  theme(strip.text.y = element_text(angle = 0, size=10), strip.text.x = element_text(size=14), axis.text.x = element_text(angle = 270, size=14)) +
  theme(legend.position="none") +
  ggtitle("AJ Overview of 'Minor FIRMICUTE' OTUs Rel.Abundances (<1%)")



## RE THINK ORDINATION: aj_SLV PROB NOT APPROP RE: ABUNDANCES
## ordinate same scale as tested in LEfSe, as introduces least amount of bullshit:
# aj_lef <- aj
# otu_table(aj_lef)<-otu_table(apply(otu_table(aj_lef),2,function(x)x/1000000),taxa_are_rows = TRUE)
## O R D I N A T E
# NMDS: stress is nearly zero- insufficient data?
plot_ordination(aj,ordinate(aj,'DCA'),color='Phylum',shape='State',title='CCA of AJ',label='Replicate',type='biplot')
plot_ordination(aj_h,ordinate(aj,'CCA'),color='Phylum',shape='State',title='CA of AJ',label='Replicate',type='biplot')
# can use transfomred values, but affects the result very strongly. not advisable as translated to X-sq distances?

dev.off()

z.rm <- psmelt(subset_taxa(aj_s.01.23,Phylum=='Firmicutes'))
ggplot(z.rm, aes(x=Replicate,y=Abundance,fill=Genus)) +
  geom_bar(aes(fill=Genus), stat ="identity", position="stack",  colour="black") +
  facet_grid(Family~State,scales='freex',space='freex') + #, switch='y'
  theme(strip.text.y = element_text(angle = 0, size=10), strip.text.x = element_text(angle = 0, size=12), axis.text.x = element_text(angle = 270, size=12)) +
  ggtitle("AJ OTUs with Rel.Abundances above 1% by Class/Genus") +
  theme(strip.background = element_blank()) #+

