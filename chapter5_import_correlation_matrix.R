
# mv_cor.R

library('heatmaply')
library('ape')
library('phyloseq')
library('ggplot2')
library('DESeq2')
library('scales')

## IMPORT COR & PP MATRICES
mv_cor <- read.table('C:/Users/Blip the Tall/Downloads/mv_sparcc/mv_cor.txt',sep='\t',header=TRUE,row.names=1)
mv_pp <- read.table('C:/Users/Blip the Tall/Downloads/mv_sparcc/mv_two_sided_pp.txt',sep='\t',header=TRUE,row.names=1)

# sort phylo object and then pass to decendants 
# mv_sort <- merge_phyloseq(otu_table(mv),tax_table(mv),sample_data(mv)) ## if necessary to remove tree
mv_sort <- mv
tax_table(mv_sort) <- tax_table(mv_sort)[with(as.data.frame(tax_table(mv_sort)), order(Kingdom,Phylum,Class,Order,Family,Genus)),]
otu_table(mv_sort) <- otu_table(mv_sort)[,with(as.data.frame(tax_table(mv_sort)), order(Kingdom,Phylum,Class,Order,Family,Genus))]
tax_table(mv_sort)[45:55,1:5]
# DOES NOT WORK if includes a tree

# propagate subset of 2300x2300
mv_sort_id <- rownames(tax_table(mv_sort)) ; head(mv_sort_id)
mv_cor2 <- mv_cor[mv_sort_id,mv_sort_id] ;   mv_cor2[1:10,1:10]
mv_pp2 <- mv_pp[mv_sort_id,mv_sort_id] ;     mv_pp2[1:10,1:10]

## NUKE CORR via PP >0.1  - doing early avoids confounding later due to autocorelation=1
z.pNA <- apply(mv_pp2,c(1,2),function(x)x>0.01)    # logical matrix of values >0.01
# consider a more stringent value to reduce size of set
mv_cor2[z.pNA] <- c(0)                              # index cor_10 via log.mat and apply NA


# ## REMOVE NON-INFORMATIVE COLUMNS/ROWS   -   BEFORE you de-phylo the names!
# # need to do this before calculating the dendrogram for correct number of tips
#   z.1 <- Filter(function(x) all(x < 0.80 & x > -0.80), mv_cor2) ; dim(z.1)         # identify the shit ones
#   z.2 <- !(rownames(mv_cor2) %in% colnames(z.1))
#   mv_.80 <- subset_taxa(mv_sort, z.2 )                                             # remove shit ones from phylo
#   cor_.80 <- mv_cor2[ z.2 , z.2 ]    # remvoe shit ones from COR-MATR
#   dim(cor_.80)# ; heatmaply(cor_.80) 
#   
#   ## NOW mutilate/propagate names from subset phylobject
#   rownames(cor_.80) <- paste( sub('mvsv', '', rownames(tax_table(mv_.80)) ) , tax_table(mv_.80)[,6] ) #, 1:nrow(cor_10) )   # need unique, meaningful names - mvsv_ID and Genus
#   colnames(cor_.80) <- paste( sub('mvsv', '', rownames(tax_table(mv_.80)) ) , tax_table(mv_.80)[,6] ) 
#   # as values nuked above, no longer necessary to propagate mv_pp names etc.
#   cor_.80[1:10,1:10]  
#   
#   ##
#   heatmaply(cor_.80, Rowv=FALSE, col=c("#440154FF","#21908CFF","#27AD81FF","#21908CFF","#FDE725FF"))
#   # shows clusters, but arguably not as well as would symmetric heatmap (dual hclust'ing)
# #heatmaply(cor_.80)


#=====
# further subsets - DA from DESeq2 and SV from mv_02
#mv_.80 <- subset_taxa(mv_sort, z.2 )                      # subset DESeq DAs
z.cor1 <- mv_cor2[rownames(tax_table(mv_01)), rownames(mv_in_r1_sig)]
rownames(z.cor1) <- paste( sub('mvsv', '', rownames(tax_table(mv_01)) ) , tax_table(mv_01)[,6] ) 
colnames(z.cor1) <- paste( sub('mvsv', '', rownames(mv_in_r1_sig))  , tax_table(mv)[rownames(mv_in_r1_sig),6] ) #, 1:nrow(cor_10) ) 
# NOW remove non-informative rows, helping x-axis legibility
# or... necessary? 
# dont have great confidence in your filtering code from before.
    # z.1 <- Filter(function(x) all(x < 0.70 & x > -0.70), mv_cor2) ; dim(z.1)    # prob incorrect     
    # z.2 <- z.cor1[ , -c(colnames(z.cor1) %in% colnames(z.1)) ] ; dim(z.2)  #gives 0 cols

# heatmaply_cor(z.cor1, col=c("#440154FF","#21908CFF","#21908CFF","#27AD81FF","#21908CFF","#21908CFF","#FDE725FF") ) #, file='mv_co1a_stretch.html')
# dim(z.cor1)

# cols 38:50 are dead weight
z.1 <- Filter(function(x) all(x < 0.75 & x > -0.75), z.cor1) ; dim(z.1)         # identify the shit ones
z.2 <- !(colnames(z.cor1) %in% colnames(z.1))
z.co1a <- z.cor1[  , z.2]
heatmaply(z.co1a, col=c("#440154FF","#21908CFF","#21908CFF","#27AD81FF","#21908CFF","#21908CFF","#FDE725FF")) #, file='mv_co1a_stretch.html')

#==
# simple OTU tbale correlation  
z.sv <- c( rownames(mv_in_r1_sig)) ; z.sv <- unique(z.sv) #rownames(tax_table(mv_02)) , 
z.sv <- as.data.frame(as.matrix( otu_table(mv_h)[,z.sv] ))
#rownames(z.sv) <- paste( sub('mvsv', '', rownames(mv_in_r1_sig) ) , tax_table(mv)[rownames(mv_in_r1_sig),6] )
colnames(z.sv) <- paste( sub('mvsv', '', rownames(mv_in_r1_sig) ) , tax_table(mv)[rownames(mv_in_r1_sig),6] )

heatmaply(z.sv)  
  
# presumably ES:
# 31
# 40
# 19
# 20
# ES: 31 40 19 192 146 70 82 19 161 177 20 24 75 113 159 103

# 
# presumably IS:
# 14
# 52
# IS: 14 244 76 52 192 145 66 114 33 214 127 110 45 41 59 276 

#
# additional cohorts:
# #Gelria/Treponema:
# z.gt <- c('mvsv0183','mvsv0212','mvsv0106','mvsv0067','mvsv0152','mvsv0255')
# #W5/NB1-n:
# z.WN <- c('mvsv0031','mvsv0099','mvsv0158','mvsv0204','mvsv0040','mvsv0129','mvsv0150')




####   = = =    = = =   = = =    = = =   = = =    = = =   = = =    = = =
####   = = =    = = =   = = =    = = =   = = =    = = =   = = =    = = =
####   = = =    = = =   = = =   A B A N D O N W A R E =   = = =    = = =
####   = = =    = = =   = = =    = = =   = = =    = = =   = = =    = = =
####   = = =    = = =   = = =    = = =   = = =    = = =   = = =    = = =



# z.cor2 <- mv_cor[ rownames(mv_in_r1_sig) , rownames(mv_in_r1_sig) ]
# rownames(z.cor2) <- paste( sub('mvsv', '', rownames(mv_in_r1_sig) ) , tax_table(mv)[rownames(mv_in_r1_sig),6] ) 
# colnames(z.cor2) <- paste( sub('mvsv', '', rownames(mv_in_r1_sig) ) , tax_table(mv)[rownames(mv_in_r1_sig),6] ) 
# heatmaply(z.cor2, col=c("#440154FF","#21908CFF","#27AD81FF","#21908CFF","#FDE725FF"))
# 
# 
# z.cor3 <- mv_cor[ rownames(tax_table(mv_02)), rownames(tax_table(mv_02)) ]
# rownames(z.cor3) <- paste( sub('mvsv', '', rownames(tax_table(mv_02)) ) , tax_table(mv_02)[,6] ) 
# colnames(z.cor3) <- paste( sub('mvsv', '', rownames(tax_table(mv_02)) ) , tax_table(mv_02)[,6] ) 
# heatmaply(z.cor3, col=c("#440154FF","#21908CFF","#27AD81FF","#21908CFF","#FDE725FF"))
# 
# summary(z.cor4)
# z.cor4 <- mv_cor[ rownames(tax_table(mv)[1:100,]), rownames(tax_table(mv)[1:100,]) ]
# rownames(z.cor4) <- paste( sub('mvsv', '', rownames(tax_table(mv)[1:100,]) ) , tax_table(mv)[1:100,6] ) 
# colnames(z.cor4) <- paste( sub('mvsv', '', rownames(tax_table(mv)[1:100,]) ) , tax_table(mv)[1:100,6] ) 
# heatmaply(z.cor4, Rowv=FALSE, Colv=FALSE, col=c("#440154FF","#21908CFF","#27AD81FF","#21908CFF","#FDE725FF"))
# 
# 
# mv_note <- c( rownames(mv_in_r1_sig), colnames(otu_table(mv_02)) )
# mv_note <- unique(mv_note)
# mv_noted <- mv_cor2[mv_note, mv_note]
# rownames(mv_noted) <- paste( sub('mvsv', '', mv_note) , tax_table(mv)[mv_note,6] ) 
# colnames(mv_noted) <- paste( sub('mvsv', '', mv_note) , tax_table(mv)[mv_note,6] ) 
# heatmaply(mv_noted)
