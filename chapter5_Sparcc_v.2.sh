#! /usr/bin/zsh

## next iteration of this, with rlog-transformed MV SV table
## called using screen to detach command
# ssh -t jfg@143.239.200.154 screen -dm -S mv_sparcc /media/cart/Dropbox/SilentGeno/scripts/Sparcc_mv_2.sh

cd /home/jfg/BioApps/sparcc/mv_sparcc2


## 1: Test OTU Table for Covariance, with default iterations and xiterations (10,20)
time python2.7  /home/jfg/BioApps/sparcc/SparCC.py  /home/jfg/BioApps/sparcc/mv_sparcc2/mv_sv_rlog.txt --cor_file= /home/jfg/BioApps/sparcc/mv_sparcc2/mv_cor.txt --cov_file= /home/jfg/BioApps/sparcc/mv_sparcc2/mv_cov.txt


## 2: Iterate a load (500!) of shuffled OTU tables: (fast!)
time python2.7  /home/jfg/BioApps/sparcc/MakeBootstraps.py  /home/jfg/BioApps/sparcc/mv_sparcc2/mv_sv_rlog.txt -n 500 -t mv_permd_#.txt -p  /home/jfg/BioApps/sparcc/mv_sparcc2/mv_boots/


## 3.a: change directory for the next step to work: 
cd mv_sparcc2/mv_boots

## 3.b: SparCC on each bootstrap table made above (super. fucking. slow.):
time for i in $(ls *.txt); do python2.7  /home/jfg/BioApps/sparcc/SparCC.py $i --cor_file cor_"$i" --cov_file cov_"$i"; done


## 4: Generate Pseudo-P Values for Correlation using Iterated, Correlated Datasets:
time python2.7  /home/jfg/BioApps/sparcc/PseudoPvals.py '/home/jfg/BioApps/sparcc/mv_sparcc2/mv_cor.txt' '/home/jfg/BioApps/sparcc/mv_sparcc2/mv_boots/cor_mv_permd_#.txt' 500 -o /home/jfg/BioApps/sparcc/mv_sparcc2/mv_two_sided_pp.txt -t two_sided
