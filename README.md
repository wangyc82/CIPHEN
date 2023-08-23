**CPIHN** 

**C**ompound-**D**rug **I**nteraction using **H**eterogeneous **N**etwork. 

# Status

Active development

# Introduction

Early prediction of compound-protein interactions (CPIs) is a critical step in drug discovery. Among various obstacles hindering clinical translation, lacking effective methods for multi-type network integration has become a bottleneck. CPIHN detected CPIs across the entire human protein space through integrating protein-protein interaction (PPI) network into a heterogeneous network (HN).

# Usage

1. Installation

   Prerequisites of DeepDRK includes the following: 

   - R is properly installed; 

   - Rscript is available in your system path ($PATH);

   - git (2.21.1)

    Installation of DeepDR includes the following steps:

    - step 1: git clone https://github.com/wangyc82/CPIHN;

    - step 2: download the example data (combination_data.RData) from https://wanglab.shinyapps.io/DeepDRK/, and put it in the CPIHN folder.

    Dependencies of CPIHN includes the following: 

    - word2vec package and its dependencies.
    - randomForest package and its dependencies.

    Testing of successful installation by running the following commands in R:
     
       > library(word2vec)
       > library(randomForest)


2. Preparation of the input files

Acd (adjacent matrix representing the interactions between new compounds and known ligands)
Adp (adjacent matrix representing the interactions between known ligands and their receptors)
App (adjacent matrix representing PPI interactions)

take example data for instance

     > load('~/CPIHN/example.RData')

example data includes the MACCs fingerprints for new compounds and known ligands, MACCs.fp.mat.chembl and new.chem.MACCs.mat
adjacent matrix representing ligands and their recptors,comp_tar_mat2
adjacent matrix representing PPI interaction, PPI_adj_mat

preparing Acd

     > sim.MACCs.mat.new<-matrix(0,nrow(new.chem.MACCs.mat),nrow(MACCs.fp.mat.chembl))
     > rownames(sim.MACCs.mat.new)<-rownames(new.chem.MACCs.mat)
     > colnames(sim.MACCs.mat.new)<-rownames(MACCs.fp.mat.chembl)
     > for (i in 1:nrow(new.chem.MACCs.mat)) {
     >   for (j in 1:nrow(MACCs.fp.mat.chembl)) {
     >     sim.MACCs.mat.new[i,j]<-fpSim(new.chem.MACCs.mat[i,], MACCs.fp.mat.chembl[j,], method="Tanimoto")
     >      }
     >     }

     > Acd.full<-matrix(0,nrow(sim.MACCs.mat.new),nrow(MACCs.fp.mat.chembl))
     > rownames(Acd.full)<-rownames(sim.MACCs.mat.new)
     > colnames(Acd.full)<-colnames(sim.MACCs.mat.new)
     > Acd.full[which(sim.MACCs.mat.new>=0.3)]<-1
threshold 0.3 determined by the distribution of as.vector(sim.MACCs.mat.new) and make sure length(which(rowSums(Acd.full)==0))==0
     
     > Acd<-Acd.full[,which(colSums(Acd.full)!=0)]

prepare Adp

     > Adp.full<-comp_tar_mat2[colnames(Acd),intersect(colnames(comp_tar_mat2),rownames(PPI_adj_mat))]
     > Adp<-Adp.full[which(rowSums(Adp.full)!=0),]
     > Adp<-Adp[,which(colSums(Adp)!=0)]
     > Acd<-Acd[,rownames(Adp)]

prepare App

    > App.full<-PPI_adj_mat[colnames(Adp),]
    > App<-App.full[,-which(is.na(colnames(App.full)))]
    > App<-App[,which(colSums(App)!=0)]

make sure compound and protein name has no any delimeter

set len1, len2, len3, len4, len5 accoding the dimension of Acd, Adp, App
for this example case, we set len1=10000;len2=100000;len3=1e6;len4=1e+7;len5=1e+8;

3. Running CPIHN

The main function of DeepDRK is DeepDRKpredictor.R. Get your input files prepared, and run it like this:

Usage example:

    > source('~/CPIHN/CPIHN.R')
    > predictions<-CPIHN(Acd,Adp,App)
    
In the following, we will show how to get MACCs fingerprints from compund sdf file

    > library(ChemmineR)
    > library(ChemmineOB)
    > compounds<-read.SDFset("~/CPIHN/compounds.sdf")
    > MACCs.fp.mat<-fingerprintOB(compounds,"MACCS")@fpma


# Contact

For technical issues please send an email to wangyongcui@mail.kib.ac.cn.
