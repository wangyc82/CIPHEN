# this procedure is to perform CPIHN
# the training data came from CHEMBL all target data
# prepare the input data
## Acd (adjacent matrix representing the interactions between new compounds and known ligands)
## Adp (adjacent matrix representing the interactions between known ligands and their receptors)
## App (adjacent matrix representing PPI interactions)

# take example data for instance
# example data includes the MACCs fingerprints for new compounds and known ligands, MACCs.fp.mat.chembl and new.chem.MACCs.mat
# adjacent matrix representing ligands and their recptors,comp_tar_mat2
# adjacent matrix representing PPI interaction, PPI_adj_mat

# preparing Acd
# sim.MACCs.mat.new<-matrix(0,nrow(new.chem.MACCs.mat),nrow(MACCs.fp.mat.chembl))
# rownames(sim.MACCs.mat.new)<-rownames(new.chem.MACCs.mat)
# colnames(sim.MACCs.mat.new)<-rownames(MACCs.fp.mat.chembl)
# for (i in 1:nrow(new.chem.MACCs.mat)) {for (j in 1:nrow(MACCs.fp.mat.chembl)) {sim.MACCs.mat.new[i,j]<-fpSim(new.chem.MACCs.mat[i,], MACCs.fp.mat.chembl[j,], method="Tanimoto")}}
# Acd.full<-matrix(0,nrow(sim.MACCs.mat.new),nrow(MACCs.fp.mat.chembl))
# rownames(Acd.full)<-rownames(sim.MACCs.mat.new)
# colnames(Acd.full)<-colnames(sim.MACCs.mat.new)
# Acd.full[which(sim.MACCs.mat.new>=0.3)]<-1, threshold 0.3 determined by the distribution of as.vector(sim.MACCs.mat.new) and make sure length(which(rowSums(Acd.full)==0))==0
# Acd<-Acd.full[,which(colSums(Acd.full)!=0)]

# prepare Adp
# Adp.full<-comp_tar_mat2[colnames(Acd),intersect(colnames(comp_tar_mat2),rownames(PPI_adj_mat))]
# Adp<-Adp.full[which(rowSums(Adp.full)!=0),]
# Adp<-Adp[,which(colSums(Adp)!=0)]
# Acd<-Acd[,rownames(Adp)]

# prepare App
# App.full<-PPI_adj_mat[colnames(Adp),]
# App<-App.full[,-which(is.na(colnames(App.full)))]
# App<-App[,which(colSums(App)!=0)]
# make sure compound and protein name has no any delimeter

# set len1, len2, len3, len4, len5 accoding the dimension of Acd, Adp, App
# for this example case, we set len1=10000;len2=100000;len3=1e6;len4=1e+7;len5=1e+8;
CIPHEN<-function(Acd, Adp, App,len1,len2,len3,len4,len5) {
  Cd<-lapply(1:nrow(Acd), function(x) colnames(Acd)[which(Acd[x,]==1)])
  names(Cd)<-rownames(Acd)
  Dc<-lapply(1:ncol(Acd), function(x) rownames(Acd)[which(Acd[,x]==1)])
  names(Dc)<-colnames(Acd)
  Dp<-lapply(1:nrow(Adp), function(x) colnames(Adp)[which(Adp[x,]==1)])
  names(Dp)<-rownames(Adp)
  P<-lapply(1:ncol(Adp), function(x) rownames(Adp)[which(Adp[,x]==1)])
  names(P)<-colnames(Adp)
  P1<-lapply(1:nrow(App), function(x) colnames(App)[which(App[x,]==1)])
  names(P1)<-rownames(App)
  P2<-lapply(1:ncol(App), function(x) rownames(App)[which(App[,x]==1)])
  names(P2)<-colnames(App)
  
  #then generate the random walk based on each meta-path
  #for C-D-C meth-path
  #set len1 is a pre-defined parameter that to determine how many random walks will be generated based on this meta-path
  I=rownames(Acd)
  L1<-character()
  for (i in 1:len1) {
    p<-sample(I,1)
    m<-Cd[[p]]
    n<-sample(m,1)
    s<-Dc[[n]]
    L1[i]<-paste(p,n,sample(s,1),sep = " ")
    rm(p,m,n,s)
  }
  
  #for D-P-D meth-path
  #set len2 is a pre-defined parameter that to determine how many random walks will be generated based on this meta-path
  I=rownames(Adp)
  L2<-character()
  for (i in 1:len2) {
    p<-sample(I,1)
    m<-Dp[[p]]
    n<-sample(m,1)
    s<-P[[n]]
    L2[i]<-paste(p,n,sample(s,1),sep = " ")
    rm(p,m,n,s)
  }
  
  #for C-D-P meth-path
  #set len3 is a pre-defined parameter that to determine how many random walks will be generated based on this meta-path
  I=rownames(Acd)
  L3<-character()
  for (i in 1:len3) {
    p<-sample(I,1)
    m<-Cd[[p]]
    n<-sample(m,1)
    s<-Dp[[n]]
    L3[i]<-paste(p,n,sample(s,1),sep = " ")
    rm(p,m,n,s)
  }
  
  #if adding PPI network, for P-P-P meth-path
  #set len4 is a pre-defined parameter that to determine how many random walks will be generated based on this meta-path
  I=rownames(App)
  L4<-character()
  for (i in 1:len4) {
    p<-sample(I,1)
    m<-P1[[p]]
    n<-sample(m,1)
    s<-P2[[n]]
    L4[i]<-paste(p,n,sample(s,1),sep = " ")
    rm(p,m,n,s)
  }
  
  #if adding PPI network, for C-D-P-P-D-C meth-path
  #set len5 is a pre-defined parameter that to determine how many random walks will be generated based on this meta-path
  I=rownames(Acd)
  L5<-character()
  for (i in 1:len5) {
    p<-sample(I,1)
    m<-Cd[[p]]
    n<-sample(m,1)
    s<-Dp[[n]]
    o<-sample(s,1)
    w<-P1[[o]]
    e<-sample(w,1)
    f<-P2[[e]]
    a<-sample(f,1)
    b<-P[[a]]
    c<-sample(b,1)
    d<-Dc[[c]]
    L5[i]<-paste(p,n,o,e,a,c,sample(d,1),sep = " ")
  }
  library("word2vec")
  model <- word2vec(x = c(L1,L2,L3,L4,L5), type = "skip-gram", dim = 50, iter = 20)
  embedding <- as.matrix(model)
  labC<-match(rownames(Acd),rownames(embedding))
  labD<-match(rownames(Adp),rownames(embedding))
  labP<-match(colnames(Adp),rownames(embedding))
  labP1<-match(colnames(App),rownames(embedding))
  Xcomp<-embedding[labC,]
  Xdrug<-embedding[labD,]
  Xtar<-embedding[labP,]
  Xpro<-embedding[labP1,]
  A<-which(Adp==1,arr.ind = TRUE)
  B<-which(Adp==0,arr.ind = TRUE)
  Xpos<-cbind(Xdrug[A[,1],],Xtar[A[,2],])
  Xneg<-cbind(Xdrug[B[sample(nrow(B),nrow(A)),1],],Xtar[B[sample(nrow(B),nrow(A)),2],])
  Xrn<-rbind(Xpos,Xneg)
  Yrn<-rep(c(1,0),c(nrow(A),nrow(A)))
  Ast<-matrix(1,nrow(Acd),ncol(App))
  Cst<-which(Ast==1,arr.ind = TRUE)
  Xst<-cbind(Xcomp[Cst[,1],],Xpro[Cst[,2],])
  library(randomForest)
  model.RF<-randomForest(Xrn,as.factor(Yrn),xtest=Xst)
  prob_test_RF<-model.RF$test$votes[,2]
  prdnewdrug<-matrix(prob_test_RF,ncol=ncol(Ast),byrow=FALSE)
  rownames(prdnewdrug)<-rownames(Acd)
  colnames(prdnewdrug)<-colnames(App)
  return(prdnewdrug)
}
