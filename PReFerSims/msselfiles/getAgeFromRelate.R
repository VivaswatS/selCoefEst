## R script to plot estimated ages across the genome using Relate
## haplotypes simulated using mssel 
## Jan 30, 2022

library(latex2exp)

setwd("/Users/vivaswatshastry/selCoefEst/PReFerSims/msselfiles")

# read in the third line
posn <- readLines("haps10.0.ms",n=3)[3]

# remove the first field 
posn <- as.numeric(strsplit(posn," ")[[1]][-1])

# read in estimated ages from the Relate output
rel.out <- read.table("example.mut",header=T,sep=';')

# getting an estimate of age by looking at SNPs around the selected region
(reg.age<-pmin(0.5*(rel.out$age_begin+rel.out$age_end),
     2*rel.out$age_begin)[posn[rel.out$pos_of_snp]>0.495
                          & posn[rel.out$pos_of_snp]<0.505])

mean(reg.age)

# plot the results from Relate
plot(posn[rel.out$pos_of_snp],0.5*(rel.out$age_begin+rel.out$age_end),
     ylab='est. age (gens)',xlab='pos (in mssel units)',ylim=c(0,20000),
     main=TeX('$\\rho=200,nsites=20'),col='grey',pch=8,xlim=c(0.49,0.51))
points(posn[rel.out$pos_of_snp],pmin(0.5*(rel.out$age_begin+rel.out$age_end),2*rel.out$age_begin),
     col='red',pch=20)
abline(v=0.5,col='grey40',lty=2)
abline(h=14900,col='red',lty=5)
legend("topright",pch=c(8,20),legend=c(TeX('$T_{mid}$'),TeX('$min T_{mid},2T_{low}$')))


# how does the Relate estimate compare to the heuristic from Simons et al 2022?
plot(0.5*(rel.out$age_begin+rel.out$age_end),
     pmin(0.5*(rel.out$age_begin+rel.out$age_end),2*rel.out$age_begin),
     col='grey80',pch=20)
text(0.5*(rel.out$age_begin+rel.out$age_end),
     pmin(0.5*(rel.out$age_begin+rel.out$age_end),2*rel.out$age_begin),
     col='grey80',labels=rel.out$pos_of_snp)
abline(0,1)

nt<-read.tree('/Users/vivaswatshastry/selCoefEst/PReFerSims/msselfiles/relfiles/tempfiles/sim100.0_15980626.newick')
bid<-243
start.end<-matrix(NA,nrow=100,ncol=2)
for(i in 1:100){
    start.end[i,1]<-max(node.depth.edgelength(nt[[i]]))-node.depth.edgelength(nt[[i]])[bid[2]]
    start.end[i,2]<-max(node.depth.edgelength(nt[[i]]))-node.depth.edgelength(nt[[i]])[bid[1]]
}

nddep <- rep(0,399)
for(n in 202:399){
        nddep[n] <- max(node.depth.edgelength(nt))-node.depth.edgelength(nt)[n]
}

brlen <- rep(0,398)
for(ie in 1:398){
        brlen[ie]<-node.depth.edgelength(nt)[nt$edge[ie,2]]-node.depth.edgelength(nt)[nt$edge[ie,1]]
}
