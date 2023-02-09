## R master script to run entire workflow from PReFerSim -> mssel -> Relate
## Jan 30, 2023
library(latex2exp)
setwd("/Users/vivaswatshastry/selCoefEst/PReFerSims/")
options(scipen=999)

gamma <- 100 # change Mutation rate accordingly (200 for -10, 400 for -100)

## this is for a single selection coefficient 

#----------------PReFerSim----------------#

simfile <- read.delim("simfiles/ParameterFilesConstant.txt",header=F,sep=":")
simfile[1,2] <- paste0(" ",400)
simfile[3,2] <- paste0(" ",gamma*0.5/10000)
simfile[8,2] <- paste0(" outfiles/ConstantSize",gamma,".0")
write.table(simfile,"simfiles/ParameterFilesConstant.txt",
            sep=":",row.names=F,quote=F,col.names=F)

system("GSL_RNG_SEED=100496 GSL_RNG_TYPE=mrg ../../PReFerSim/PReFerSim \\
       simfiles/ParameterFilesConstant.txt 1")

system(paste0("perl GetListOfRunsWhereFrequencyMatches.pl 0 1 \\
              outfiles/ConstantSize",gamma,".0.1.full_out.txt MiniTest/Alleles",gamma,".0.txt"))

meta.prf<-read.csv(paste0("outfiles/ConstantSize",gamma,".0.1.full_out.txt"),header=F,sep="\t")

alleles <- read.csv(paste0("MiniTest/Alleles",gamma,".0.txt"),header=F)

# randomly subset 100 alleles for tracking (since 10k is too much)
if(length(alleles[,1]>1500)){
    alleles <- sample(alleles[,1],1500)
    # write.table(alleles,paste0("MiniTest/Alleles",gamma,".0.txt"),
    #             row.names=F,quote=F,col.names=F)
}

alleles <- matrix(alleles,ncol=1)

# edit the FreqTraj file (no need to do this unless running a new gamma)

simfile <- read.delim("simfiles/ParameterFilesFreqTraj.txt",header=F,sep=":")
simfile[1,2] <- paste0(" ", 400)
simfile[3,2] <- paste0(" ",gamma*0.5/10000)
simfile[7,2] <- paste0(" outfiles/ConstantSize",gamma,".0")
simfile[8,2] <- paste0(" MiniTest/Alleles",gamma,".0.txt")
simfile[9,2] <- paste0(" MiniTest/Traj",gamma,".0.txt")
write.table(simfile,"simfiles/ParameterFilesFreqTraj.txt",
            sep=":",row.names=F,quote=F,col.names=F)

system("GSL_RNG_SEED=100496 GSL_RNG_TYPE=mrg ../../PReFerSim/PReFerSim \\
       simfiles/ParameterFilesFreqTraj.txt 1",)

meta.prf<-read.csv(paste0("outfiles/ConstantSize",gamma,".0.1.full_out.txt"),header=F,sep="\t")

#----------------mssel----------------#

setwd("~/selCoefEst/PReFerSims/")

# selsite <- rep(0,length(alleles[,1]))
cnt <- 1

## dry run
al <- alleles[sample(200,1),1]
system(paste0("sed -n '/^",al,"/p' MiniTest/Traj",gamma,".0.txt > MiniTest/TrajMini",gamma,".0.txt"))
system(paste0("perl TrajToMsselFormat.pl MiniTest/TrajMini",gamma,".0 20000 msselfiles/trajfiles/TrajMsselLike",gamma,".0.txt ",al," 0 1"))
system(paste0("cat msselfiles/trajfiles/TrajMsselLike",gamma,".0.txt | \\
                  ~/mssel/stepftn > msselfiles/trajfiles/CurrTraj",gamma,".0.txt"))
nder <- round(meta.prf[meta.prf[,5]==al,2]*200)
system(paste0("~/mssel/mssel3 ",200," 1 ",200-nder," ",nder," \\
                    msselfiles/trajfiles/CurrTraj",gamma,".0.txt 500000 -r 200 1000000 -t 200\\
                  > msselfiles/relfiles/haps",gamma,".0_",al,".ms"))
# selsite[cnt] <- as.numeric(strsplit(readLines(paste0("msselfiles/relfiles/haps",gamma,".0_",al,".ms"),n=6)[6]," ")[[1]][2])
# posn <- readLines(paste0("msselfiles/relfiles/haps",gamma,".0_",al,".ms"),n=8)[8]
# print(c(al,selsite[cnt]))
# posn <- as.numeric(strsplit(posn," ")[[1]][-1])
# selsite[cnt] <- round(posn[selsite[cnt]]*1000000)
system(paste0("sed -i '' '1,6d' msselfiles/relfiles/haps",gamma,".0_",al,".ms"))
system(paste0("echo \"// selsite: ",selsite[cnt],"\n$(cat msselfiles/relfiles/haps",gamma,".0_",al,".ms)\" > \\
                  msselfiles/relfiles/haps",gamma,".0_",al,".ms"))

## full run
for(al in alleles[,1]){
    system(paste0("sed -n '/^",al,"/p' MiniTest/Traj",gamma,".0.txt > MiniTest/TrajMini",gamma,".0.txt"))
    system(paste0("perl TrajToMsselFormat.pl MiniTest/TrajMini",gamma,".0 20000 msselfiles/trajfiles/TrajMsselLike",gamma,".0.txt ",al," 0 1"))
    system(paste0("cat msselfiles/trajfiles/TrajMsselLike",gamma,".0.txt | \\
                  ~/mssel/stepftn > msselfiles/trajfiles/CurrTraj",gamma,".0.txt"))
    nder <- round(meta.prf[meta.prf[,5]==al,2]*200)
    system(paste0("~/mssel/mssel3 ",200," 1 ",200-nder," ",nder," \\
                    msselfiles/trajfiles/CurrTraj",gamma,".0.txt 500000 -r 200 1000000 -t 200\\
                  > msselfiles/relfiles/haps",gamma,".0_",al,".ms"))
    # selsite[cnt] <- as.numeric(strsplit(readLines(paste0("msselfiles/relfiles/haps",gamma,".0_",al,".ms"),n=6)[6]," ")[[1]][2])
    # posn <- readLines(paste0("msselfiles/relfiles/haps",gamma,".0_",al,".ms"),n=8)[8]
    # print(c(al,selsite[cnt]))
    # posn <- as.numeric(strsplit(posn," ")[[1]][-1])
    # selsite[cnt] <- round(posn[selsite[cnt]]*1000000)
    system(paste0("sed -i '' '1,6d' msselfiles/relfiles/haps",gamma,".0_",al,".ms"))
    system(paste0("echo \"// selsite: ",selsite[cnt],"\n$(cat msselfiles/relfiles/haps",gamma,".0_",al,".ms)\" > \\
                  msselfiles/relfiles/haps",gamma,".0_",al,".ms"))
    cnt <- cnt+1
}

#----------------Relate----------------# 
setwd("~/selCoefEst/PReFerSims/msselfiles/relfiles")

for(al in alleles[,1]){
    # segsites <- as.numeric(strsplit(readLines(paste0("haps",gamma,".0_",al,".ms"),n=2)[2],": ")[[1]][2])
    # fm <- read.csv("fake.map",sep=" ")
    # fm[2,1] <- segsites
    # write.table(fm,"fake.map",sep=" ",row.names=F,quote=F,col.names=F)
    
    system(paste0("Rscript ../ms2haps.R haps",gamma,".0_",al,".ms \\
                  infiles/rel",gamma,".0_",al," ",1000000))
    
    system(paste0("~/relate_v1.1.9_MacOSX_M1/bin/Relate --mode All -m 1e-8 \\
                  -N 20000 --map fake.map --haps infiles/rel",gamma,".0_",al,".haps \\
                  --sample infiles/rel",gamma,".0_",al,".sample -o out",gamma,".0_",al),ignore.stdout=T,ignore.stderr=T)
    system(paste0("mv out",gamma,"* outfiles/"))
}

setwd("~/selCoefEst/PReFerSims/msselfiles/gevafiles")

#----------------GEVA#----------------#

cnt <- 1
for(al in alleles[,1]){
    ##system(paste0("python3 ../ms2vcf.py ../relfiles/haps",gamma,".0_",al,".ms 200 1000000"))
    system(paste0("Rscript ../ms2vcf2.R ../relfiles/haps",gamma,".0_",al,".ms 1000000"))
    system(paste0("~/geva/geva_v1beta --vcf infiles/haps",gamma,".0_",al,".vcf \\
                  --out infiles/ingeva",gamma,".0_",al),ignore.stdout=T)
    system(paste0("~/geva/geva_v1beta -i infiles/ingeva",gamma,".0_",al,".bin -o outfiles/out",gamma,".0_",al," --position 500000 \\
    --Ne 10000 --mut 1e-8 --hmm ~/geva/hmm/hmm_initial_probs.txt ~/geva/hmm/hmm_emission_probs.txt"),ignore.stdout=T)
    system(paste0("Rscript ~/geva/estimate.R outfiles/out",gamma,".0_",al,".pairs.txt 10000"),ignore.stdout=T)
    cnt <- cnt+1
}

#----------------processing of results----------------#
# getting ages from Relate & GEVA output

setwd("~/selCoefEst/PReFerSims/msselfiles/")

rel.est <- matrix(0,length(alleles[,1]),2)
geva.est <- matrix(0,length(alleles[,1]),3)
rel.true <- rep(0,length(alleles[,1]))

cnt<-1
for(al in alleles[,1]){
    rel.out <- read.table(paste0("relfiles/outfiles/out",gamma,".0_",al,".mut"),header=T,sep=';')
    rel.est[cnt,1] <- rel.out$age_begin[rel.out$pos_of_snp==500000]
    rel.est[cnt,2] <- rel.out$age_end[rel.out$pos_of_snp==500000]
     
    geva.est[cnt,] <- tryCatch(
        read.table(paste0("gevafiles/outfiles/out",gamma,".0_",al,".sites2.txt"),header=T,sep=' ')[,5],
        error = function(e) NA)
    
    rel.true[cnt] <- 80000+1-meta.prf[meta.prf[,5]==al,4]
    cnt<-cnt+1
}
 
plot(rel.true,0.5*(rel.est[,1]+rel.est[,2]),frame=F,log='xy',main=TeX('$\\gamma=-100$'),
     ylab='est. age',xlab='true age (PReFerSim)',col='grey80',pch=20)
segments(x0=rel.true,x1=rel.true,y0=rel.est[,1]+1,y1=rel.est[,2],col='grey80')
points(rel.true,pmin(0.5*(rel.est[,1]+rel.est[,2]),2*rel.est[,1]+1),col='black',pch=20)
# points(rel.true,geva.est[,2],col='blue',pch=20)
abline(0,1,col='grey20')
# abline(485,0.55,col='grey80',lty=2)
legend("topleft",col=c('grey','black'),bty='n',
       pch=20,legend=c(TeX('$T_{mid}$'),TeX('$min T_{mid},2T_{low}$')))

## plotting mechanism in which you sort the true ages and plot range on x-axis & (lb,ub) on y-axis
plot(1:length(alleles[,1]),sort(rel.true),col='red',pch=4,frame=F,log='y',main=TeX('$\\gamma=-10$'),
     ylab='age (gens)',xlab='rank (increasing)')
segments(x0=1:length(alleles[,1]),x1=1:length(alleles[,1]),
         y0=rel.est[order(rel.true),1]+1,y1=rel.est[order(rel.true),2],
         col='grey80')

plot(geva.est[,3],pmin(0.5*(rel.est[,1]+rel.est[,2]),2*rel.est[,1]),log='xy',
     col='black',pch=19,xlab='GEVA est.',ylab='Relate est.',frame=F)
text(600,13000,paste0('cor = ',signif(cor(
    geva.est[,3],pmin(0.5*(rel.est[,1]+rel.est[,2]),2*rel.est[,1]),use='complete.obs',method='spearman'),2)))
abline(0,1,col='grey80')
