## R master script to run entire workflow from PReFerSim -> mssel -> Relate
## Jan 30, 2023

setwd("/Users/vivaswatshastry/selCoefEst/PReFerSims/")

gamma <- 5

## this is for a single selection coefficient 

#----------------PReFerSim----------------#

simfile <- read.delim("simfiles/ParameterFilesConstant.txt",header=F,sep=":")
simfile[3,2] <- paste0(" ",gamma*0.5/10000)
simfile[8,2] <- paste0("outfiles/ConstantSize",gamma,".0")
write.table(simfile,"simfiles/ParameterFilesConstant.txt",
            sep=":",row.names=F,quote=F,col.names=F)

system(paste0("perl GetListOfRunsWhereFrequencyMatches.pl 0 1 outfiles/ConstantSize",gamma,".0.1.full_out.txt MiniTest/Alleles",gamma,".0.txt"))

alleles<-read.csv("MiniTest/Alleles5.0.txt",header=F)

# edit the FreqTraj file

simfile <- read.delim("simfiles/ParameterFilesFreqTraj.txt",header=F,sep=":")
simfile[3,2] <- paste0(" ",gamma*0.5/10000)
simfile[7,2] <- paste0(" outfiles/ConstantSize",gamma,".0")
simfile[8,2] <- paste0(" MiniTest/Alleles",gamma,".0.txt")
simfile[9,2] <- paste0(" MiniTest/Traj",gamma,".0.txt")
write.table(simfile,"simfiles/ParameterFilesFreqTraj.txt",
            sep=":",row.names=F,quote=F,col.names=F)

system("GSL_RNG_SEED=100496 GSL_RNG_TYPE=mrg ../../PReFerSim/PReFerSim 
       simfiles/ParameterFilesFreqTraj.txt 1")

meta.prf<-read.csv(paste0("outfiles/ConstantSize",gamma,".0.1.full_out.txt"),header=F,sep="\t")

for(al in alleles[,1]){
    system(paste0("perl TrajToMsselFormat.pl MiniTest/Traj",gamma,".0 20000 msselfiles/TrajMsselLike",gamma,".0.txt ",al," 0 1"))
    system(paste0("cat msselfiles/TrajMsselLike",gamma,".0.txt | ~/mssel/stepftn > msselfiles/CurrTraj",gamma,".0_",al,".txt"))
    nder <- round(meta.prf[meta.prf[,5]==6514609,2]*200)
    system(paste0("~/mssel/mssel3 ",200," 1 ",nder," ",200-nder," msselfiles/CurrTraj",gamma,".0_",al,".txt 0.5 -r 200 20 -t 200 > msselfiles/haps",gamma,".0_",al,".ms"))
    system(paste0("sed -i -e '1,6d' < haps",gamma,".0_",al,".ms"))
    system(paste0("echo \"//\n$(cat haps",gamma,".0_",al,".ms)\" > haps",gamma,".0_",al,".ms"))
}

setwd("./msselfiles")

for(al in alleles[,1]){
    segsites<-
}
