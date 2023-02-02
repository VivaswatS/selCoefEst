## R master script to run entire workflow from PReFerSim -> mssel -> Relate
## Jan 30, 2023

setwd("/Users/vivaswatshastry/selCoefEst/PReFerSims/")

gamma <- 0

## this is for a single selection coefficient 

#----------------PReFerSim----------------#

simfile <- read.delim("simfiles/ParameterFilesConstant.txt",header=F,sep=":")
simfile[3,2] <- paste0(" ",gamma*0.5/10000)
simfile[8,2] <- paste0("outfiles/ConstantSize",gamma,".0")
write.table(simfile,"simfiles/ParameterFilesConstant.txt",
            sep=":",row.names=F,quote=F,col.names=F)

system("GSL_RNG_SEED=100496 GSL_RNG_TYPE=mrg ../../PReFerSim/PReFerSim \\
       simfiles/ParameterFilesConstant.txt 1")

system(paste0("perl GetListOfRunsWhereFrequencyMatches.pl 0 1 \\
              outfiles/ConstantSize",gamma,".0.1.full_out.txt MiniTest/Alleles",gamma,".0.txt"))

alleles <- read.csv(paste0("MiniTest/Alleles",gamma,".0.txt"),header=F)

# randomly subset 100 alleles for tracking (since 10k is too much)
if(length(alleles[,1]>100)){
    alleles <- sample(alleles[,1],100)
    write.table(alleles,paste0("MiniTest/Alleles",gamma,".0.txt"),
                row.names=F,quote=F,col.names=F)
}

alleles <- read.csv(paste0("MiniTest/Alleles",gamma,".0.txt"),header=F)

# edit the FreqTraj file

simfile <- read.delim("simfiles/ParameterFilesFreqTraj.txt",header=F,sep=":")
simfile[3,2] <- paste0(" ",gamma*0.5/10000)
simfile[7,2] <- paste0(" outfiles/ConstantSize",gamma,".0")
simfile[8,2] <- paste0(" MiniTest/Alleles",gamma,".0.txt")
simfile[9,2] <- paste0(" MiniTest/Traj",gamma,".0.txt")
write.table(simfile,"simfiles/ParameterFilesFreqTraj.txt",
            sep=":",row.names=F,quote=F,col.names=F)

system("GSL_RNG_SEED=100496 GSL_RNG_TYPE=mrg ../../PReFerSim/PReFerSim \\
       simfiles/ParameterFilesFreqTraj.txt 1")

meta.prf<-read.csv(paste0("outfiles/ConstantSize",gamma,".0.1.full_out.txt"),header=F,sep="\t")

#----------------mssel----------------#

selsite <- rep(0,length(alleles[,1]))
cnt <- 1

for(al in alleles[,1]){
    # system(paste0("perl TrajToMsselFormat.pl MiniTest/Traj",gamma,".0 20000 msselfiles/trajfiles/TrajMsselLike",gamma,".0.txt ",al," 0 1"))
    # system(paste0("cat msselfiles/trajfiles/TrajMsselLike",gamma,".0.txt | \\
    #               ~/mssel/stepftn > msselfiles/trajfiles/CurrTraj",gamma,".0_",al,".txt"))
    # nder <- round(meta.prf[meta.prf[,5]==al,2]*200)
    # system(paste0("~/mssel/mssel3 ",200," 1 ",200-nder," ",nder," \\
    #                 msselfiles/trajfiles/CurrTraj",gamma,".0_",al,".txt 500000 -r 200 1000000 -t 200\\
    #               > msselfiles/relfiles/haps",gamma,".0_",al,".ms"))
    # selsite[cnt] <- as.numeric(strsplit(readLines(paste0("msselfiles/relfiles/haps",gamma,".0_",al,".ms"),n=6)[6]," ")[[1]][2])
    # posn <- readLines(paste0("msselfiles/relfiles/haps",gamma,".0_",al,".ms"),n=3)[3]
    # posn <- as.numeric(strsplit(posn," ")[[1]][-1])
    # selsite[cnt] <- round(posn[473]*1000000)
    # system(paste0("sed -i '' '1,6d' msselfiles/relfiles/haps",gamma,".0_",al,".ms"))
    # system(paste0("echo \"// selsite: ",selsite[cnt],"\n$(cat msselfiles/relfiles/haps",gamma,".0_",al,".ms)\" > \\
    #               msselfiles/relfiles/haps",gamma,".0_",al,".ms"))
    # cnt <- cnt+1
}

setwd("./msselfiles/relfiles")

#----------------Relate----------------# 

for(al in alleles[,1]){
    # segsites <- as.numeric(strsplit(readLines(paste0("haps",gamma,".0_",al,".ms"),n=2)[2],": ")[[1]][2])
    # fm <- read.csv("fake.map",sep=" ")
    # fm[2,1] <- segsites
    # write.table(fm,"fake.map",sep=" ",row.names=F,quote=F,col.names=F)
    
    system(paste0("Rscript ../ms2haps.R haps",gamma,".0_",al,".ms \\
                  rel",gamma,".0_",al," ",1000000))
    
    system(paste0("~/relate_v1.1.9_MacOSX_M1/bin/Relate --mode All -m 1e-8 \\
                  -N 20000 --map fake.map --haps infiles/rel",gamma,".0_",al,".haps \\
                  --sample infiles/rel",gamma,".0_",al,".sample -o outfiles/out",gamma,".0_",al))
}

setwd("../gevafiles")

#----------------GEVA#----------------#

cnt <- 1
for(al in alleles[1:2,1]){
    ##system(paste0("python3 ../ms2vcf.py ../relfiles/haps",gamma,".0_",al,".ms 200 1000000"))
    system(paste0("Rscript ../ms2vcf2.R ../relfiles/haps",gamma,".0_",al,".ms 1000000"))
    system(paste0("~/geva/geva_v1beta --vcf infiles/haps",gamma,".0_",al,".vcf \\
                  --out infiles/ingeva",gamma,".0_",al))
    system(paste0("~/geva/geva_v1beta -i infiles/ingeva",gamma,".0_",al,".bin -o outfiles/out",gamma,".0_",al," --position ",selsite[cnt],"\\
    --Ne 10000 --mut 1e-8 --hmm ~/geva/hmm/hmm_initial_probs.txt ~/geva/hmm/hmm_emission_probs.txt"))
    system(paste0("Rscript ~/geva/estimate.R outfiles/out",gamma,".0_",al,".pairs.txt 10000"))
    cnt <- cnt+1
}

#----------------processing of results----------------#
# getting ages from relate output

setwd("../relfiles/outfiles")

cnt<-1
for(al in alleles[,1]){
    rel.out <- read.table(paste0(""),header=T,sep=';')
}

