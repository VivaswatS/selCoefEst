## R master script to run entire workflow from PReFerSim -> mssel -> Relate
## Jan 30, 2023
library(latex2exp)
# setwd("/Users/vivaswatshastry/selCoefEst/PReFerSims/")
options(scipen=999)

gamma <- -10 # change theta accordingly (200 for -10, 400 for -100)
theta <- 200
    
## this is for a single selection coefficient 

#----------------PReFerSim----------------#

simfile <- read.delim("simfiles/ParameterFilesConstant.txt",header=F,sep=":")
simfile[1,2] <- paste0(" ",theta)
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
simfile[1,2] <- paste0(" ", theta)
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
system(paste0("perl TrajToMsselFormat.pl MiniTest/TrajMini",gamma,".0 20000 msselfiles/trajfiles/TrajMsselLike",gamma,".0.txt 1 0 1"))
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

#----------------GEVA----------------#

cnt <- 1
for(al in alleles[,1]){
    ##system(paste0("python3 ../ms2vcf.py ../relfiles/haps",gamma,".0_",al,".ms 200 1000000"))
    system(paste0("Rscript ../ms2vcf2.R ../relfiles/haps",gamma,".0_",al,".ms 1000000"))
    system(paste0("~/geva/geva_v1beta --vcf infiles/haps",gamma,".0_",al,".vcf \\
                  --out infiles/ingeva",gamma,".0_",al))
    system(paste0("~/geva/geva_v1beta -i infiles/ingeva",gamma,".0_",al,".bin -o outfiles/out",gamma,".0_",al," --position 500000 \\
    --Ne 10000 --mut 1e-8 --hmm ~/geva/hmm/hmm_initial_probs.txt ~/geva/hmm/hmm_emission_probs.txt"))
    system(paste0("Rscript ~/geva/estimate.R outfiles/out",gamma,".0_",al,".pairs.txt 10000"))
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
     
    # geva.est[cnt,] <- tryCatch(
    #     read.table(paste0("gevafiles/outfiles/out",gamma,".0_",al,".sites2.txt"),header=T,sep=' ')[,5],
    #     error = function(e) NA)
    
    rel.true[cnt] <- 80000+1-meta.prf[meta.prf[,5]==al,4]
    cnt<-cnt+1
}

## reading in RelConstantSize files from midway3
df<-read.table('msselfiles/relfiles/sumfiles/RelConstantSize0_id20.9.full_out.txt',header=F,sep='\t')
df <- rbind(df,read.table('msselfiles/relfiles/sumfiles/RelConstantSize100_id21.9.full_out.txt',header=F,sep='\t'))
rel.true <- 80000+1-df[,5]
rel.est <- matrix(cbind(df[,7],df[,8]),ncol=2)
 
plot(rel.true,0.5*(rel.est[,1]+rel.est[,2]),log='xy',frame=F,main=TeX('$\\gamma=-1$'),
     ylab='est. age (Relate)',xlab='true age (PReFerSim)',col='grey80',pch=20,lwd=4,ylim=c(1,100000))
segments(x0=rel.true,x1=rel.true,y0=rel.est[,1]+1,y1=rel.est[,2],col='grey80',lwd=3)
points(rel.true,pmin(0.5*(rel.est[,1]+rel.est[,2]),2*rel.est[,1]+1),col='black',pch=20)
# points(rel.true,geva.est[,2],col='blue',pch=20)
# abline(lsfit(rel.true,0.5*(rel.est[,1]+rel.est[,2])),col='grey60',lty=2)
abline(lsfit(log10(rel.true),log10(0.5*(rel.est[,1]+rel.est[,2]))),col='grey60',lty=2,lwd=3)
abline(0,1,col='grey20',lwd=3)
text(5,4000,TeX('$y \\approx 0.06 x + 2300$'))
text(5,2000,TeX('$R^2 \\approx 0.06$'))
legend("topleft",col=c('grey','black',"blue"),bty='n',
       pch=20,legend=c(TeX('$T_{mid}$'),TeX('$min T_{mid},2T_{low}$'),'GEVA'))

# -100
# 1.5174 + 0.3649x, (no log) 169.6137567 + 0.8649518x, resid = 146
# 0
# 0.7724 + 0.7765x, (no log) 2005.3094596 + 0.6372335x, resid = 156

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

plot(df[,3]*200,rel.true,pch=20,frame=F,col='grey50',
     xlab='sample freq.',ylab='true age (PReFerSim)',)
points(df[,3]*200,rel.true,pch=20,col='grey80')
legend('bottomright',c(TeX('$\\gamma=-100$'),TeX('$\\gamma=-1$')),col=c('grey80','grey50'),pch=20)


## reading in Newick trees output by mssel
library(ape)
# tt<-read.tree('/Users/vivaswatshastry/selCoefEst/PReFerSims/msselfiles/haps5.0_r0.newick')
tt <- read.tree(text='((5:0.123237,6:0.123237):0.866739,(3:0.322997,(1:0.221121,(2:0.058565,4:0.058565):0.162556):0.101876):0.666980);')

## plotting thresholded SFS
compute_sfs <- function(n, gamma) {
    x <- seq(1/(2*n), 1-1/(2*n), length.out=n-1)
    sfs <- (1-exp(-2*gamma*(1-x))) / (x*(1-x)*(1-exp(-2*gamma)))
    return(sfs / sum(sfs))
}

# Set parameters
n <- 10000  # sample size
gamma_values <- c(0, -1, -5, -10, -50, -100)  # selection strengths
colors <- gray.colors(length(gamma_values), start=0.2, end=0.8)

# Compute SFS for each gamma value
sfs_list <- lapply(gamma_values, function(g) compute_sfs(n, g))

# Prepare data for plotting
plot_data <- data.frame(
    frequency = rep(seq(1/(2*n), 1-1/(2*n), length.out=n-1), length(gamma_values)),
    density = unlist(sfs_list),
    gamma = factor(rep(gamma_values, each=n-1))
)

# Create the plot
ggplot(plot_data, aes(x=frequency, y=density, color=gamma)) +
    geom_line(size=2) +
    scale_color_manual(values=colors) + scale_x_log10() + scale_y_log10(limits=c(1e-10,1)) +
    labs(x="Allele Frequency", y="Density", 
         title="Expected Site Frequency Spectrum under Selection",
         color="γ", linetype="γ") +
    theme_classic() + 
    theme(legend.position="right") + 
    geom_vline(xintercept=0.025,linetype='dashed',color='dodgerblue4',size=1.5) + 
    geom_abline(slope=-1,color='black')





