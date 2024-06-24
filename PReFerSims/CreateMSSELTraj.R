## Script to convert a trajectory from PReFerSim to mssel traj
## vivaswat, Jul 15 2022

args<-commandArgs(trailingOnly=TRUE)
fileName<-args[1]
outFile<-args[2]
twoN0<-20000 #numeric(args[3])

traj<-read.table(fileName, header=F)
allids<-unique(traj[,1])

# picking an allele id that is older than 20 gens
xx<-sample(allids,1)
while(nrow(traj[traj[,1]==xx,])<20){
    xx<-sample(allids,1)
}
print(paste0("Sampling allele ID: ",xx))

# writing out file to mssel format
sink(outFile)
cat("ntraj: 1\n")
cat("maxpop: 1\n")
cat(paste0("n: ",nrow(traj[traj[,1]==xx,]),"\n"))
for(gen in 1:nrow(traj[traj[,1]==xx,])){
    cat(paste0(gen/twoN0," ",traj[traj[,1]==xx,][gen,2],"\n"))
}
sink()
