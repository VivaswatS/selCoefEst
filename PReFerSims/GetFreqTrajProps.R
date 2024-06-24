args <- commandArgs(trailingOnly = TRUE)

fl <- args[1]

alls <- read.csv(fl, header=F)

lts <- numeric(length(alls[,1]))
afs <- numeric(length(alls[,1]))

for(i in 1:length(alls[,1])){
    lts[i] <- as.numeric(system(sprintf("grep '%d' %s | wc -l",alls[i,1],sprintf("MiniTest/Traj%s",gsub(".*Alleles(.+)", "\\1", fl))),intern=T))
    afs[i] <- as.numeric(strsplit(system(sprintf("grep '%d' %s | tail -n 1",alls[i,1],sprintf("MiniTest/Traj%s",gsub(".*Alleles(.+)", "\\1", fl))),intern=T),"\t")[[1]][2])
}

print("Top 5 longest trajectories: ")
cbind(alls[order(lts,decreasing=T)[1:5],1],sort(lts,decreasing=T)[1:5],afs[order(lts,decreasing=T)[1:5]])

print("Bottom 3 shortest trajectories: ")
cbind(alls[order(lts)[1:3],1],sort(lts)[1:3],afs[order(lts)[1:3]])
