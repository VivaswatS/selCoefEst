options(scipen=999)
# Assumes that input is of form
#  //(some text)
#  segsites: (num_segsites)
#  position: (positions)
#  10101110 etc

#Input filename
infile <- commandArgs(trailingOnly = T)[1]
outfile <- paste0("infiles/",tail(strsplit(strsplit(infile,'.ms')[[1]][1],'/')[[1]],1),".vcf")
#if nsites argument is specified, multiply this to pos, otherwise multiply 1
if(length(commandArgs(trailingOnly = T)) == 2){
    nsites       <- as.numeric(commandArgs(trailingOnly = T)[2])
}else if(length(commandArgs(trailingOnly = T)) > 2 || length(commandArgs(trailingOnly = T)) < 2){
    
    cat("###################################################\n")
    cat("Usage: Rscript ms2vcf.R infile.ms nsite\n\n")
    cat("infile.ms: Input filename with file extension.\n")
    cat("outfile: Output filename withouth file extension.\n")
    cat("nsites: (Optional) Number of simulated sites. Default value is 1. This is multiplied to the positions.\n")
    cat("###################################################\n")
    
    quit(status=1)
}else{
    nsites <- 1
}

# read file as vector of strings
infile_lines <- readLines(infile)
startlines   <- which(startsWith(infile_lines, prefix = "//"))

if(length(startlines) == 0){
    cat("No input in ms format.\n")
    
    quit(status = 1)
}

pos  <- unlist(strsplit(infile_lines[startlines+2], split = " "))
pos  <- round(as.numeric(as.matrix(pos[-1])) * nsites)

seq  <- infile_lines[(startlines+3):length(infile_lines)]
seq  <- seq[seq!=""] #remove empty lines at the end
# convert seq to a matrix
seq  <- t(matrix(as.numeric(unlist(strsplit(seq, split = ""))), ncol = length(pos), byrow = T))

# find sites with multiple mutations
dupl <- which(duplicated(pos))
if(length(dupl) > 0){
    dupl <- sort(unique(c(dupl, dupl-1)))
    pos  <- pos[-dupl]
    seq  <- seq[-dupl,]
}
if(length(pos) == 0){
    cat("No segsites\nBP have to be integers! (Use third argument)")
    quit(status = 1)
}

N <- dim(seq)[2]/2
L <- dim(seq)[1]

vcf.seq <- matrix(0,L,N)
for(i in 1:L){
    for(j in 1:N){
        vcf.seq[i,j] <- paste0(seq[i,2*j-1],"|",seq[i,2*j])
    }
}

##### write outfile.haps #####
fileConn<-file(outfile)
writeLines(c('##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="All filters passed">\n##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##contig=<ID=1,length=2147483647>',
             paste("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
                   paste0("UNR",1:N,collapse='    '),sep='\t')), fileConn)
close(fileConn)
haps <- cbind(rep(1,L), pos, paste("SNP", 1:L, sep = ""), "A", "T", ".", "PASS", ".", "GT", vcf.seq)
write.table(haps, outfile, sep='\t', quote = F, row.names = F, col.names = F, append=T)


