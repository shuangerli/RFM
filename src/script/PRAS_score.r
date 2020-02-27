#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This script is used for calculating the PRAS score using the peak intensity and position.

## Scoring function of PRAS_exp
PRAS_exp_gene_score=function(binding_table,gene_name,d0){
    gene_name=as.character(gene_name)
    A=binding_table[,c('reads','Gene','dis')]
    if(!gene_name %in% A$Gene){return(0)}
    else{
        Ag=A[A$Gene == gene_name,]
        Ag$score=exp(-Ag$dis/d0)*Ag$reads
        S=aggregate(Ag$score,list(Ag$Gene),sum)
        return(S$x)
    }
}




# Five inputs: 
# 1. Annotated binding table
# 2. Interested gene list
# 3. 5/3 from 5' end or 3' end is the reference site (based on RBP function)
# 4. d0
# 5. PRAS scores output filename

# read binding table
B=read.delim(args[1],header=T)
# read interested gene list
gene_ls = read.delim(args[2],header = F)
gene_ls = gene_ls[,1]
# get corresponding columns considering reference site.
if(args[3]=='5') {
  B=B[,c('reads','end5_nearest_gene','dis5')]
} else if (args[3]=='3') {
  B=B[,c('reads','end3_nearest_gene','dis3')]
}
colnames(B)=c('reads','Gene','dis')
# calculate score of genes using PRAS_exp
PRAS_exp_score = unlist(lapply(gene_ls,FUN = PRAS_exp_gene_score,binding_table=B,d0=as.integer(args[4])))
# generate PRAS table
PRAS_tab = data.frame(cbind(data.frame(gene_ls),PRAS_exp_score))
colnames(PRAS_tab)=c("Gene","PRAS")
PRAS_tab$log2PRAS = log2(PRAS_tab$PRAS + 1)
PRAS_tab$rank = rank(-PRAS_tab$PRAS)
PRAS_tab = PRAS_tab[order(PRAS_tab$rank),]

write.table(PRAS_tab,file=args[5],quote=F,col.names=T,row.names=F,sep="\t")

# Command:
# Rscript --vanilla PRAS_score.r binding_tab gene_list 5/3 d0 PRAS_tab

