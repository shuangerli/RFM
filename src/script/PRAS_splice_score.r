#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This script is used for calculating the PRAS splice score using the peak intensity and position.

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
# 1. Annotated binding table for exons and introns
# 2. Interested gene list
# 4. d0
# 5. PRAS scores output filename

# read binding table
B_exon=read.delim(args[1],header=T)
B_intron = read.delim(args[2],header=T)
# read interested gene list
gene_ls = read.delim(args[3],header = F)
gene_ls = gene_ls[,1]
# get corresponding columns considering reference site.
## 5' splice site
B_5_exon = B_exon[,c('reads','end3_nearest_gene','dis3')]
B_5_intron = B_intron[,c('reads','end5_nearest_gene','dis5')]
B_3_exon = B_exon[,c('reads','end5_nearest_gene','dis5')]
B_3_intron = B_intron[,c('reads','end3_nearest_gene','dis3')]

colnames(B_5_exon)=colnames(B_5_intron)=colnames(B_3_exon)=colnames(B_3_intron)=c('reads','Gene','dis')
# calculate score of genes using PRAS_exp
PRAS_exp_5_exon_score = unlist(lapply(gene_ls,FUN = PRAS_exp_gene_score,binding_table=B_5_exon,d0=as.integer(args[4])))
PRAS_exp_5_intron_score = unlist(lapply(gene_ls,FUN = PRAS_exp_gene_score,binding_table=B_5_intron,d0=as.integer(args[4])))
PRAS_exp_3_exon_score = unlist(lapply(gene_ls,FUN = PRAS_exp_gene_score,binding_table=B_3_exon,d0=as.integer(args[4])))
PRAS_exp_3_intron_score = unlist(lapply(gene_ls,FUN = PRAS_exp_gene_score,binding_table=B_3_intron,d0=as.integer(args[4])))
# generate PRAS table
PRAS_tab = data.frame(cbind(data.frame(gene_ls),PRAS_exp_5_exon_score,PRAS_exp_5_intron_score,PRAS_exp_3_exon_score,PRAS_exp_3_intron_score))
colnames(PRAS_tab)=c("Gene","PRAS_5_SS_exon","PRAS_5_SS_intron","PRAS_3_SS_exon","PRAS_3_SS_intron")
PRAS_tab$PRAS_5_SS = PRAS_tab$PRAS_5_SS_exon + PRAS_tab$PRAS_5_SS_intron
PRAS_tab$PRAS_3_SS = PRAS_tab$PRAS_3_SS_exon + PRAS_tab$PRAS_3_SS_intron
PRAS_tab$PRAS_SS_exon = PRAS_tab$PRAS_5_SS_exon + PRAS_tab$PRAS_3_SS_exon
PRAS_tab$PRAS_SS_intron = PRAS_tab$PRAS_5_SS_intron + PRAS_tab$PRAS_3_SS_intron
PRAS_tab$PRAS_SS = PRAS_tab$PRAS_SS_exon + PRAS_tab$PRAS_SS_intron
clnms = colnames(PRAS_tab)

PRAS_log2_tab = log2(PRAS_tab[,2:10]+1)
PRAS_tab = cbind(PRAS_tab,PRAS_log2_tab)
PRAS_tab$rank = rank(-PRAS_tab$PRAS_SS)
PRAS_tab = PRAS_tab[order(PRAS_tab$rank),]

colnames(PRAS_tab)[11:19]=paste("log2",clnms[2:10],sep="")
write.table(PRAS_tab,file=args[5],quote=F,col.names=T,row.names=F,sep="\t")

# Command:
# Rscript --vanilla PRAS_score.r binding_tab_exon binding_tab_intron gene_list d0 PRAS_tab

