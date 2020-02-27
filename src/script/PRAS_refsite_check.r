#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This script is used for calculating the binding signal accumulation around the potential reference sites. PRAS will plot out the binding profiles around the 5' end and 3' end of the selected genomic region and let user decide which one to use.


# Four inputs:
# 1. Annotated binding tables in the order of 5'UTR, CDS, 3'UTR.
# 2. Interested gene list
# 3. Range of distance for examine the profiles
# 4. output filename.

# read binding table
B_5utr = read.delim(args[1],header=T)
B_cds = read.delim(args[2],header=T)
B_3utr = read.delim(args[3],header=T)

# read interested gene list
gene_ls = read.delim(args[4],header = F)
gene_ls = gene_ls[,1]
n = NROW(gene_ls)

# define accumulate reads function
row_contr = function(x,range_){
    v = rep(0,length(range_))
    hw = as.integer(x[1]/2)
    avg_ = x[2]/(x[1]+1)
    s = max(0,x[3] - hw)
    if (s> max(range_)){
        return(v)
    }else {
        e = min(max(range_)+1,x[3] + hw)
        v[s:e] = avg_
        return(v)
    }
}

# get distance matrix
get_dis = function(B,gene_ls,end_,l){
    # determine whether is from 5'end or 3'end
    dis = paste("dis",end_,sep="")
    # get peak length
    B$peak_len = B$peak_end - B$peak_start
    # get corresponding columns considering reference site.
    B=B[B$end5_nearest_gene %in% gene_ls,c('peak_len','reads',dis)]

    # generate distance range to examine.
    x = as.matrix(seq(0,l,by=1))
    colnames(x) = "dis"
    # generate contribution of each row
    apldf = apply(B,1,row_contr,range_ = seq(0,l,by=1))
    B_pos = as.matrix(rowSums(apldf))
    B_prof = as.data.frame(cbind(x,B_pos))
    B_prof[is.na(B_prof)]=0
    colnames(B_prof)[2]="value"
    return(B_prof)
}


l = as.numeric(args[5])
# collect binding profiles
ref1_down = get_dis(B_5utr,gene_ls,5,l)
ref2_up = get_dis(B_5utr,gene_ls,3,l)
ref2_down = get_dis(B_cds,gene_ls,5,l)
ref3_up = get_dis(B_cds,gene_ls,3,l)
ref3_down = get_dis(B_3utr,gene_ls,5,l)
ref4_up = get_dis(B_3utr,gene_ls,3,l)
# get the largest value of the intensities within the defined region
ymax = max(c(max(ref1_down$value),max(ref2_up$value),max(ref2_down$value),max(ref3_up$value),max(ref3_down$value),max(ref4_up$value)))/n


# plot the upstream and downstream profile
up_down = function(up,down,l,ymax,lb){
    plot(-up[,1],up[,2]/n,xlab = "Distance",ylab="",yaxt = 'n',xaxt = 'n',type='l',lwd=2,col="blue",main = "",xlim=c(-l,l),ylim=c(0,ymax))
    lines(down[,1],down[,2]/n,col='blue',lwd=2)
    axi_s = round(seq(-l,l,length.out=5))
    axis(1,at = axi_s,labels=c(axi_s[1:2],lb,axi_s[4:5]))
}

# plot out the profiles
pref = args[6]
ref1_file = paste(pref,"_",l,"nt_around_TSS_TIS_TTS_TES.pdf",sep="")

# transcription start site
pdf(ref1_file,width=15,heigh=4)
layout(matrix(c(1,2,2,3,3,4), 1, 6, byrow = TRUE))
plot(ref1_down[,1],ref1_down[,2]/n,xlab = "Distance", ylab="Mean binding signals per gene",xaxt = 'n',type='l',lwd=2,col="blue",main = "",ylim=c(0,ymax))
axi_s = round(seq(0,l,length.out=5))
axis(1,at = axi_s,labels=c("TSS",axi_s[2:5]))
# translation initiation site
up_down(ref2_up,ref2_down,l,ymax,"TIS")
# translation termination site
up_down(ref3_up,ref3_down,l,ymax,"TTS")
# transcription end site
plot(-ref4_up[,1],ref4_up[,2]/n,xlab = "Distance",ylab="", yaxt='n',xaxt = 'n',type='l',lwd=2,col="blue",main = "",ylim=c(0,ymax))
axi_s = round(seq(-l,0,length.out=5))
axis(1,at = axi_s,labels=c(axi_s[1:4],"TES"))
dev.off()


# determine the GR and RS for PRAS by genomic region RPK sum and binding profile in defined region.
GR_candidates = c("5UTR","CDS","3UTR")
ref_site_candidates = c("5","3")
utr5_sum = sum(ref1_down$value)+sum(ref2_up$value)
cds_sum = sum(ref2_down$value)+sum(ref3_up$value)
utr3_sum = sum(ref3_down$value)+sum(ref4_up$value)
The_GR = GR_candidates[which.max(c(utr5_sum,cds_sum,utr3_sum))]

# Also estimate d0 using the proximal value of examined region.
estimate_d0 = function(y){
        l_en = NROW(y)
        r_eg = round(l_en/100)#examine the proximal region of 1% length of the entire checked region.
        y_reg = mean(y[(l_en-r_eg+1):l_en])
        d0_est = l_en/(-log(y_reg/max(y)))
        d0_est = 10*(round(d0_est/10))
        return(d0_est)
}


if (The_GR=="5UTR"){
    The_RS = ref_site_candidates[which.max(c(sum(ref1_down$value),sum(ref2_up$value)))]
    if (The_RS == "5"){
        y = ref1_down$value
        d0_est = estimate_d0(y)
    }else {
        y = ref2_up$value
        d0_est = estimate_d0(y)
    }
} else if (The_GR=="CDS"){
    The_RS = ref_site_candidates[which.max(c(sum(ref2_down$value),sum(ref3_up$value)))]
    if (The_RS == "5"){
        y = ref2_down$value
        d0_est = estimate_d0(y)
    }else {
        y = ref3_up$value
        d0_est = estimate_d0(y)
    }
} else if (The_GR=="3UTR"){
    The_RS = ref_site_candidates[which.max(c(sum(ref3_down$value),sum(ref4_up$value)))]
    if (The_RS == "5"){
        y = ref3_down$value
        d0_est = estimate_d0(y)
    }else {
        y = ref4_up$value
        d0_est = estimate_d0(y)
    }
}


# give suggestion for running PRAS
#cat(paste("One can use -s",The_GR,"and -r",The_RS,"for PRAS running after running this 'check' mode."))
ofilename = paste(pref,"_",l,"nt_PRAS_option_suggestions.txt",sep="")
write.table(matrix(c(The_GR,The_RS,d0_est),ncol=3),file=ofilename,row.names=F,col.names=c("-s","-r","-d"),quote=F,sep="\t")



# Command:
# Rscript --vanilla PRAS_refsite_check.r binding_tab gene_list distance plot_filename
# The script will give suggestion on the "genomic region", "reference site", and decay parameter "d0" as options for PRAS score mode running.

