#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This script is used for calculating the binding signal accumulation around the splicing sites. PRAS will plot out the binding profiles around the 5' splicing and 3' splicing sites and let user decide whether to use splicing site as the reference site.


# Four inputs:
# 1. Annotated binding tables in the order of exons and introns.
# 2. Interested gene list
# 3. Range of distance for examine the profiles
# 4. output filename.

# read binding table
B_exon = read.delim(args[1],header=T)
B_intron = read.delim(args[2],header=T)

# read interested gene list
gene_ls = read.delim(args[3],header = F)
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


l = as.numeric(args[4])
# collect binding profiles
# to note, the 5' or 3' splice site are defined by the intron's orientation, so 5' splice site exon reads should be the reads on the exons upstream to this splice site, which is the 3' end of the upstream exon
splice_5_exon = get_dis(B_exon,gene_ls,3,l)
splice_5_intron = get_dis(B_intron,gene_ls,5,l)
splice_3_exon = get_dis(B_exon,gene_ls,5,l)
splice_3_intron = get_dis(B_intron,gene_ls,3,l)
# get the largest value of the intensities within the defined region
ymax = max(c(max(splice_5_exon$value),max(splice_5_intron$value),max(splice_3_exon$value),max(splice_3_intron$value)))/n


# plot the upstream and downstream profile
up_down = function(up,down,l,ymax,lb){
    plot(-up[,1],up[,2]/n,xlab = "Distance",ylab="Mean binding signals per gene",xaxt = 'n',type='l',lwd=2,col="blue",main = "",xlim=c(-l,l),ylim=c(0,ymax))
    lines(down[,1],down[,2]/n,col='blue',lwd=2)
    axi_s = round(seq(-l,l,length.out=5))
    axis(1,at = axi_s,labels=c(axi_s[1:2],lb,axi_s[4:5]))
}

# plot out the profiles
pref = args[5]
ref1_file = paste(pref,"_splice.pdf",sep="")
# 5' splice site
pdf(ref1_file,width=15,heigh=4)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
up_down(splice_5_exon,splice_5_intron,l,ymax,"5' SS")
# 3' splice site
up_down(splice_3_intron,splice_3_exon,l,ymax,"3' SS")
dev.off()




# Command:
# Rscript --vanilla PRAS_splice_site_check.r binding_tab gene_list distance plot_filename

