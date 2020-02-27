## This script extracted the 5'UTR, CDS, and 3'UTR annotation from the transcript annotation in genepred format. Then output the files in genepred format. 
## non-coding transcripts will be automatically removed.

## original programmer: Jianan Lin (jianan.lin@jax.org)
## Last edit: Jan 28th, 2019
## Contact by email for issues.

import sys
import os

### function to get the exon index of a position ###
def get_idx(pos_ls,pos):
    pl = map(int, pos_ls)
    p = int(pos)
    k = 0
    for i,item in enumerate(pl):
        if p >= item:
            k = i
    return k

### split a given region to exons ###
def segment_transcript(start,end,exon_starts,exon_ends):
    s = int(start)
    e = int(end)
    ess = map(int,exon_starts)
    ees = map(int,exon_ends)
    i = get_idx(ess,s)
    j = get_idx(ess,e)
    if i < j:
        end_s = [ees[m] for m in range(i,j)]+[e]
        start_s = [s]+ [ess[n] for n in range((i+1),(j+1))]
        k = j - i + 1
    elif i==j:
        end_s = [e]
        start_s = [s]
        k = 1
    else:
        start = 0
        end = 0
        start_s = [0]
        end_s = [0]
        k = 0
        print "Error! Please check annotation file!"
    return start,end,start_s,end_s,k

def separate_gr(genepred):
    ifile = open(genepred)
    ofile_5utr = open(genepred+".utr5.genepred","w")
    ofile_cds = open(genepred+".cds.genepred","w")
    ofile_3utr = open(genepred+".utr3.genepred","w")
    for line in ifile:
        tmp = line.strip().split("\t")
        if tmp[5] == tmp[6]:
            continue
        elif int(tmp[5])<int(tmp[6]):
            starts = tmp[8].strip(",").split(",")
            ends = tmp[9].strip(",").split(",")
            ## upstream on chr location
            start,end,start_s,end_s,k = segment_transcript(tmp[3],tmp[5],starts,ends)
            l2w = "\t".join(map(str,[tmp[0],tmp[1],tmp[2],start,end,start,end,k,",".join(map(str,start_s))+",",",".join(map(str,end_s))+","]))
            if tmp[2] == "+":
                ofile_5utr.write(l2w+"\n")
            elif tmp[2] == "-":
                ofile_3utr.write(l2w+"\n")
            ## cds 
            start,end,start_s,end_s,k = segment_transcript(tmp[5],tmp[6],starts,ends)
            l2w = "\t".join(map(str,[tmp[0],tmp[1],tmp[2],start,end,start,end,k,",".join(map(str,start_s))+",",",".join(map(str,end_s))+","]))
            ofile_cds.write(l2w+"\n")
            ## downstream on chr location
            start,end,start_s,end_s,k = segment_transcript(tmp[6],tmp[4],starts,ends)
            l2w = "\t".join(map(str,[tmp[0],tmp[1],tmp[2],start,end,start,end,k,",".join(map(str,start_s))+",",",".join(map(str,end_s))+","]))
            if tmp[2] == "-":
                ofile_5utr.write(l2w+"\n")
            elif tmp[2] == "+":
                ofile_3utr.write(l2w+"\n")
        elif int(tmp[5])>int(tmp[6]):
            print "Error! Check annotation file!"
            continue
    ofile_5utr.close()
    ofile_cds.close()
    ofile_3utr.close()
    ifile.close()


args = sys.argv

separate_gr(args[1])




