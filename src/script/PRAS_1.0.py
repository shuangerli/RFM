#!/usr/bin/env python
'''
This script is the first version of PRAS.
Original Programmer: jianan.lin@jax.org
Version: 1.0
last edit date: Feb 7th, 2019
'''
import sys
import os
import multiprocessing
from subprocess import *
import argparse
from argparse import RawTextHelpFormatter
import time
import os.path
import math
import pandas as pd
from math import log

# produce annotation file in PRAS style
def produce_anno(gtf_file,id_file):
    cmd = "./gtf2PRAS.sh {0} {1}".format(gtf_file,id_file)
    FNULL = open(os.devnull, 'w')
    Popen(cmd,shell=True, stdout=FNULL, stderr=STDOUT).communicate()


# split annotations to exons
def split_exons(annotation_file):
    ifile = open(annotation_file)
    ofile = open(annotation_file[:-4]+'_split.bed','w')
    for line in ifile:
        [ge_ne,tr_an,chr,strand,start,end,num_exon,leng_th,st_arts] = line.strip('\r\n').split('\t')
        for i in range(int(num_exon)):
            st_art=st_arts.split(',')[i]
            e_nd=int(st_art)+int(leng_th.split(',')[i])
            ofile.write(chr+'\t'+st_art+'\t'+str(e_nd)+'\t'+ge_ne+'_'+str(i)+'\t'+tr_an+'_'+str(i)+'\t'+strand+'\n')
    ofile.close()
    ifile.close()
    cmd = "cat {0} | cut -f1 | sort | uniq > {1}".format(annotation_file,"_all_genels.txt")
    Popen(cmd,shell=True,stdout=PIPE).communicate()

# split entire transcript region to introns
def split_introns(annotation_file):
    ifile = open(annotation_file)
    ofile = open(annotation_file[:-4]+'_intron_split.bed','w')
    for line in ifile:
        [ge_ne,tr_an,chr,strand,start,end,num_exon,leng_th,st_arts] = line.strip('\r\n').split('\t')
        for i in range(int(num_exon)-1):
            ex_st_art = st_arts.split(',')[i]
            intr_start = int(ex_st_art)+int(leng_th.split(',')[i])
            intr_end = st_arts.split(',')[i+1] 
            ofile.write(chr+'\t'+str(intr_start)+'\t'+str(intr_end)+'\t'+ge_ne+'_'+str(i)+'\t'+tr_an+'_'+str(i)+'\t'+strand+'\n')
    ofile.close()
    ifile.close()



# extend cross-linking site to n bp upstream and downstream
def ext_n(filename,hw):
    ifile1 = open(filename)
    hw = int(hw)
    ofile = open(filename[:-4]+"_ext_win"+str(2*hw+1)+".bed",'w')
    for line in ifile1:
        tmp=line.strip("\r\n").split("\t")
        ofile.write('\t'.join([tmp[0],str(int(tmp[1])-hw),str(int(tmp[2])+hw),".",str(abs(float(tmp[4]))),tmp[5]])+"\n")
    ifile1.close()
    ofile.close()


# assign peaks to the transcript
def assign_peaks_to_the_tran_id_mergedpeak(annotation_file,binding_file,assignment_file):
    ifile = open(annotation_file)
    d = {}
    for line in ifile:
        tmp=line.strip('\r\n').split('\t')
        d[tmp[1]]=[tmp[4],tmp[5],tmp[6],tmp[7],tmp[8],tmp[0]]
    ifile.close()
    ifile = open(binding_file)
    ofile = open(assignment_file,'w')
    ofile.write('chr\tpeak_start\tpeak_end\tcl\treads\tstrand\ttranscript\tgenename\tgstart\tgend\texon_num\texon_lens\texon_gstarts\tlocate_exon\tdis_to_utr_start\tdis_to_utr_end\tdis5\tdis3\n')
    for line in ifile:
        tmp=line.strip('\r\n').split('\t')
        tran_start,chr,hit_position_s,hit_position_e,cl,reads,strand,tran_info=int(tmp[1]),tmp[6],tmp[7],tmp[8],tmp[9],abs(float(tmp[10])),tmp[11],tmp[4]
        hit_position=str((int(hit_position_s)+int(hit_position_e))/2)#center of peak
        tmp2=tran_info.split('_')
        tran_id='_'.join(tmp2[:-1])
        transcript_to_use=d[tran_id]
        gstart=transcript_to_use[0]
        gend=transcript_to_use[1]
        gname=transcript_to_use[5]
        tran_exon_num=transcript_to_use[2]
        tran_exon_lens=map(int,transcript_to_use[3].split(','))
        tran_exon_starts=map(int,transcript_to_use[4].split(','))
        exon_index=str(tran_exon_starts.index(tran_start))
        dis_to_utr_start=max(0,(sum(tran_exon_lens[:int(exon_index)])+int(hit_position)-int(tran_exon_starts[int(exon_index)])))
        dis_to_utr_end=max(0,(sum(tran_exon_lens[int(exon_index):])-(int(hit_position)-int(tran_exon_starts[int(exon_index)]))))
        if strand=='+':
            dis5=dis_to_utr_start
            dis3=dis_to_utr_end
        elif strand=='-':
            dis3=dis_to_utr_start
            dis5=dis_to_utr_end
        ofile.write(chr+'\t'+hit_position_s+'\t'+hit_position_e+'\t'+cl+'\t'+str(reads)+'\t'+strand+'\t'+tran_id+'\t'+gname+'\t'+str(gstart)+'\t'+str(gend)+'\t'+str(tran_exon_num)+'\t'+','.join(map(str,tran_exon_lens))+'\t'+','.join(map(str,tran_exon_starts))+'\t'+exon_index+'\t'+str(dis_to_utr_start)+'\t'+str(dis_to_utr_end)+'\t'+str(dis5)+'\t'+str(dis3)+'\n')
    ofile.close()
    ifile.close()

# assign each peak to the nearest annotation
def assign_nearest_annotation(assignment_file,filter_file):
    ifile1=open(assignment_file)
    d={}
    next(ifile1)
    for line in ifile1:
        tmp=line.strip('\r\n').split('\t')
        chr,peak_start,peak_end,cl,reads,strand,transcript,genename,dis5,dis3=tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[16],tmp[17]
        try:
                info_tran=d[chr+peak_start+peak_end+reads+strand]
                if int(dis5)< int(info_tran[8]):
                        d[chr+peak_start+peak_end+reads+strand][6:9]=[transcript,genename,dis5]
                if int(dis3)< int(info_tran[11]):
                        d[chr+peak_start+peak_end+reads+strand][9:12]=[transcript,genename,dis3]
        except KeyError:
                tran5=transcript
                gene5=genename
                tran3=transcript
                gene3=genename
                d[chr+peak_start+peak_end+reads+strand]=[chr,peak_start,peak_end,cl,reads,strand,tran5,gene5,dis5,tran3,gene3,dis3]
    ifile1.close()
    ofile=open(filter_file,'w')
    ofile.write("chr\tpeak_start\tpeak_end\tcl\treads\tstrand\tend5_nearest_transcript\tend5_nearest_gene\tdis5\tend3_nearest_transcript\tend3_nearest_gene\tdis3\n")
    for key in d:
        ofile.write('\t'.join(d[key])+"\n")
    ofile.close()

# assign each peak to the nearest splice site.
# We calculate exon and intron separately.
def assign_nearest_splice(binding_file,filter_file):
    ifile = open(binding_file)
    d_peak = {}
    for line in ifile:
        tmp = line.strip().split("\t")
        chrm,peak_start,peak_end,cl,reads,strand,anno_start,anno_end,anno_gene,anno_trans = tmp[6],tmp[7],tmp[8],tmp[9],tmp[10],tmp[11],tmp[1],tmp[2],"_".join(tmp[3].split("_")[:-1]),"_".join(tmp[4].split("_")[:-1])
        peak_center = (int(peak_start)+int(peak_end))/2
        dis_left = max(0,peak_center - int(anno_start))
        dis_right = max(0,int(anno_end) - peak_center)
        peak_id = "__".join([chrm,peak_start,peak_end,cl,reads,strand])
        try:
            [old_dis_left,old_dis_right] = d_peak[peak_id][:2]
            if dis_left < old_dis_left:
                d_peak[peak_id][0]=dis_left
                d_peak[peak_id][2]=anno_gene
                d_peak[peak_id][3]=anno_trans
            if dis_right < old_dis_right:
                d_peak[peak_id][1]=dis_right
                d_peak[peak_id][4]=anno_gene
                d_peak[peak_id][5]=anno_trans
        except KeyError:
            d_peak[peak_id] = [dis_left,dis_right,anno_gene,anno_trans,anno_gene,anno_trans]
    ifile.close()
    ofile=open(filter_file,'w')
    ofile.write("chr\tpeak_start\tpeak_end\tcl\treads\tstrand\tend5_nearest_transcript\tend5_nearest_gene\tdis5\tend3_nearest_transcript\tend3_nearest_gene\tdis3\n")
    for key in d_peak:
        [dis_left,dis_right,anno_gene_l,anno_trans_l,anno_gene_r,anno_trans_r]=d_peak[key]
        peak_id = key.split("__")
        if peak_id[5]=="+":
            info_ls = peak_id+[anno_trans_l,anno_gene_l,dis_left,anno_trans_r,anno_gene_r,dis_right]
        elif peak_id[5]=="-":
            info_ls = peak_id+[anno_trans_r,anno_gene_r,dis_right,anno_trans_l,anno_gene_l,dis_left]
        ofile.write('\t'.join(map(str,info_ls))+"\n")
    ofile.close()








# check the binding profile around the candidate reference site.
def check_prof(assignment_file_ls,interested_genels,dis,plot_filename):
    cmd = "Rscript --vanilla PRAS_refsite_check.r {0} {1} {2} {3} {4} {5}".format(assignment_file_ls[0],assignment_file_ls[1],assignment_file_ls[2],interested_genels,dis,plot_filename)
    Popen(cmd,shell=True,stdout=PIPE).communicate()

# check the binding profile around the splicing site
def check_splice_prof(assignment_file_ls,interested_genels,dis,plot_filename):
    cmd = "Rscript --vanilla PRAS_splice_site_check.r {0} {1} {2} {3} {4}".format(assignment_file_ls[0],assignment_file_ls[1],interested_genels,dis,plot_filename)
    Popen(cmd,shell=True,stdout=PIPE).communicate()


# This is the main function when the reference site is not given and we only consider the transcriptomic reads
def run_order_check(gtf_file,id_file,binding_file,hw,filter_file,interested_genels,refsite,d0,PRAS_score_tab,check_dis,check_filename):
    print "Start PRAS 'check' mode."
    # choose annotation of genomic region
    annotation_file_5utr = gtf_file+".utr5.PRASanno"
    annotation_file_cds = gtf_file+".cds.PRASanno"
    annotation_file_3utr = gtf_file+".utr3.PRASanno"
    anno_ls = [annotation_file_5utr,annotation_file_cds,annotation_file_3utr]
    # generate annotation in PRAS style if not exists.
    if not (os.path.isfile(annotation_file_5utr) and os.path.isfile(annotation_file_cds) and os.path.isfile(annotation_file_3utr)):
        print "Generating PRAS annotations..."
        produce_anno(gtf_file,id_file)
        print "Finished annotation generation."
    # Annotate the peaks
    print "Annotating the binding peaks..."
    filter_file_ls=[]
    for file in anno_ls:
        gr = file.split(".")[-2]
        split_exons(file)
        annotation_file1 = file[:-4]+"_split.bed"
        binding_file1=binding_file[:-4]+"_"+gr+".bed"
        cmd = "bedtools intersect -a {0} -b {1} -u -s > {2}".format(binding_file,annotation_file1,binding_file1)
        Popen(cmd,shell=True,stdout=PIPE).communicate()
        # Extend the binding region to hw upstream and downstream
        ext_n(binding_file1,hw)
        ext_binding_file = binding_file1[:-4]+"_ext_win"+str(2*hw+1)+".bed"
        intersect_file = annotation_file1[:-4]+"_intersect_win"+str(2*hw+1)+".txt"
        cmd="bedtools intersect -a {0} -b {1} -wa -wb -s > {2}".format(annotation_file1,ext_binding_file,intersect_file)
        Popen(cmd,shell=True,stdout=PIPE).communicate()
        # Assign peaks to transcripts
        assignment_file = ext_binding_file[:-4]+"_transcriptome.txt"
        assign_peaks_to_the_tran_id_mergedpeak(file,intersect_file,assignment_file)
        filter_file = binding_file[:-4] +"_"+ gr+".assign.txt"
        assign_nearest_annotation(assignment_file,filter_file)
        filter_file_ls.append(filter_file)
        print "Finished binding peak annotation in "+gr
        # remove the tmp files
        cmd = "rm "+"\nrm ".join([assignment_file,intersect_file,binding_file1,annotation_file1,ext_binding_file])+"\n"
        Popen(cmd,shell=True,stdout=PIPE).communicate()
    # reference site check
    # plot out the binding profile around the reference site candidates
    print "Generating the binding plot..."
    check_prof(filter_file_ls,interested_genels,check_dis,check_filename)
    print "Plot done."
    ifilename = check_filename + "_" + str(check_dis) + "nt_PRAS_option_suggestions.txt"
    ifile = open(ifilename)
    next(ifile)
    for line in ifile:
        tmp = line.strip().split("\t")
    print "The {0}' end of {1} is the most enriched of CLIP-seq signals. One can use -s {2}, -r {3} and -d {4} for PRAS 'score' mode after running 'check' mode.".format(tmp[1],tmp[0],tmp[0],tmp[1],tmp[2])
    ifile.close()     

# This is the main function when the reference site is not given and we consider all pre-mRNA level reads.
def run_order_check_splice(gtf_file,id_file,binding_file,hw,filter_file,interested_genels,refsite,d0,PRAS_score_tab,check_dis,check_filename):
    print "Start PRAS 'check' mode for splicing site."
    # choose annotation of genomic region
    annotation_file_transcript = gtf_file+".transcript.PRASanno"
    anno_file = annotation_file_transcript
    # generate annotation in PRAS style if not exists.
    if not (os.path.isfile(annotation_file_transcript)):
        print "Generating PRAS annotations..."
        produce_anno(gtf_file,id_file)
        print "Finished annotation generation."
    # Annotate the peaks
    print "Annotating the binding peaks..."
    split_exons(anno_file)
    split_introns(anno_file)
    annotation_file1 = anno_file[:-4]+"_split.bed"
    annotation_file2 = anno_file[:-4]+"_intron_split.bed"
    binding_file1=binding_file[:-4]+"_exon.bed"
    binding_file2=binding_file[:-4]+"_intron.bed"

    cmd = "bedtools intersect -a {0} -b {1} -u -s > {2}".format(binding_file,annotation_file1,binding_file1)
    Popen(cmd,shell=True,stdout=PIPE).communicate()
    cmd = "bedtools intersect -a {0} -b {1} -u -s > {2}".format(binding_file,annotation_file2,binding_file2)
    Popen(cmd,shell=True,stdout=PIPE).communicate()
    # Extend the EXON binding region to hw upstream and downstream
    ext_n(binding_file1,hw)
    ext_binding_file1 = binding_file1[:-4]+"_ext_win"+str(2*hw+1)+".bed"
    intersect_file1 = annotation_file1[:-4]+"_intersect_win"+str(2*hw+1)+".txt"
    cmd="bedtools intersect -a {0} -b {1} -wa -wb -s > {2}".format(annotation_file1,ext_binding_file1,intersect_file1)
    Popen(cmd,shell=True,stdout=PIPE).communicate()
    # Extend the INTRON binding region to hw upstream and downstream
    ext_n(binding_file2,hw)
    ext_binding_file2 = binding_file2[:-4]+"_ext_win"+str(2*hw+1)+".bed"
    intersect_file2 = annotation_file2[:-4]+"_intersect_win"+str(2*hw+1)+".txt"
    cmd="bedtools intersect -a {0} -b {1} -wa -wb -s > {2}".format(annotation_file2,ext_binding_file2,intersect_file2)
    Popen(cmd,shell=True,stdout=PIPE).communicate()

    # Assign peaks to splice sites
    filter_file1 = binding_file1[:-4] + ".assign.txt"
    assign_nearest_splice(intersect_file1,filter_file1)
    filter_file2 = binding_file2[:-4] + ".assign.txt"
    assign_nearest_splice(intersect_file2,filter_file2)
    filter_file_ls = [filter_file1,filter_file2]
    print "Finished binding peak annotation in Exons and Introns"
    # remove the tmp files
    cmd = "rm "+"\nrm ".join([intersect_file1,intersect_file2,binding_file1,binding_file2,annotation_file1,annotation_file2,ext_binding_file1,ext_binding_file2])+"\n"
    Popen(cmd,shell=True,stdout=PIPE).communicate()
    # reference site check
    # plot out the binding profile around the 5' and 3' splice site 
    print "Generating the binding plot..."
    check_splice_prof(filter_file_ls,interested_genels,check_dis,check_filename)
    print "Plot done."



# This is the main function when the reference site is given
# Single or multiple known reference sites can be provided by gr and refsite arguments. comma separated.
def run_order_known_ref(gtf_file,id_file,gr,binding_file,hw,filter_file,interested_genels,refsite,d0,PRAS_score_tab):
    print "Start PRAS 'score' mode."
    gr_ls = gr.split(",")
    refsite_ls = refsite.split(",")
    if len(gr_ls)!=len(refsite_ls):
        print "Please check your input -s and -r options, they should be the same length separated by comma."
        sys.exit()
    if 'splice' in gr_ls:
        print "'splice' cannot appear in multiple genomic region input."
        sys.exit()
    score_file_ls = []
    for gr,refsite in zip(gr_ls,refsite_ls):
        # choose annotation of genomic region
        if gr == "5UTR":
            annotation_file = gtf_file+".utr5.PRASanno"
        elif gr == "CDS":
            annotation_file = gtf_file+".cds.PRASanno"
        elif gr == "3UTR":
            annotation_file = gtf_file+".utr3.PRASanno"
        elif gr == "transcript":
            annotation_file = gtf_file+".transcript.PRASanno"
        # generate annotation in PRAS style if not exists.
        if not os.path.isfile(annotation_file):
            print "Generating PRAS annotations..."
            produce_anno(gtf_file,id_file)
            print "Finished annotation generation."
        # Annotate the peaks
        print "Annotating the binding peaks..."
        split_exons(annotation_file)
        annotation_file1 = annotation_file[:-4]+"_split.bed"
        binding_file1=binding_file[:-4]+"_"+gr+".bed"
        cmd = "bedtools intersect -a {0} -b {1} -u -s > {2}".format(binding_file,annotation_file1,binding_file1)
        Popen(cmd,shell=True,stdout=PIPE).communicate()
        # Extend the binding region to hw upstream and downstream
        ext_n(binding_file1,hw)
        ext_binding_file = binding_file1[:-4]+"_ext_win"+str(2*hw+1)+".bed"
        intersect_file = annotation_file1[:-4]+"_intersect_win"+str(2*hw+1)+".txt"
        cmd="bedtools intersect -a {0} -b {1} -wa -wb -s > {2}".format(annotation_file1,ext_binding_file,intersect_file)
        Popen(cmd,shell=True,stdout=PIPE).communicate()
        # Assign peaks to transcripts
        assignment_file = ext_binding_file[:-4]+"_transcriptome.txt"
        assign_peaks_to_the_tran_id_mergedpeak(annotation_file,intersect_file,assignment_file)
        assign_nearest_annotation(assignment_file,filter_file)
        print "Finished binding peak annotation."
        # calculate the PRAS scores
        print "Calculating PRAS scores..."
        cmd = "Rscript --vanilla PRAS_score.r {0} {1} {2} {3} {4}".format(filter_file,interested_genels,refsite,d0,PRAS_score_tab+"."+gr+"."+refsite+".tmp")
        Popen(cmd,shell=True,stdout=PIPE).communicate()
        print "Finished PRAS score calculation for "+gr+" from the "+refsite+"' site."
        # remove the tmp files
        cmd = "rm "+"\nrm ".join([assignment_file,intersect_file,binding_file1,annotation_file1,ext_binding_file])+"\n"
        Popen(cmd,shell=True,stdout=PIPE).communicate()
        score_file_ls.append(PRAS_score_tab+"."+gr+"."+refsite+".tmp")
    # sum all score file together for multiple refsite
    if len(gr_ls)==1:
        cmd = "mv {0} {1}".format(PRAS_score_tab+"."+gr+"."+refsite+".tmp",PRAS_score_tab)
        Popen(cmd,shell=True,stdout=PIPE).communicate()
    else:
        d_score = {}
        for score_file in score_file_ls:
            ifile = open(score_file)
            next(ifile)
            for line in ifile:
                tmp = line.strip().split("\t")
                try:
                    d_score[tmp[0]]+=float(tmp[1])
                except KeyError:
                    d_score[tmp[0]]=float(tmp[1])
            ifile.close()
            cmd = "rm "+ score_file
            Popen(cmd,shell=True,stdout=PIPE).communicate()
        # dictionary to dataframe
        df = pd.DataFrame.from_dict(d_score,orient="index")
        df.columns = ["PRAS"]
        df["log2PRAS"]=df["PRAS"].apply(lambda x: log(x+1,2))
        df["rank"]=df["log2PRAS"].rank(ascending=False)
        df.sort_values("rank", inplace = True)
        df.to_csv(PRAS_score_tab,sep="\t",index_label="Gene")
        for gr,refsite in zip(gr_ls,refsite_ls):
            score_file = PRAS_score_tab+"."+gr+"."+refsite+".tmp"
            cmd = "rm "+score_file+"\n"
            Popen(cmd,shell=True,stdout=PIPE).communicate()

def run_order_known_splice(gtf_file,id_file,gr,binding_file,hw,filter_file,interested_genels,refsite,d0,PRAS_score_tab):
    print "Start PRAS 'score' mode on splice site."
    score_file_ls = []
    # choose annotation of genomic region
    anno_file = gtf_file+".transcript.PRASanno"
    # generate annotation in PRAS style if not exists.
    if not os.path.isfile(anno_file):
        print "Generating PRAS annotations..."
        produce_anno(gtf_file,id_file)
        print "Finished annotation generation."
    # Annotate the peaks
    print "Annotating the binding peaks..."
    split_exons(anno_file)
    split_introns(anno_file)
    annotation_file1 = anno_file[:-4]+"_split.bed"
    annotation_file2 = anno_file[:-4]+"_intron_split.bed"
    binding_file1=binding_file[:-4]+"_exon.bed"
    binding_file2=binding_file[:-4]+"_intron.bed"

    cmd = "bedtools intersect -a {0} -b {1} -u -s > {2}".format(binding_file,annotation_file1,binding_file1)
    Popen(cmd,shell=True,stdout=PIPE).communicate()
    cmd = "bedtools intersect -a {0} -b {1} -u -s > {2}".format(binding_file,annotation_file2,binding_file2)
    Popen(cmd,shell=True,stdout=PIPE).communicate()
    # Extend the EXON binding region to hw upstream and downstream
    ext_n(binding_file1,hw)
    ext_binding_file1 = binding_file1[:-4]+"_ext_win"+str(2*hw+1)+".bed"
    intersect_file1 = annotation_file1[:-4]+"_intersect_win"+str(2*hw+1)+".txt"
    cmd="bedtools intersect -a {0} -b {1} -wa -wb -s > {2}".format(annotation_file1,ext_binding_file1,intersect_file1)
    Popen(cmd,shell=True,stdout=PIPE).communicate()
    # Extend the INTRON binding region to hw upstream and downstream
    ext_n(binding_file2,hw)
    ext_binding_file2 = binding_file2[:-4]+"_ext_win"+str(2*hw+1)+".bed"
    intersect_file2 = annotation_file2[:-4]+"_intersect_win"+str(2*hw+1)+".txt"
    cmd="bedtools intersect -a {0} -b {1} -wa -wb -s > {2}".format(annotation_file2,ext_binding_file2,intersect_file2)
    Popen(cmd,shell=True,stdout=PIPE).communicate()

    # Assign peaks to splice sites
    filter_file1 = binding_file1[:-4] + ".assign.txt"
    assign_nearest_splice(intersect_file1,filter_file1)
    filter_file2 = binding_file2[:-4] + ".assign.txt"
    assign_nearest_splice(intersect_file2,filter_file2)
    print "Finished binding peak annotation in Exons and Introns"
    # remove the tmp files
    cmd = "rm "+"\nrm ".join([intersect_file1,intersect_file2,binding_file1,binding_file2,annotation_file1,annotation_file2,ext_binding_file1,ext_binding_file2])+"\n"
    Popen(cmd,shell=True,stdout=PIPE).communicate()
    # calculate the PRAS scores
    print "Calculating PRAS splice scores..."
    cmd = "Rscript --vanilla PRAS_splice_score.r {0} {1} {2} {3} {4}".format(filter_file1,filter_file2,interested_genels,d0,PRAS_score_tab)
    Popen(cmd,shell=True,stdout=PIPE).communicate()
    print "Finished PRAS score calculation for splicing sites."



def garparser():
    argparser = argparse.ArgumentParser(description="PRAS: Protein-RNA Association Strength: predicting functional targets of RNA-binding proteins based on CLIP-seq peaks\nVersion: 1.0", formatter_class=RawTextHelpFormatter)
    argparser.add_argument("-g","--gtf",dest = "gtf_file", type = str, required = True, help = "Input annotation file in gtf format. The standard gtf format can be downloaded from UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables).")
    argparser.add_argument("-t","--id",dest = "id_file",type = str, required = True, help = "Input the transcript and gene ID file. This file has to have two columns, the first column is the transcript ID, and the second column is the gene id or gene name. Please check the instruction page for file example.")
    argparser.add_argument("-i","--input",dest = "binding_file", type = str,required = True, help = "Input binding file in bed format (6 columns). Please refer to the instruction webpage for detailed example.")
    argparser.add_argument("-a","--assignment",dest="assign_file", type = str,required=True, help = "The filename of one of the output files, peak annotation file. This file output the peaks with their annotations.")
    argparser.add_argument("-m","--mode",dest = "PRAS_mode", type = str, required = False, default = "check",help= "Input the PRAS mode 'check' or 'score' to run. When the user don't have pre-knowledge of the RBP's binding preference, 'check' mode will check the binding profiles around selected reference site. When the user has the pre-knowledge or has got the suggestion from 'check' mode, 'score' mode can be used for PRAS score calculation. In the 'check' mode, PRAS will plot out the binding profiles around both 5' and 3' end of the reference sites (or SS) instead of calculating the PRAS score. If 'check' mode is running without enabling splice site check (via -s splice), PRAS will give a suggestion for selection of -s, -r, and -d for running the 'score' mode afterwards. However, when '-s splice' is enabled, PRAS will not give suggestion for parameter selection. In the 'score' mode, PRAS score will be calculated directly based on parameters -s and -r.")
    argparser.add_argument("-s","--segment",dest = "genomic_region",type = str, required = False, default="transcript", help = "Input the genomic region(s). Comma separated multiple genomic regions are acceptable. For single genomic region, there are five options: '5UTR' is the 5'UTR, 'CDS' is the coding region, '3UTR' is the 3'UTR, 'transcript' is the entire transcript, and 'splice' is for splicing site (SS). If no input provided, PRAS will take the transcript as the default input. As for multiple genomic regions, any combination of '5UTR', 'CDS', and '3UTR' is accepted. 'splice' cannot appear in the multiple genomice region input.")
    argparser.add_argument("-w","--halfwindow",dest="hw",type = int,required=False,default=10,help = "Number of nt to extend towards upstream and downstream of the binding sites. If the binding file (given via -i) includes the peak regions, '-w k' will extend the peak region by k nt from both boudaries of the peak region. If the binding file (given via -i) includes the significant cross-linking sites, '-w k' will extract a region of 2*k+1 nt centered at each cross-linking site. We suggest to use 0 for this option if your binding file is peak region instead of the cross-linking sites.")
    argparser.add_argument("-l","--genelist",dest="interested_gl",type = str, required=False,default="_all_genels.txt",help = "Input an interested gene list file (1 column file without header). If no input, PRAS will take all the genes in the annotation file as the gene list.")
    argparser.add_argument("-r","--refsite",dest="reference_site",type=str,required=False,default="5",help = "Input comma separated list of 5 or 3 to tell whether the reference site is on the 5' end or the 3' end of the selected annotation regions. The input comma separated list should match the length of -s option. Only used in the 'score' mode. Users don't have to set this parameter if splice site (-s splice) is enabled, because both 5' and 3' SS score will be calculated. Single reference site can be input directly. For example, -r 3 together with -s 3UTR means calculating PRAS score from the 3' end of the 3'UTR. Multiple reference sites are acceptable via comma separation. For example, -r 3,5 together with -s 5UTR,3UTR means calculating combined PRAS score from the 3' end of the 5'UTR and 5' end of the 3'UTR.")
    argparser.add_argument("-d","--d0",dest="distance_parameter",type=int, required=False,default=1000,help = "The parameter d0 in the PRAS score calculation. This is the decay rate in the PRAS score formula. A smaller d0 stands for a faster decay. 1000nt is just an experiential selection. It is better to select d0 based on the CLIP-seq read distribution.")
    argparser.add_argument("-o","--output",dest="PRAS_score_tab",type=str,required=False,default="binding",help = "prefix of the filename of the PRAS score output table. If not provided, PRAS will generate a file named 'binding_PRAS.txt'.")
    argparser.add_argument("-c","--check_dis",dest="check_dis",type=int,required=False,default=1000,help = "If 'check' mode is enabled, this is the distance from the candidate reference site to be examined. The default is 1000nt, which means PRAS will plot out the binding profiles at the distance of 0nt to 1000nt from both the 5' end and 3' end of the selected genomic region.")
    argparser.add_argument("-p","--check_filename",dest="check_filename_prefix",type=str,required=False,default="check_profile",help = "If 'check' mode is enabled, this is the output pdf name of the plot. This plot shows the binding profiles at the distance range (defined by -c) from the reference sites. If splice sites (SS) check is enabled via -s splice, PRAS will plot the binding profile around both the 5' SS and the 3' SS. If other genomic region is given via -s, PRAS will plot the binding profile around the transcription start site, translation initiation site, translation termination site, and transcription stop site.")
    return(argparser)


def mainfun():
    argsp = garparser()
    args = argsp.parse_args()
    mode = args.PRAS_mode
    anno = args.gtf_file
    idf = args.id_file
    gr = args.genomic_region
    bind = args.binding_file
    hw = args.hw
    assignfile = args.assign_file
    gene_ls = args.interested_gl
    ref_site = args.reference_site
    d0 = args.distance_parameter
    ofile = args.PRAS_score_tab+"_PRAS.txt"
    check_dis=args.check_dis
    check_filename=args.check_filename_prefix
    start_time = time.time()
    print "Starting PRAS program ..."
    if mode == "check":
        if gr == "splice":
            run_order_check_splice(anno,idf,bind,hw,assignfile,gene_ls,ref_site,d0,ofile,check_dis,check_filename)
        else:
            run_order_check(anno,idf,bind,hw,assignfile,gene_ls,ref_site,d0,ofile,check_dis,check_filename)
    elif mode == "score":
        if gr == "splice":
            run_order_known_splice(anno,idf,gr,bind,hw,assignfile,gene_ls,ref_site,d0,ofile)
        else:
            run_order_known_ref(anno,idf,gr,bind,hw,assignfile,gene_ls,ref_site,d0,ofile)
    else:
        print "Please provide the PRAS mode to run. 'check' for examine binding profile or 'score' for PRAS score calculation."
    print "Program Finished!"
    print("--- The running time is %.2f seconds ---" % (time.time() - start_time))

mainfun()






