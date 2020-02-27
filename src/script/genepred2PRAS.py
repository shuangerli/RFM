## gtf file downloaded from UCSC can be transferred to genePred format using the gtfToGenepred tool from UCSC.
## This script transfer the genePred format to the file format used by PRAS.
## original programmer: Jianan Lin (jianan.lin@jax.org)
## Last edit: Jan 28th, 2019
## Contact by email for issues.

import sys
import os

def get_ID_d(ifilename):
    ifile = open(ifilename)
    d = {}
    for line in ifile:
        tmp = line.strip().split()
        d[tmp[0]] = tmp[1] # transcript ID as key, gene name as value
    ifile.close()
    return d

def reformat_anno(genepredfile, d, ofilename):
    ifile = open(genepredfile)
    ofile = open(ofilename,"w")
    for line in ifile:
        tmp = line.strip().split("\t")
        estarts = map(int,tmp[8].strip(",").split(","))
        eends = map(int,tmp[9].strip(",").split(","))
        if len(estarts) == len(eends):
            elens = [x-y for x,y in zip(eends, estarts)]
            try:
                gname = d[tmp[0].split("_dup")[0]]
                l2w = [gname,tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[7],",".join(map(str,elens)),",".join(map(str,estarts))]
                ofile.write("\t".join(l2w)+"\n")
            except KeyError:
                print "Check transcript's annotation:"+tmp[0]
        else:
            print "Exon starts and Exon ends differ in length. Please check annotation."
    ofile.close()
    ifile.close()

def genepred2PRAS(id_file,genepredfile,ofilename):
    d = get_ID_d(id_file)
    reformat_anno(genepredfile,d,ofilename)

args = sys.argv
genepred2PRAS(args[1],args[2],args[3])

