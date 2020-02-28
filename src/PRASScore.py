import sys

if __name__ == "__main__":
	bedFile = "intersect"
	inFile = open("OutData_500nt_PRAS_option_suggestions.txt", 'r')
	A = inFile.read().splitlines()[1].split("	")
	outFile = open("ScoreCommand.txt", 'w')
	outFile.write("python PRAS_1.0.py -g annotation.gtf -t IDs.txt -m score -s %s -i %s.bed -a %s_%s.assign.txt -w 10 -r %s -d %s" % (A[0], bedFile, bedFile, A[0], A[1], A[2]))
	outFile.close()
