import sys
from collections import Counter

def incpg(cpgloc, mutid):
	x = mutid.split(':')
	chr = x[0]
	pos = int(x[1])
	locs = cpgloc[chr]
	for s,e in locs:
		if pos >= s and pos <= e:
			return True
	return False


def main():
	transition = ["A>G","G>A","C>T","T>C"]
	if len(sys.argv) < 5:
		print("Usage: python stat.py input_mutations_file cpgfile famfile outfile")
		exit(1)
	infile = sys.argv[1]
	#incl_file = sys.argv[2]
	cpgfile = sys.argv[2]
	famfile = sys.argv[3]
	outfile = sys.argv[4]
	cpgloc = {}
	#incloc = [l.strip() for l in open(incl_file)]

	clines = [l.strip() for l in open(cpgfile)]
	for l in clines:
		print(l)
		if l[0] != "#":
			x = l.split('\t')
			if x[0] not in cpgloc.keys():
				cpgloc[x[0]] = [(int(x[1]), int(x[2]))]
			else:
				cpgloc[x[0]].append((int(x[1]), int(x[2])))


	famlines = [l.strip() for l in open(famfile)]
	names = []
	for l in famlines:
		if '\t' in x:
			x = l.split('\t')
		else:
			x = l.split()
		names.append(x[1])


	lines = [l.strip() for l in open(infile)]
	inmap = {}

	for n in names:
		inmap[n] = [0 for _ in range(28)]

	for l in lines[1:]:
		x = l.split('\t')
		mutid = x[1]
		type = x[2]
		chg = x[3]
		dbsnp = x[4]
		eff = x[6]
		mutinds =[v for v in  x[-1].split(';') if v in names]
		z = Counter(mutinds)
		zk = list(z.keys())
		zv = list(z.values())
		notz = list(set(names)-set(zk))
		for v in notz:
			inmap[v][3] += 1
		for k in zk:
			id = k
			if z[k] == 2:
				if type == "INDEL":
					inmap[id][26] += 1
					inmap[id][2] += 2
				elif type == "SNP":
					inmap[id][5] += 1
			elif z[k] == 1:
				if type == "INDEL":
					inmap[id][27] += 1
					inmap[id][2] += 1
				elif type == "SNP":
					inmap[id][4] += 1
					if chg in transition:
						inmap[id][0] += 1
					else:
						inmap[id][1] += 1
					if dbsnp == '.':
						inmap[id][6] += 1
						if incpg(cpgloc, mutid):
							inmap[id][11] += 1
						if chg in transition:
							inmap[id][14] += 1
							if incpg(cpgloc,mutid) == False:
								inmap[id][18] += 1
						else:
							inmap[id][15] += 1
							if incpg(cpgloc,mutid) == False:
								inmap[id][19] += 1
					else:
						inmap[id][7] += 1
						if incpg(cpgloc, mutid):
							inmap[id][10] += 1
						if chg in transition:
							inmap[id][12] += 1
							if incpg(cpgloc, mutid) == False:
								inmap[id][16] += 1
						else:
							inmap[id][13] += 1
							if incpg(cpgloc, mutid) == False:
                                                        	inmap[id][17] += 1
					if 'synonymous' in eff:
						if chg in transition:
							inmap[id][20] += 1
						else:
							inmap[id][21] += 1
					elif 'missense' in eff:
						if chg in transition:
							inmap[id][22] += 1
						else:
							inmap[id][23] += 1
					elif 'stop_gained' in eff:
						if chg in transition:
							inmap[id][24] += 1
						else:
							inmap[id][25] += 1
	of = open(outfile,"w")
	of.write("ID\tR/A_Ts\tR/A_Tv\tIndel\tRef\tHet\tHom\tnovelSNP\tknownSNP\tA/D_Ts\tA/D_Tv\tknownCpG\tnovelCpG\tknownTs\tknownTv\tnovelTs\tnovelTv\tnCpG-K_Ts\tnCpG-K_Tv\tnCpG-N_Ts\tnCpG-N_Tv\tsynonymousTs\tsynonymousTv\tmissenseTs\tmissenseTv\tnonsenseTs\tnonsenseTv\thomINDEL\thetINDEL\n")
	for i in inmap.keys():
		s = i
		for v in inmap[i]:
			s += '\t' + str(v)
		s += '\n'
		of.write(s)
	of.close()
	
if __name__=="__main__":
	main()
