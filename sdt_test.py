import sys
import os
from scipy.stats import binomtest as btest

def main():
	if len(sys.argv) < 4:
        print("Usage: python3 generate_sdt_counts.py input_mutations_file famfile outputfile")
		exit(1)
		
	infile = sys.argv[1]
	famfile = sys.argv[2]
	outfile = sys.argv[3]
	
	fam = [l.strip() for l in open(famfile)]
	cases = []
	controls = []
	for l in fam:
		x = l.split()
		if x[5] == '2':
			cases.append(x[1])
		elif x[5] == '1':
			controls.append(x[1])
            
	Ncase = len(cases)
	Ncontrol = len(controls)
	Ntotal = Ncase + Ncontrol
    
	singletonP = Ncase / Ntotal
	doubletonP = 2 * Ncase/Ntotal * Ncontrol/Ntotal
	tripletonP = 1 - (Ncase/Ntotal)**3 - (Ncontrol/Ntotal)**3

	scase = [0,0]
	scon = [0,0]
	dcase = [0,0]
	dcon = [0,0]
	dsh = [0,0]
	tcase = [0,0]
	tcon = [0,0]
	tsh = [0,0]

	lines = [l.strip() for l in open(infile)]
	for l in lines[1:]:
		x = l.split('\t')
		if x[2] == "INDEL":
			continue
		y = x[23].split(';')
		ncase = len([v for v in y if v in cases])
		ncontrol = len([v for v in y if v in controls])
		csq = x[6]
		ll = ncase + ncontrol
		if ll > 3:
			continue
		elif ll == 1:
			if csq == 'synonymous_variant':
				if ncase == 1:
					scase[0] += 1
				elif ncontrol == 1:
					scon[0] += 1
			elif csq == 'missense_variant':
				if ncase == 1:
					scase[1] += 1
				elif ncontrol == 1:
					scon[1] += 1
		elif ll == 2:
			if csq == 'synonymous_variant':
				if ncase == 2:
					dcase[0] += 1
				elif ncontrol == 2:
					dcon[0] += 1
				elif ncase == 1 and ncontrol == 1:
					dsh[0] += 1
			elif csq == 'missense_variant':
				if ncase == 2:
					dcase[1] += 1
				elif ncontrol == 2:
					dcon[1] += 1
				elif ncase == 1 and ncontrol == 1:
					dsh[1] += 1
		elif ll == 3:
			if csq == 'synonymous_variant':
				if ncase == 3:
					tcase[0] += 1
				elif ncontrol == 3:
					tcon[0] += 1
				elif (ncase == 2 and ncontrol == 1) or (ncase == 1 and ncontrol == 2):
					tsh[0] += 1
			elif csq == 'missense_variant':
				if ncase == 3:
					tcase[1] += 1
				elif ncontrol == 3:
					tcon[1] += 1
				elif (ncase == 2 and ncontrol == 1) or (ncase == 1 and ncontrol == 2):
					tsh[1] += 1
					
    
    syn_singleton_pval = btest(scase[0], scase[0]+scon[0], singletonP, alternative="greater").pvalue
    syn_doubleton_pval = btest(dsh[0], dcase[0]+dcon[0]+dsh[0], doubletonP, alternative="less").pvalue
    syn_tripleton_pval = btest(tsh[0], tcase[0]+tcon[0]+tsh[0], tripletonP, alternative="less").pvalue
    
    mis_singleton_pval = btest(scase[1], scase[1]+scon[1], singletonP, alternative="greater").pvalue
    mis_doubleton_pval = btest(dsh[1], dcase[1]+dcon[1]+dsh[1], doubletonP, alternative="less").pvalue
    mis_tripleton_pval = btest(tsh[1], tcase[1]+tcon[1]+tsh[1], tripletonP, alternative="less").pvalue
    
    
    of = open(outfile, "w")
    of.write("Synonymous\n")
    of.write("Singleton\t"+str(scase[0]) + '\t' + str(scon[0]) + '\t\t' + str(singletonP) + '\t' + str(syn_singleton_pval) + '\n')
    of.write("Doubleton\t"+str(dcase[0]) + '\t' + str(dcon[0]) + '\t' + str(dsh[0]) + '\t' + str(doubletonP) + '\t' + str(syn_doubleton_pval) + '\n')
    of.write("Tripleton\t"+str(tcase[0]) + '\t' + str(tcon[0]) + '\t' + str(tsh[0]) + '\t' + str(tripletonP) + '\t' + str(syn_tripleton_pval) + '\n')
    
    of.write("Missense\n")
    of.write("Singleton\t"+str(scase[1]) + '\t' + str(scon[1]) + '\t\t' + str(singletonP) + '\t' + str(syn_singleton_pval) + '\n')
    of.write("Doubleton\t"+str(dcase[1]) + '\t' + str(dcon[1]) + '\t' + str(dsh[1]) + '\t' + str(doubletonP) + '\t' + str(syn_doubleton_pval) + '\n')
    of.write("Tripleton\t"+str(tcase[1]) + '\t' + str(tcon[1]) + '\t' + str(tsh[1]) + '\t' + str(tripletonP) + '\t' + str(syn_tripleton_pval) + '\n')
    
    of.close()
    
if __name__=="__main__":
	main()
