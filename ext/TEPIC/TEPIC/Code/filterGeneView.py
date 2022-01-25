import sys

#Removes genes for which no score good be computed

def isValidAffinity(lineSplit):
	for i in range(1,len(lineSplit)):
		if ("ENSG" in lineSplit[i]):
			return False
	for i in range(1,len(lineSplit)):
		if (float(lineSplit[i]) != 0.0):
			return True	
	return False


def check_tfs(tfs,ensg_counts,total_sample_number,gene_symbol_lengths,ensg_sym_map,cutoff,outfile_1):

	outfile_tpm_values = outfile_1.replace(".txt","_TPM_values.txt")
	output_tpm_values = open(outfile_tpm_values, "w")
	header = "ENSG\tSYMBOL\tTPM\tRPK\tGENE_LENGTH"
	header = header + "\n"
	output_tpm_values.write(header)

	per_million=total_sample_number/(1000000*1.0)
	tfs=tfs[1:]

	all_rpk=0
	for ec in ensg_counts:
		ec_counts=ensg_counts.get(ec)
		ec_lengths=gene_symbol_lengths.get(ec)

		if(ec_lengths is None or ec_counts is None or ec_lengths==0):
			continue

		ec_lengths=ec_lengths/(1000*1.0)
		current_rpk= ec_counts/(ec_lengths*1.0)
		all_rpk=all_rpk+current_rpk



	scaling_factor=all_rpk/(1000000*1.0)

	should_write=[]
	should_write.append(1)

	for tf in tfs:

		if(':' in tf):
			split_tf_2=tf.split("::")
			wanna_print=False
			for s in split_tf_2:
				tf_ensg=ensg_sym_map.get(s.upper())
				ensg_count=ensg_counts.get(tf_ensg)
				ensg_length=gene_symbol_lengths.get(tf_ensg)
				if(ensg_length is None or ensg_count is None or ensg_length==0):
					continue
				else:
					ensg_length=ensg_length/(1000*1.0)
					rpk=ensg_count/(ensg_length*1.0)
					tpm = rpk / (scaling_factor * 1.0)
					if(tpm>cutoff):
						wanna_print=True
					tpm_values_line=tf_ensg+"\t"+s+"\t" + str(tpm) + "\t" + str(rpk) +"\t"+str(ensg_length)+ "\n"
					output_tpm_values.write(tpm_values_line)
			if(wanna_print):
				should_write.append(1)
			else:
				should_write.append(0)
			continue

		split_tf=tf.split("_")
		tf_ensg=ensg_sym_map.get(split_tf[0].upper())
		ensg_count=ensg_counts.get(tf_ensg)
		ensg_length=gene_symbol_lengths.get(tf_ensg)
		if(ensg_length is None):
			should_write.append(0)
			continue
		ensg_length=ensg_length/(1000*1.0)
		rpk=ensg_count/(ensg_length*1.0)
		tpm=rpk/(scaling_factor*1.0)
		if(tpm>=cutoff):
			should_write.append(1)
		else:
			should_write.append(0)
		tpm_values_line = tf_ensg + "\t" + split_tf[0] + "\t" + str(tpm) + "\t" + str(rpk) + "\t" + str(ensg_length) + "\n"
		output_tpm_values.write(tpm_values_line)

	output_tpm_values.close()

	return should_write


def main():
	#Checking Affinity
	infile=open(sys.argv[1],"r")
	outfile_1=sys.argv[1].replace(".txt","_Filtered.txt")
	output=open(outfile_1,"w")
	#Copy header line
	header=infile.readline()
	output.write(header)
	#Check individual lines
	for l in infile:
		if (isValidAffinity(l.split())):
			output.write(l)
	infile.close()
	output.close()

	if(sys.argv[2]!="NOT_SET" and sys.argv[2]!=""):
		infile_counts=open(sys.argv[2],"r")
		infile_counts.readline()
		counts=[]
		for l in infile_counts:
			counts.append(int(l.replace("\n","")))
		infile_counts.close()

		infile_annot_counts = open(sys.argv[4],"r")
		infile_annot_counts.readline()
		ensg_numbers=[]
		for l in infile_annot_counts:
			ensg_numbers.append(l.replace("\n",""))
		infile_annot_counts.close()

		total_sample_number=0
		ensg_counts={}
		i=0
		for ensg in ensg_numbers:
			ensg_counts[ensg]=counts[i]
			total_sample_number=total_sample_number+counts[i]
			i=i+1

		infile_ensg_sym = open(sys.argv[5],"r")
		infile_ensg_sym.readline()
		ensg_sym_map={}
		for l in infile_ensg_sym:
			split = l.split()
			if(len(split)>1):
				ensg_sym_map[split[1].replace("\n","").upper()]=split[0].replace("\n","")
		infile_ensg_sym.close()
		infile_gtf=open(sys.argv[6],"r")
		gene_symbol_lengths={}
		for l in infile_gtf:
			if(l.startswith("#")):
				continue
			split = l.split()
			if(split[2]=="transcript"):
				start=int(split[3])
				end=int(split[4])
				diff=end-start+1
				ensg_name=split[9].replace("\"","")
				ensg_name=ensg_name.replace(";","")
				ensg_name=ensg_name.split(".")
				gene_symbol_lengths[ensg_name[0]]=diff
		infile_gtf.close()
		infile_tpm=open(outfile_1,"r")
		tfs=infile_tpm.readline()
		tfs_splitted=tfs.split()
		valid_columns=check_tfs(tfs_splitted,ensg_counts,total_sample_number,gene_symbol_lengths,ensg_sym_map,float(sys.argv[3]), outfile_1)
		outfile_2=outfile_1.replace(".txt","_TPM.txt")
		output_tpm = open(outfile_2, "w")
		header=""
		c=0
		for i in valid_columns:
			if(i==1):
				if(c>0):
					header+="\t"+tfs_splitted[c]
				else:
					header+=tfs_splitted[c]
			c=c+1

		header=header+"\n"
		output_tpm.write(header)
		for l in infile_tpm:
			split=l.split()
			header=""
			c = 0
			for i in valid_columns:
				if (i == 1):
					if (c > 0):
						header += "\t" + split[c]
					else:
						header += split[c]
				c = c + 1
			output_tpm.write(header+"\n")

		output_tpm.close()
		infile_tpm.close()





main()
