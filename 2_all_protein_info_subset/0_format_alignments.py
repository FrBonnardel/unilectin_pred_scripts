domain='not yet any value'
protein='not yet any value'
score = 0
begin = 0
align_ref=''
align_match=''
outfile = open("python_output/align_clean.txt", 'w')
with open("hmmer_output/lectinpred_align_raw_refseq.txt", "r") as file:
	for line in file:
		line = line.rstrip('\n')
		if line == '':
			continue
		if 'Query:' in line or '>> ' in line or  '  == domain ' in line:
			if score != 0:
				outfile.write(domain+'\t'+protein+'\t'+score+'\t'+begin+'\t'+align_ref+'\t'+align_match+'\n')
				align_ref=''
				align_match=''
				score = 0
				begin = 0
		if 'Query:' in line:
			domain = line.split()[1]
		elif '>> ' in line:
			protein = line.split()[1]
		elif '  == domain ' in line:
			score = line.split()[4]
		elif domain in line and len(line.split()) < 5:
			align_ref += line.split()[2]
		elif protein in line and len(line.split()) < 5:
			align_match += line.split()[2]
			if begin == 0:
				begin = line.split()[1]
#Write last
if score != 0:
	outfile.write(domain+'\t'+protein+'\t'+score+'\t'+begin+'\t'+align_ref+'\t'+align_match+'\n')

domain='not yet any value'
protein='not yet any value'
score = 0
begin = 0
align_ref=''
align_match=''
with open("hmmer_output/lectinpred_align_raw_uniprot.txt", "r") as file:
	for line in file:
		line = line.rstrip('\n')
		if line == '':
			continue
		if 'Query:' in line or '>> ' in line or  '  == domain ' in line:
			if score != 0:
				outfile.write(domain+'\t'+protein+'\t'+score+'\t'+begin+'\t'+align_ref+'\t'+align_match+'\n')
				align_ref=''
				align_match=''
				score = 0
				begin = 0
		if 'Query:' in line:
			domain = line.split()[1]
		elif '>> ' in line:
			protein = line.split()[1]
			if '|' in protein:
				protein = protein.split('|')[1]
		elif '  == domain ' in line:
			score = line.split()[4]
		elif domain in line and len(line.split()) < 5:
			align_ref += line.split()[2]
		elif protein in line and len(line.split()) < 5:
			align_match += line.split()[2]
			if begin == 0:
				begin = line.split()[1]
#Write last
if score != 0:
	outfile.write(domain+'\t'+protein+'\t'+score+'\t'+begin+'\t'+align_ref+'\t'+align_match+'\n')
outfile.close()
