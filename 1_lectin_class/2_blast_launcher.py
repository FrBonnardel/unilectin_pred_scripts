import os

command = ("C:/NCBI/blast/bin/makeblastdb.exe -in 1_unilectin_seq_clean.fasta  -dbtype prot -out lectin ")
os.system(command)

command = ("C:/NCBI/blast/bin/blastp.exe -db lectin -query 1_unilectin_seq_clean.fasta -out 2_blast_results.txt -word_size 2 -outfmt 7 -evalue 1.0e-3 -num_threads 3 ")
os.system(command)
