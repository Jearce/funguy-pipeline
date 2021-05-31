from Bio import SeqIO

def eggnog_get_fasta(annofile, prot_file):
	with open(annofile,"r") as file:
		all_anno={}
		for line in file:
			if not line.startswith('#'):
				all_anno[line.split('\t')[0]]=line.split('\t')[7]
				current_anno=line.split('\t')[0]
	with open(prot_file,'r') as orig, open('eggnog_out','w') as output:
		records=SeqIO.parse(orig,'fasta')
		for record in records:
			#print(record.id+'\t'+record.seq)
			for k, v in all_anno.items():
				if k == record.id:
					record.id='eggnog_'+record.id+'_'+v
					record.description=''	  #'eggnog_'+record.description+'_'+v
			SeqIO.write(record, output,'fasta')


rule create_database_eggnog: #emapper.py --list_taxa > taxa_eggnog.txt to get list of taxa. 4751 is fungi
	output:
		"eggnog_database/fungi.dmnd",
		"eggnog_database/eggnog.db",
		"eggnog_database/e5.proteomes.faa"
	run:
		shell("if [ -d eggnog_database ]; then echo eggnog_database folder already exists; else mkdir eggnog_database; fi")
		shell("download_eggnog_data.py --data_dir eggnog_database/ -y")
		shell("create_dbs.py -y --taxa Fungi --data_dir eggnog_database --dbname fungi")
		
rule eggnog_anno: ##assume this fungi and use diamond to annotate
					# assume input is protein
	input:
		prot_file="{species}/predicted_prot.fa", ## need predicted gene file name here
		database="eggnog_database/fungi.dmnd",
		eggnog_data="eggnog_database/eggnog.db",
		proteome_data="eggnog_database/e5.proteomes.faa"
	output:
		"{species}/eggnog_anno/eggnog_out.fa"
	threads: 4
	run:
		shell("if [ -d {wildcards.species}/eggnog_anno ]; then echo eggnog_anno folder already exists; else mkdir {wildcards.species}/eggnog_anno; fi")
		shell("emapper.py -m diamond -i {input.prot_file} --dmnd_db {input.database} "
		"-o eggnog_{wildcards.species} --output_dir {wildcards.species}/eggnog_anno "
		"--cpu {threads} --override --data_dir eggnog_database/") 	# can add diamond sensitivity, tax_scope (more specific taxa) and 
										# what level to retrieve annotation ( eg tax scope is gammabact but anno is at bacteria level)
		eggnog_get_fasta(wildcards.species+"/eggnog_anno/eggnog_ecoli.emapper.annotations",input.prot_file)
		shell("mv eggnog_out {output}")
