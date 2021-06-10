from Bio import SeqIO

configfile: "../config/annotation/config.yaml"

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

def find_eggnog_path(filename): #look up a filename in eggnog db in cluster location defined in config, if not found return local path to be build
	if os.path.exists(config["eggnog_db_path"]["cluster_path"]+'{}'.format(filename)):
		return config["eggnog_db_path"]["cluster_path"]+'{}'.format(filename)
	else:
		return config["eggnog_db_path"]["local_path"]+'{}'.format(filename)

rule create_database_eggnog: #emapper.py --list_taxa > taxa_eggnog.txt to get list of taxa
	output:
		"eggnog_database/fungi.mmseqs/fungi.mmseqs",
		"eggnog_database/eggnog.db",
		"eggnog_database/e5.proteomes.faa"
	run:
		shell("if [ -d eggnog_database ]; then echo eggnog_database folder already exists; else mkdir eggnog_database; fi")
		shell("download_eggnog_data.py -M -D --data_dir eggnog_database/ -y")
		shell("create_dbs.py -y -m mmseqs --taxa Fungi --data_dir eggnog_database --dbname fungi")
		
rule eggnog_anno: ##assume this fungi and use mmseqs to annotate
					# assume input is protein
					#if cluster database is specified in config, scripts will find relevant inputs. If inputs not there, will create local db
	input:
		prot_file="{species}/augustus.hints.aa", ## need predicted gene file name here
		database=find_eggnog_path('fungi.mmseqs/fungi.mmseqs'),
		eggnog_data=find_eggnog_path('eggnog.db'),
		proteome_data=find_eggnog_path('e5.proteomes.faa')
	output:
		"{species}/eggnog_anno/eggnog_out.fa"
	threads: 20
	run:
		shell("if [ -d {wildcards.species}/eggnog_anno ]; then echo eggnog_anno folder already exists; else mkdir {wildcards.species}/eggnog_anno; fi")
		shell("emapper.py -m mmseqs -i {input.prot_file} --mmseqs_db {input.database} "
		"-o eggnog_{wildcards.species} --output_dir {wildcards.species}/eggnog_anno "
		"--cpu {threads} --override --dbmem --data_dir eggnog_database/")
										# what level to retrieve annotation ( eg tax scope is gammabact but anno is at bacteria level)
		eggnog_get_fasta(wildcards.species+"/eggnog_anno/eggnog_ecoli.emapper.annotations",input.prot_file)
		shell("mv eggnog_out {output}")
