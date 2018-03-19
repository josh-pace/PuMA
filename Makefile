PUMA_OUT = puma-out

e2bs: clean
	./puma_run.py -i BPV2_REF.gb -f genbank -s e2bs --blastdb_dir=blast_database

e1e7: clean
	./puma_run.py -i BPV2_REF.gb -f genbank -s e7,e1 --blastdb_dir=blast_database

e5: clean
	./puma_run.py -i BPV2_REF.gb -f genbank -s e5 --blastdb_dir=blast_database

e1bs: clean
	./puma_run.py -i BPV2_REF.gb -f genbank -s e1bs --blastdb_dir=blast_database

all: clean
	./puma_run.py -i BPV2_REF.gb -f genbank -s all --blastdb_dir=blast_database

l1: clean
	./puma_run.py -i BPV2_REF.gb -f genbank -s l1 --blastdb_dir=blast_database

ky: clean
        ./ky_run.py -i BPV2_REF.gb -f genbank -s e2bs --blastdb_dir=blast_database

simple:
	./run.sh -f BPV2_REF.fa

clean:
	rm -rf $(PUMA_OUT)
