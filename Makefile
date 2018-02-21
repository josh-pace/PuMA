PUMA_OUT = puma-out

ky: clean
	./ky_run.py -i BPV2_REF.gb -f genbank -s e2bs --blastdb_dir=blast_database

e2bs: clean
	./puma_run.py -i BPV2_REF.gb -s e2bs --blastdb_dir=blast_database

e1e7: clean
	./puma_run.py -i BPV2_REF.gb -s e7,e1 --blastdb_dir=blast_database

clean:
	rm -rf $(PUMA_OUT)
