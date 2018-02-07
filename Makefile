PUMA_OUT = puma-out

test: clean
	./puma_run.py -i BPV2_REF.gb -opt e7 e1 --blastdb_dir=blast_database

clean:
	rm -rf $(PUMA_OUT)
