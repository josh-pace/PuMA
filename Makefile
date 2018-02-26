PUMA_OUT = puma-out

e2bs: clean
	./puma_run.py -i BPV2_REF.gb -s e2bs --blastdb_dir=blast_database
e1e7: clean
        ./puma_run.py -i BPV2_REF.gb -s e7,e1 --blastdb_dir=blast_database
e5: clean
        ./puma_run.py -i BPV2_REF.gb -s e5 --blastdb_dir=blast_database
e1bs: clean
        ./puma_run.py -i BPV2_REF.gb -s e1bs --blastdb_dir=blast_database
all: clean
        ./puma_run.py -i BPV2_REF.gb -s all --blastdb_dir=blast_database

backup_e2bs: clean
	./puma_run_backup.py -i BPV2_REF.gb -s e2bs --blastdb_dir=blast_database

backup_e1e7: clean
	./puma_run_backup.py -i BPV2_REF.gb -s e7,e1 --blastdb_dir=blast_database
backup_e5: clean
	./puma_run_backup.py -i BPV2_REF.gb -s e5 --blastdb_dir=blast_database
backup_e1bs: clean
	./puma_run_backup.py -i BPV2_REF.gb -s e1bs --blastdb_dir=blast_database
backup_all: clean
	./puma_run_backup.py -i BPV2_REF.gb -s all --blastdb_dir=blast_database
ky: clean
        ./ky_run.py -i BPV2_REF.gb -f genbank -s e2bs --blastdb_dir=blast_database
jp: clean
	./learning_ky_run.py -i BPV2_REF.gb -f genbank -s e2bs --blastdb_dir=blast_database

clean:
	rm -rf $(PUMA_OUT)
