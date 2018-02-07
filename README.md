# PuMA

Pipeline between Genbank and PaVE

Authors: Koenraad Van Doorslaer, Josh Pace

University of Arizona

# Dependencies

Please install the following:

* Python 3.x
* BioPython
* NCBI BLAST+ 2.7.x
* MEME, FIMO (http://meme-suite.org/)

# To run:

## Make BLAST db

Run `make` in "blast_database":

    $ (cd blast_database && make)
   

Then type:

    $ make e2bs
