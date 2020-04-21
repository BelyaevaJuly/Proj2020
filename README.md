# PARS
Tree&Distance_matrices to random protein alignment
arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        Input dir with alignment database
  -w WORKING, --working WORKING
                        Dir with dist matrixes and aligments for work
  -b, --block           Put this flag if you need to get block-alignment
  -p PERCENT, --percent PERCENT
                        Percent of block-cutting
  -t TREES, --trees TREES
                        Input file with names of files with trees for
                        alignment-generation// .nwk
  -n NUMBER, --number NUMBER
                        The number of random alignments you want to be
                        generated
  -l LENGTH, --length LENGTH
                        The length of output sequences in random alignment
  -o OUTPUT, --output OUTPUT
                        Output dir for generated random alignments
                        
HOW TO RUN (if you already have multiple protein alignments without gaps and corresponding distance matrices):

$ python random_aln_generator.v3.py -d . -w ./path-to-Alignments-database-and-Distance-matrices/ -o . -t ./path-to-tree-in-nwk-format/
 -n <number of replics> -l <length>
