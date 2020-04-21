'''
Tree to random alignment

To generate random alignment in your working directory should be:
a) a directory with an alignment bank
b) a directory with the distance matrices of these alignments
An output directory with random alignments will be created automatically
Example:
> python random_aln_generator.v2.py -d ./database -w ./blocks -o . -b -p 0.65 -t ./trees.nwk -n 1 -l 100     
> python random_aln_generator.v2.py -d Sample2/ -w blockalignments/ -o version2/ -t trees.txt -n 1 -l 50  
'''
import argparse
from ete3 import Tree
import os
import glob
import numpy as np
from Bio import AlignIO 
import random
import copy
from copy import deepcopy

#from nd.v2 import nd, T

parser = argparse.ArgumentParser(description='Tree to random alignment')
parser.add_argument('-d', '--database', type=str, help='Input dir with alignment database')
parser.add_argument('-w', '--working', type=str, help='Dir with dist matrixes and aligments for work')
parser.add_argument('-b', '--block', action="store_true", help='Put this flag if you need to get block-alignment')
parser.add_argument('-p', '--percent', type=float, help='Percent of block-cutting')
parser.add_argument('-t', '--trees', type=str, help='Input file with names of files with trees for alignment-generation// .nwk')
parser.add_argument('-n', '--number', type=int, help='The number of random alignments you want to be generated')
parser.add_argument('-l', '--length', type=int, help='The length of output sequences in random alignment')
parser.add_argument('-o', '--output', type=str, help='Output dir for generated random alignments')
args = parser.parse_args()

tree_files = []
with open(args.trees, 'r') as inp:
	for string in inp:
		tree_files.append(string.strip())

current_dir = os.getcwd().replace("\\", '/')
print('Current dir:	', current_dir)
database = current_dir + "/" + args.database
workdir = current_dir + "/" + args.working
outdir = current_dir + "/" + args.output

if args.block:
	print('Getting block-alignment..')
	os.chdir(database)
	for file in glob.glob("*.fasta"):
    		matrix = np.array([np.array(rec.seq) for rec in AlignIO.read(open(file, 'r'), "fasta")])
    		names = [rec.name for rec in AlignIO.read(open(file, 'r'), "fasta")]
    		per = args.percent
    		indexes = [num for num in range(matrix.shape[1]) if list(matrix[:,num]).count("-") > len(matrix[:,num]) * per]
    		matrix = np.delete(matrix, indexes, axis=1)
    		with open(workdir+"/"+file.split(".")[0]+".fasta", 'w') as fname:
        		for name, seq in zip(names, list(matrix)):
            			fname.write(">" + name + '\n' + ''.join(list(seq)) + '\n')
	os.chdir(current_dir)
	print('Done')

class nd:
    def __init__(self, dist):
        self.dist = dist
        
    def seq(self, seq):
        self.s = seq
    
    def p(self, p):
        self.p = p
        
    def children(self, children):
        self.children = children

    def name(self, name):
        self.name = name
        
    def names(self,names):
        self.names = names

class T():
    def __init__(self, root, y):
        self.root = root
        self.y = y
        super().__init__()

    #Choose random alignment (with seq_len more than self.x)
    def random_aln(self):
        import random
        import glob
        from Bio import AlignIO 
        from sys import stderr, exit

        maxn = 1000
        count = 0
        randomfile = random.choice([file for file in glob.glob(workdir + "/*.fasta")])
        aln = AlignIO.read(open(randomfile, 'r'), "fasta")
        while len(str(aln[0].seq)) < self.y and count < maxn:
            count += 1
            randomfile = random.choice([file for file in glob.glob(workdir + "/*.fasta")])
            aln = AlignIO.read(open(randomfile, 'r'), "fasta")
        if count == maxn:
            #stderr.write("Not enough alignments of length greater than {}\nAborted\n".format(self.x))
            stderr.write("Not enough alignments of length greater than {}\nAborted\n".format(self.y))
            exit(1)
        self.file = randomfile.replace("\\", "/")
        self.alignment_dict = {record.id:str(record.seq) for record in aln}
    
    #Make root sequence and pin columns   
    def root_seq(self):
        column_id = [random.randint(0, self.y) for _ in range(self.y)]
        a = AlignIO.read(open(self.file, 'r'), "fasta")
        names, seq = [],''
        for num in range(self.y):
            record = random.randint(0, len(a)-1)
            names.append(a[record].id)
            seq += str(a[record].seq)[column_id[num]] 
        self.columns = column_id
        self.fdseq = seq
        self.fdnames = names
        
    #Make all_dists and all_id lists
    def dists(self):
        import re
        all_id, all_dists, values, one_id = [], [], [], ''
        with open(self.file[:-6]+".dist", 'r') as fname:
            for string in fname:
                if re.search("[A-Za-z]+", string) != None:
                    all_id.append(one_id)
                    all_dists.append(values)
                    one_id, values = string.split()[0], string.split()[1:]
                else:
                    values += string.split()
            all_id.append(one_id)
            all_dists.append(values)
        self.all_id = all_id[1:]
        self.all_dists = all_dists[1:]
        self.check_point_dists = copy.deepcopy(self.all_dists)

    def build(self, dist, children):
        self.elem = nd(dist)
        self.root = self.elem
        self.elem.name = 'root'
        self.elem.p = None
        self.elem.children = children
        self.elem.seq = self.fdseq
        self.elem.names = self.fdnames
        parent = [self.elem]
        prev = self.root
        while parent != []:
            p = parent.pop(0)
            new_children = []
            
            if p.children != []:
                for child in p.children:
                    self.x = nd(child.dist)
                    self.x.name = child.name
                    self.x.children = child.children
                    self.x.p = p
                    new_children.append(self.x)
                    parent.append(self.x)
                p.children = new_children
            
            if p != self.root:
                dist = p.dist
                self.names = copy.deepcopy(p.p.names)
                #print(self.names)
                seq = ''
                names = []
                #less, more, less_index, more_index = 100, 100, [], [] 
                less_index, more_index = [], []
                self.all_dists = copy.deepcopy(self.check_point_dists)  

                def eq(x,y): 
                    if x + y > 0: return x/(x + y)
                    return 0.0
                for al_num in range(len(self.names)):
                    values = self.all_dists[self.all_id.index(self.names[al_num])]
                    less, more = 100, 100 
                    for i in range(len(values)):
                        if float(values[i]) <= dist and abs(float(values[i])-dist) < abs(float(less)-dist):
                            less = values[i]
                        elif float(values[i]) > dist and abs(float(values[i])-dist) < abs(float(more)-dist):
                            more = values[i]
                    if less != 100 and more != 100:
                        diff1 = (round(abs(dist-float(less)), 4), less)
                        diff2 = (round(abs(dist-float(more)),4), more)
                        b1,b2 = min(diff1[0], diff2[0]), max(diff1[0],diff2[0])
                        prob = eq(b1, b2)
                        our_dist = np.random.choice([b1,b2], 1, p = [1-prob, prob])
                        if our_dist == diff1[0]:
                            our_dist = less
                        elif our_dist == diff2[0]:
                            our_dist = more
                        #print(our_dist, dist, prob, 1-prob, less, more)
                    else:
                        if less != 100: our_dist = less
                        elif more != 100: our_dist = more
                    ind = random.choice([i for i, x in enumerate(values) if x == our_dist])
                    self.all_dists[self.all_dists.index(values)][ind] = '1000'
                    names.append(self.all_id[ind])
                    #print(self.all_id[ind])
                    #print(self.alignment_dict)
                    seq += self.alignment_dict[self.all_id[ind]][self.columns[al_num]]
                p.seq = copy.deepcopy(seq) 
                p.names = copy.deepcopy(names)
                #print(p.seq, '\t', p.dist)
# Alignments generation
print('Generating alignments..')
try:
	os.mkdir(outdir)
except:
	pass

os.chdir(workdir)
for tree in tree_files:
	with open(current_dir + "/" + tree, 'r') as fname:
		t = Tree(fname.readline().strip())
		#print('Current tree:\n', t)
		for i in range(args.number):
			with open(outdir + '/' + str(tree).split(".")[0] + "_" + str(i) + '.fasta', 'w') as out:
				my_tree = T(nd(None), args.length)
				my_tree.random_aln()
				my_tree.root_seq()
				my_tree.dists()
				my_tree.build(t.dist, t.children)
				output = my_tree.root.children
				while output != []:
					child = output.pop(0)
					if child.children:
						for i in child.children:
                					output.append(i)
					if child.name:
						out.write('>'+child.name + '\n' + child.seq + '\n')
os.chdir(current_dir)
print('Done')
			
