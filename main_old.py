import sys
import os
from FastaSeq import *
from FastaGroup import *

# Usage python main.py group_1_dir group_1_out_dir group_2_dir group_2_out_dir

def main():
	if not os.path.exists(sys.argv[2]):
		os.mkdir(sys.argv[2])# make dir to hold data for spearman
		group_1 = FastaGroup(sys.argv[1]) # load all fastas into FastaGroup object
		group_1.spearmanPrep(sys.argv[1], sys.argv[2]) # Prep data for Spearman
	
	else:
		group_1 = FastaGroup(sys.argv[1]) # load all fastas into FastaGroup object
		group_1.spearmanPrep(sys.argv[1], sys.argv[2]) # Prep data for Spearman


	if not os.path.exists(sys.argv[4]):
		os.mkdir(sys.argv[4])# make dir to hold data for spearman
		group_2 = FastaGroup(sys.argv[3]) # load all fastas into FastaGroup object
		group_2.spearmanPrep(sys.argv[3], sys.argv[4]) # Prep data for Spearman
	
	else:
		group_2 = FastaGroup(sys.argv[3]) # load all fastas into FastaGroup object
		group_2.spearmanPrep(sys.argv[3], sys.argv[4]) # Prep data for Spearman


if __name__ == "__main__":
	main()