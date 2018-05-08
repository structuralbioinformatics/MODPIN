import sys
import argparse
import os
import re
import ConfigParser
import shutil
import subprocess
import gzip

# Add '.' to sys.path
src_path = os.path.abspath(os.path.dirname(__file__))+"/../../scripts"
sys.path.append(src_path)

# Read configuration file
config = ConfigParser.ConfigParser()
config_file = os.path.join(src_path, 'config.ini')
config.read(config_file)

# Add SBI library path to sys.path
src_path = os.path.join(config.get('Paths','modppi_path'),config.get('Paths', 'sbi_library_path'))
sys.path.append(src_path)

# Add functions path to sys.path
src_path = os.path.join( config.get('Paths','modppi_path'), config.get('Paths', 'functions_path') )
sys.path.append(src_path)

# Import functions

from SBI.structure import PDB
import SBI.structure.contacts as cn
from functions import *
from SeqIO import *


def main():
   #Initialize
   options=parse_user_arguments()
   verbose=options.show
   pdb_path=os.path.join( config.get('Paths','modppi_path'),config.get('Paths', 'pdb_path') )
   try:
    did_path=os.path.join(config.get('Paths','modppi_path'),config.get('Paths','3did_path'))
    data_path=os.path.join(config.get('Paths','modppi_path'),config.get('Paths','data_path'))
   except:
    did_path=options.outdir
    data_path=options.outdir

   if not os.path.exists(did_path):
    sys.stderr.write("No 3DID directory, please check your installation or INPUT\n")
    
   if not os.path.exists(data_path):
    sys.stderr.write("No DATA directory, please check your installation or INPUT\n")
     

   #Create list of interactions and FASTA sequences of 3DiD
   did_interactions=open(os.path.join(data_path,options.interactions_file),"w")
   did_fasta       =open(os.path.join(data_path,options.seq_file),"w")
   for brk in os.listdir(did_path):
       if verbose: sys.stderr.write("\t\t-- Reading %s  \n"%os.path.join(did_path,brk))
       try:
          pdb=PDB(os.path.join(did_path,brk))
          id_chain=[]
          for c in pdb.chain_identifiers:
           pdb_chain=pdb.get_chain_by_id(c)
           id_chain.append(pdb.id+"_"+c)
           printfasta( did_fasta, pdb.id+"_"+c , pdb_chain.gapped_protein_sequence)      
          did_interactions.write("%s\t%s\n"%(id_chain[0],id_chain[1]))
       except Exception as e:
          if verbose: sys.stderr.write("\t\t-- %s cannot be read\n\t\t   Error: %s\n"%(os.path.join(did_path,brk),e))
          continue
   did_interactions.close()
   did_fasta.close()


def parse_user_arguments(*args, **kwds):
    parser = argparse.ArgumentParser(
        description = 'Generate the gapped sequences and list of DDIs with the 3DiD structures',
        epilog      = '@Oliva\'s lab 2016')
    parser.add_argument('-ppi', '--interactions', dest = 'interactions_file', action = 'store',
                        help = 'Output file with a list of the interactions from 3did (default is 3did.ddi)')
    parser.add_argument('-seq', '--fasta', dest = 'seq_file', action = 'store',
                        help = 'Output file with FASTA sequences of the proteins from 3DiD (default is 3did_gapped.fasta) ')
    parser.add_argument('-v', '--verbose', dest = 'show', action = 'store_true',
                        help = 'Flag for verbose mode (default is False)')
    options = parser.parse_args()
    return options



if __name__ == '__main__':
    main()


