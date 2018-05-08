import sys
import argparse
import os
import re
import ConfigParser
import shutil
import subprocess
import time
import random
from collections import Counter

# Add '.' to sys.path
src_path = os.path.abspath(os.path.dirname(__file__))
sys.path.append(src_path)

# Read configuration file
config = ConfigParser.ConfigParser()
config_file = os.path.join(src_path, 'config.ini')
config.read(config_file)

# Add SBI library path to sys.path
src_path = os.path.join(config.get('Paths','modppi_path'),config.get('Paths','sbi_library_path'))
sys.path.append(src_path)

# Add functions path to sys.path
src_path = os.path.join(config.get('Paths','modppi_path'),config.get('Paths','functions_path'))
sys.path.append(src_path)

# Import functions
from SeqIO import *
from functions import *

def main():
    
    #initialize
    options    = parse_user_arguments()
    verbose    = options.show
    dummy_dir  = options.dummy_dir
    if not os.path.exists(options.dummy_dir): os.makedirs(options.dummy_dir)
    pdb_file   = options.pdb_file
    fasta_file = options.fasta_file
    chains_file= options.chains_file
    hydrogens  = options.hbplus
    outfile    = options.outfile

    if not fileExist(pdb_file) or not fileExist(fasta_file) or not fileExist(chains_file):
       sys.stderr.write('EXIT: Missing files, please check you have the files of PDB, FASTA and list of chains\n')
       exit(0)

    path_pdb = "/".join(os.path.abspath(pdb_file).split("/")[0:-1])
    pdb_name = pdb_file.split("/")[-1]

    proteins= {}
    if verbose: sys.stdout.write('\t\t\t-- Reading protein sequences...\n')
    for protein in FASTA_iterator(fasta_file):
            name = protein.get_identifier()
            if len(name.split("|"))>2:
               name1 = name.split("|")[1] 
               name2 = name.split("|")[2] 
            if len(name.split("|")) == 2:
               name1 = name.split("|")[1]
               name2 = name.split("|")[1]
            if len(name.split("|")) == 1:
               name1 = name.split("|")[0]
               name2 = name.split("|")[0]
            seq  = protein.get_sequence()
            if len(seq) > 0:
                try:
                    proteins.setdefault(name1, ProteinSequence(name1, seq))
                    proteins.setdefault(name2, ProteinSequence(name2, seq))
                except IncorrectSequenceLetter as e:
                    sys.stderr.write('WARNING: %s\n' %e)
                    sys.stdout.write('\t\t\t-- Skip input sequence: %s Sequence: %s\n' %(name, seq))

    sequences={}
    fd=open(chains_file,"r")
    for line in fd:
       chain,fasta_name = line.strip().split()
       if proteins.has_key(fasta_name):
          sequences.setdefault(chain,proteins[fasta_name])
          #print "%s\t%s\t%s\t%s\n"%(chain,fasta_name,proteins[fasta_name].get_identifier(),proteins[fasta_name].get_sequence())
       else:
          sys.stderr.write('EXIT: Missing FASTA sequence %s\n'%fasta_name)
          exit(0)
    fd.close()
   
    if hydrogens:
       pdb_hydro=".".join(pdb_name.split(".")[0:-1])+".h"
       if not fileExist(os.path.join(path_pdb,pdb_hydro)):
            try:
               if verbose: sys.stdout.write("\t\t\t-- Adding hydrogens and relaxing the model %s\n"%pdb_name)
               add_hydrogens(config,path_pdb,pdb_name, pdb_hydro)
            except ValueError as e:
               sys.stderr.write("WARNING %s\n"%e)
               shutil.copy(pdb_file, os.path.join(path_pdb,pdb_hydro))
       else:
            if verbose: sys.stdout.write("\t\t\t-- The model has already hydrogens %s\n"%pdb_hydro)
       input_renumber=pdb_hydro
    else:
       input_renumber=pdb_name

    if not fileExist(outfile):
     pdb_new=PDB()
     pdb_new=renumber_pdb(config,path_pdb,input_renumber,sequences,dummy_dir)
     pdb_new.write(outfile)
    else:
     if verbose: sys.stdout.write("\t\t\t-- File %s already exist\n"%outfile)
        
def parse_user_arguments(*args, **kwds):
    parser = argparse.ArgumentParser(
        description = 'Automatic Renumbering of a protein with the original sequence(s)',
        epilog      = '@Oliva\'s lab 2016')
    parser.add_argument('-l', '--list', dest = 'chains_file', action = 'store',
                        help = 'Input file with a list of pairs "CHAIN FASTA_ID" to associate the chain name in the PDB file and the name of the sequence')
    parser.add_argument('-f', '--fasta', dest = 'fasta_file', action = 'store',
                        help = 'Input file with a list of FASTA sequences of the proteins to use')
    parser.add_argument('-p', '--pdb', dest = 'pdb_file', action = 'store',
                        help = 'Input file with a list of FASTA sequences of the proteins to use')
    parser.add_argument('-o', '--out', dest = 'outfile', action = 'store',
                        help = 'PDB output file renumbered')
    parser.add_argument('-d', '--dummy_dir', dest = 'dummy_dir', action = 'store', default = '/tmp/modppi_dummy',
                        help = 'Specifies the location of the dummy folder (default is /tmp/modppi_dummy)')
    parser.add_argument('-hydro','--hydrogens', dest = 'hbplus', action = 'store_true',
                        help = 'Flag to include hydrogens (default is False)')
    parser.add_argument('-v', '--verbose', dest = 'show', action = 'store_true',
                        help = 'Flag for verbose mode (default is False)')
    options = parser.parse_args()
    return options
    
# MAIN ###########################################################################################################

if __name__ == '__main__':
    main()
