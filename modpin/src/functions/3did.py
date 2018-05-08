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
   dummy_dir=options.dummy_dir
   try:
    did_path=os.path.join(config.get('Paths','modppi_path'),config.get('Paths','3did_path'))
    data_path=os.path.join(config.get('Paths','modppi_path'),config.get('Paths','data_path'))
   except:
    did_path=options.outdir
    data_path=options.outdir

   if not os.path.exists(did_path):
    os.makedirs(did_path)
   if not os.path.exists(dummy_dir):
    os.makedirs(dummy_dir)
   if not os.path.exists(data_path):
    sys.stderr.write("No DATA directory, please check your installation or INPUT\n")

   #Parse did flat file
   did=parse_3did(options)

   #Create PDB files of 3DiD interactions
   for dd,cases in did.iteritems():
      for label in xrange(0,len(cases)):
        #Define the name of the PDB output file with domain-domain interactions
        did_file=os.path.join(did_path,dd[0]+":"+dd[1]+"#"+str(label)+".brk.gz")
        if not os.path.exists(did_file.lower()):
          did_file=os.path.join(did_path,dd[0]+":"+dd[1]+"#"+str(label)+".brk")
        if not os.path.exists(did_file.lower()):
          if verbose: sys.stderr.write("\t\t--Create %s\n"%(did_file.lower()))
          pdb_code,d1,d2 = cases[label]
          pdb_file=os.path.join(pdb_path,pdb_code[1:3].lower(),"pdb"+pdb_code+".ent")
          if not os.path.exists(pdb_file):
             pdb_file=os.path.join(pdb_path,pdb_code[1:3].lower(),"pdb"+pdb_code+".ent.gz")
          if not os.path.exists(pdb_file):
            if verbose: sys.stderr.write("\t\t\t-- %s not found\n"%pdb_file)
            continue
          try:
            pdb=PDB(pdb_file)
            brk=PDB()
            pdb_chain_A=pdb.get_chain_by_id(d1[0])
            start_A=d1[1]
            end_A  =d1[2]
            pdb_chain_B=pdb.get_chain_by_id(d2[0])
            start_B=d2[1]
            end_B  =d2[2]
            brk_chain_A=pdb_chain_A.extract(init=start_A,end=end_A)
            brk_chain_A.chain="A"
            brk.add_chain(brk_chain_A)
            brk_chain_B=pdb_chain_B.extract(init=start_B,end=end_B)
            brk_chain_B.chain="B"
            brk.add_chain(brk_chain_B)
            brk.clean()
            brk.write(did_file.lower())
          except Exception as e:
            if verbose: sys.stderr.write("\t\t\t  Error: %s\n"%e)
            continue
          

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
        description = 'Generate the datase of DDI with the 3DiD flat file',
        epilog      = '@Oliva\'s lab 2016')
    parser.add_argument('-i', '--3did', dest = 'flat_file', action = 'store',
                        help = 'Input file: 3did_flat file ', default = '3did_flat')
    parser.add_argument('-ppi', '--interactions', dest = 'interactions_file', action = 'store',
                        help = 'Output file with a list of the interactions from 3did')
    parser.add_argument('-seq', '--fasta', dest = 'seq_file', action = 'store',
                        help = 'Output file with FASTA sequences of the proteins from 3DiD')
    parser.add_argument('-o', '--output_directory', dest = 'outdir', action = 'store', default = '3did',
                        help = 'Output directory (default is 3did)')
    parser.add_argument('-d', '--dummy_dir', dest = 'dummy_dir', action = 'store', default = '/tmp/3did_dummy',
                        help = 'Location of the dummy folder (default is /tmp/3did_dummy)')
    parser.add_argument('-v', '--verbose', dest = 'show', action = 'store_true',
                        help = 'Flag for verbose mode (default is False)')
    options = parser.parse_args()
    return options

def parse_3did(options):

 did={}
 flat=open(options.flat_file,'r')
 for line in flat:
   if line.startswith("#=ID"):
      word=line.split()
      pf1=word[1]
      pf2=word[2]
      ddi=(pf1,pf2)
   if line.startswith("#=3D"):
      word=line.split()
      pdb=word[1]
      frame1=word[2]
      frame2=word[3]
      chain1=frame1.split(":")[0]
      chain2=frame2.split(":")[0]
      section1=frame1.split(":")[1]
      section2=frame2.split(":")[1]
      start1=section1.split("-")[0]
      start2=section2.split("-")[0]
      end1=section1.split("-")[1]
      end2=section2.split("-")[1]
      d1=(chain1,start1,end1)
      d2=(chain2,start2,end2)
      # Stored information in did: PDB code, and two domain fragments. 
      # Each domain fragment is defined by Chain and Aas interval (start-end)
      triad=(pdb,d1,d2)
      did.setdefault(ddi,[]).append(triad)

 return did


if __name__ == '__main__':
    main()


