import sys
import argparse
import os
import ConfigParser
import gzip
import shutil
import subprocess

# Add "." to sys.path #
src_path =  os.path.abspath(os.path.dirname(__file__))+"/../../scripts/"
sys.path.append(src_path)

# Read configuration file
config = ConfigParser.ConfigParser()
config_file = os.path.join(src_path, 'config.ini')
config.read(config_file)

# Add SBI library path to sys.path
src_path = config.get('Paths', 'sbi_library_path')
sys.path.append(src_path)

# Add functions path to sys.path
src_path = config.get('Paths', 'functions_path')
sys.path.append(src_path)

# Import functions

from SBI.structure import PDB
import SBI.structure.contacts as cn
from functions import *



def main():

    options = parse_user_arguments()
    pdb2pin(options)

def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Create a network of PPI's extracted from PDB complexes",
        epilog      = "@oliva's lab 2015")
    parser.add_argument('-i','--list_pdb_files',dest='listfile',action = 'store',
                        help = 'List of PDB files with the address location')
    parser.add_argument('-o','--output_names',dest='out',action = 'store', default='output',
                        help = 'Output rootname for files with nodes and edges(default is output)')
    parser.add_argument('--PPI_distance',dest='PPI_distance',action = 'store', type=float, default=12.0,
                        help = 'Threshold distance to calculate contacts')
    parser.add_argument('--PPI_type',dest='PPI_type',action = 'store', default='cb',
                        help = 'Type of contact: ca, cb or min')
    parser.add_argument('-v','--verbose',dest='show',action = 'store_true',
                        help = 'Verbose execution ')
    options=parser.parse_args()

    return options

def pdb2pin(options):

 filelist=options.listfile
 rootname=options.out
 verbose=options.show
 PPI_distance = options.PPI_distance
 PPI_type = options.PPI_type
 
 edges=[]
 nodes=[]

 if fileExist(filelist):
  for line in parse_list_file(filelist):
   pdb=PDB(line)
   complexes=cn.Complex(pdb,PPI_distance = PPI_distance, PPI_type = PPI_type)
   for pair in complexes.PPInterfaces:
     #if not pair.is_empty:
     if len(pair.contacts)>0:
      edges.append((pair.interactor_id,pair.protein_id))
      if verbose: sys.stdout.write("Add {0:s}\t{1:s}\n".format(pair.interactor_id,pair.protein_id))
  nodes=get_nodes(edges)
 else:
  sys.stderr.write("Missing list of PDB files\n")
  exit(0)

 ppi=open(rootname+".ppi","w")
 nds=open(rootname+".dat","w")
 for (x,y) in edges:
  ppi.write("%s\t%s\n"%(x,y))
 for x in nodes:
  nds.write("%s\n"%(x))



if  __name__ == "__main__":
    main()



