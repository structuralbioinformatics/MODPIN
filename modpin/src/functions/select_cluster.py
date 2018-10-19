import sys
import argparse
import os
import re
from Bio import SeqIO as SeqIE
from Bio import ExPASy
from Bio.Seq import Seq
from Bio import SeqRecord
from Bio.Alphabet import IUPAC
import ConfigParser
import itertools
import shutil
import subprocess
import time
import random


# Add '.' to sys.path
src_path = os.path.abspath(os.path.dirname(__file__))
sys.path.append(src_path)
script_path = src_path+"/../../scripts"
sys.path.append(script_path)

# Read configuration file
config = ConfigParser.ConfigParser()
config_file = os.path.join(script_path, 'config.ini')
config.read(config_file)

# Add SBI library path to sys.path
src_path = os.path.join(config.get('Paths','modppi_path'),config.get('Paths','sbi_library_path'))
sys.path.append(src_path)

# Add functions path to sys.path    
src_path = os.path.join(config.get('Paths','modppi_path'),config.get('Paths','functions_path'))
sys.path.append(src_path)

# Import functions
import functions
from SeqIO import *
from functions import *
from SBI.external.blast import blast_parser
from SBI.structure.contacts import Complex
from SBI.sequence import Sequence
from SBI.structure import PDB

def fileExist (file):
 return os.path.exists(file) and os.path.isfile(file)


if len(sys.argv)==1:
 sys.stdout.write("Excution is 'select_cluster <directory> <score>' where directory is the output folder of modelist.py and score is the score grouping of Rosetta to test (ddG_all by default)\n")
 exit(0)

input_dir=sys.argv[1]

if len(sys.argv)>2: score=sys.argv[2]
else: score="ddG_all"


files=os.listdir(input_dir)

ppi_sets=set()
for input_file in files:
   if input_file.endswith(".ppi"): 
      ppi_sets.add(input_file)

comparison={}
for file_ppi in ppi_sets:
    ppis=[]
    fd=open(os.path.join(input_dir,file_ppi),'r')
    for line in fd:
      pair=line.split()
      ppis.append((pair[0],pair[1]))
    fd.close()
    root_ppi=file_ppi.rstrip(".ppi")
    data=root_ppi.split("_")
    label=data[0]
    wt=(data[1],data[-1])
    ddg_file=root_ppi+"."+score+"_distribution_statistic.out"
    print "Check file  %s"%os.path.join(input_dir,ddg_file)
    if fileExist(os.path.join(input_dir,ddg_file)):
      print "Open %s"%os.path.join(input_dir,ddg_file)
      fd=open(os.path.join(input_dir,ddg_file),'r')
      for line in fd:
       pair=line.strip().split("versus")
       statistic=line.strip().split()
       if len(statistic)>1:
         if statistic[-1] == "Statistic": statistic_type=statistic[0]
       if len(pair)>1:
        ppia=pair[0].split()[0]
        ppib=pair[1].split()[0]
        p,z=ppia.split("-")
        q,conformer_a= z.split("#")
        cluster_a,interface_a= conformer_a.split("_")
        overlap_a=int(interface_a.lstrip("i"))
        x,z=ppib.split("-")
        y,conformer_b = z.split("#")
        cluster_b,interface_b= conformer_b.split("_")
        overlap_b=int(interface_b.lstrip("i"))
        if int(cluster_a) == int(cluster_b): 
          if (p,q) == wt or (q,p) == wt :
           if overlap_b > 0:
            if x>y:  interactors= x+"::"+y
            else:    interactors= y+"::"+x 
            comparison.setdefault((interactors,statistic_type),[]).append((overlap_b,cluster_b,line.strip()))
          if (x,y) == wt or (y,x) == wt  :
           if overlap_a >0:
            if p>q:  interactors= p+"::"+q
            else:    interactors= q+"::"+p 
            comparison.setdefault((interactors,statistic_type),[]).append((overlap_a,cluster_a, line.strip()))
      fd.close()

for (file_name,analysis),result in comparison.iteritems(): 
  fo=open(os.path.join(input_dir,"Select_Cluster_"+score+"_"+file_name+"_"+analysis+".out"),"w") 
  print "Write results in %s"%(os.path.join(input_dir,"Select_Cluster_"+score+"_"+file_name+"_"+analysis+".out"))     
  for overlap,cluster,data in sorted(result,key=lambda x: x[0],reverse=True):
    #print "Percentage mutation %s Cluster %s  Data %s \n"%(overlap,cluster,data)
    fo.write("Percentage of mutation in interface around %d %s   Cluster %s  Data %s \n"%(10*int(overlap),"%",cluster,data))
  fo.close()
         

    

