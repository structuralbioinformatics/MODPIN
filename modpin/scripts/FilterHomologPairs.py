import sys
import os
from Bio import PDB
import numpy as np
import scipy as sci
from scipy import spatial
import string
import argparse

def main():
 options = parse_user_arguments()
 try:
  pdbdir,interface,cluster=filter_homologs(options.input,options.potential,options.radius,options.overlap,options.show)
  #Write outputs
  fg=open(options.output+"_list.dat","w")
  for c,s in cluster.iteritems():
   pdbroot  =c.split(".")
   seed     ="_".join([str(x) for x in pdbroot[:-1]])
   seed_list=options.output+"_cluster_"+seed+".dat"
   fd=open(seed_list,"w")
   grouped="\t"
   for p in s:
    pdb=pdbdir+"/"+p
    fd.write("%s\n"%(pdb))
    if p != c:
      z=0.0
      if len(interface[c]) >0:
        z=100 * float(len(interface[c].intersection(interface[p])))/len(interface[c])
      pb=p.split(".")
      grouped += ".".join([str(x) for x in pb[1:-1]])
      grouped += " ({0:5.1f})\t".format(z)
   fd.close()
   fg.write("%s\t%s\n"%(pdbdir+"/"+c,grouped))
  fg.close()
 except IOError as e:
  print("I/O error (%s): %s" %(e.errno, e.strerror))

def parse_user_arguments(*args, **kwds):
 parser = argparse.ArgumentParser(
  description = "Cluster Similar Interactions from a list of PDB binary interactions (pairs)",
  epilog      = "@oliva's lab 2014")
 parser.add_argument('-l','--list_input_file',dest='input',action = 'store',default='input_list',
  help = 'Input file (default is input_list)')
 parser.add_argument('-o','--filtered_output_file',dest='output',action = 'store',default='output',
  help = 'Root-name for output files (default is "output")')
 parser.add_argument('-c','--distance_radius_cutoff',dest='radius',type=float,action = 'store',default='5.0',
  help = 'Distance cut-off to calculate contact residue-pairs of the interface (default 5.0 Ang)')
 parser.add_argument('-pot','--contact_atoms',dest='potential',action = 'store',default='min',
  help = 'Atoms to calculate the distance "min": any closest atom; "cb": CB-CB; "ca": CA-CA')
 parser.add_argument('-id','--percentage_identity',dest='overlap',action = 'store',type=float,default='20.0',
  help = 'Minimum percentage of identical or overlaping pairs of contacts in clusters (default 20.0)')
 parser.add_argument('-v','--verbose',dest='show',action = 'store_true',
  help = 'Verbose execution')

 options=parser.parse_args()
 return options

def fileExist (file):
 return os.path.exists(file) and os.path.isfile(file)

def printverbose(f,flag,message):
 """Define verbose to print in file 'f', if 'flag'=True, message given as 'message'"""
 if flag: f.write("%s"%(message))



def filter_homologs(input_file,potential="min",radius=5.0,overlap=20.0,show=False):

 pdbdir=""
 if fileExist(input_file):
  inp=open(input_file,"r")
  pdblist=[]
  pdbdict={}
  for line in inp:
   word=line.strip().split()
   pdbpath=word[0].split("/")
   pdbfile=pdbpath[-1]
   pdbdir="/".join([str(x) for x in pdbpath[:-1]])
   pdbroot=pdbfile.split(".")
   pdbdict.setdefault(pdbroot[0],[]).append(pdbfile)
   printverbose(sys.stdout,show,("Add in the  list[%s]: %s\n"%(pdbroot[0],pdbfile)))
   pdblist.append(pdbfile)
  inp.close()
 else:
  raise IOError("Input File not found")

 if pdbdir=="": raise IOError("Input File is empty")
 printverbose(sys.stdout,show,("PDB directory: %s\n"%(pdbdir)))

 interface={}
 for pdb in pdblist:
   try:
     pdbfile=pdbdir+"/"+pdb
     printverbose(sys.stdout,show,("PDB file checking: %s\n"%(pdbfile)))
     interface.setdefault(pdb,set()).update(contact(potential,radius,show,pdbfile))
   except Exception as e:
     print("Error: %s" %(e))

 pdbset=set(pdblist)
 cluster={}
 for pa in pdblist:
  if pa in pdbset:
   cluster.setdefault(pa,set()).add(pa)
   printverbose(sys.stdout,show,("Create Cluster centered on: %s\n"%(pa)))
   pdbset.remove(pa)
   for pb in pdblist:
    if pb in pdbset:
     z=0.0
     if len(interface[pa]) >0:
      z=100 * float(len(interface[pa].intersection(interface[pb])))/len(interface[pa])
     if z >= overlap:
      printverbose(sys.stdout,show,("Add in cluster[%s]: %s (id=%10.2f)\n"%(pa,pb,z)))
      cluster[pa].add(pb)
      pdbset.remove(pb)
 return (pdbdir,interface,cluster)

def contact(potential,radius, show,pdbfile):

 if fileExist(pdbfile):
  parser= PDB.PDBParser(QUIET=True)
  fp=open(pdbfile,"r")
  structure=parser.get_structure("complex",fp)
 else:
  raise IOError("File %s not found"%(pdbfile))

 if potential.lower() == "min":
  atoms=[]
  for a in structure.get_atoms():
   atoms.append(a)
  pairs=neighbourhood(atoms,radius,show)
 elif potential.lower() == "cb":
  atoms=[]
  for a in structure.get_atoms():
   r=PDB.Selection.unfold_entities([a], 'R')
   if a.get_name() == 'CB' or (a.get_name() == 'CA' and  r[0].get_resname() =='GLY'):
    atoms.append(a)
  pairs=neighbourhood(atoms,radius,show)
 elif potential.lower() == "ca":
  atoms=[]
  for a in structure.get_atoms():
   if a.get_name() == 'CA' :
    atoms.append(a)
  pairs=neighbourhood(atoms,radius,show)
 else:
  atoms=[]
  for a in structure.get_atoms():
   atoms.append(a)
  pairs=neighbourhood(atoms,radius,show)

 residue_pairs=set()
 for r1,a,r2,b in pairs:
   res1=str(r1.get_id()[1])+"_"+str(r1.get_resname())
   res2=str(r2.get_id()[1])+"_"+str(r2.get_resname())
   residue_pairs.add((res1,res2))
   residue_pairs.add((res2,res1))
 return residue_pairs


def neighbourhood(atoms,radius,show):
 pairs=set()
 for r1,r2 in PDB.NeighborSearch(atoms).search_all(radius,'A'):
   chain_list_1 = PDB.Selection.unfold_entities([r1], 'C')
   chain_list_2 = PDB.Selection.unfold_entities([r2], 'C')
   res_list_1 = PDB.Selection.unfold_entities([r1], 'R')
   res_list_2 = PDB.Selection.unfold_entities([r2], 'R')
   if chain_list_1[0].get_id() != chain_list_2[0].get_id():
      printverbose(sys.stdout,show,(" %s %s (Atom %s) (chain %s)  contact with %s %s (Atom %s) (chain %s) \n"%
            (res_list_1[0].get_resname(),res_list_1[0].get_id()[1],r1.get_id(),chain_list_1[0].get_id(),
             res_list_2[0].get_resname(),res_list_2[0].get_id()[1],r2.get_id(),chain_list_2[0].get_id())))
      pairs.add( ( res_list_1[0], chain_list_1[0].get_id(), res_list_2[0], chain_list_2[0].get_id() ) )
      pairs.add( ( res_list_2[0], chain_list_2[0].get_id(), res_list_1[0], chain_list_1[0].get_id() ) )
 return pairs


if __name__=="__main__":
  main()

