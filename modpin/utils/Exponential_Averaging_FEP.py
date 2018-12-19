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
import numpy as np


# Add '.' to sys.path
src_path = os.path.abspath(os.path.dirname(__file__))
sys.path.append(src_path)
script_path = src_path+"/../scripts"
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

#Parameters
#Original Boltzman constant (kcal/molK):1.9872036e-3
kBo  = 1.9872036e-3
#Temperature
T    = 300.0
#Least Square lineal fitting using Affinity Benchmark
#Maff = 0.115
#Caff = -7.326
#Maff = 1.0
Maff = 0.0003


#Boltzman constant fixed
kB   = kBo/Maff



def fileExist (file):
 return os.path.exists(file) and os.path.isfile(file)

def distribution(E):
  mean=0.0
  sigma=0.0
  nn=len([x for x in E.itervalues() if x<0])
  n=len([x for x in E.itervalues()])
  minimum=0.0
  maximum=0.0
  if n>1: minimum=min([float(x) for x in E.itervalues()])
  if n>1: maximum=max([float(x) for x in E.itervalues()])
  if n<2 : return (mean,sigma,minimum,maximum)
  for m,ee in E.iteritems():
    if ee<0:
      mean =mean +float(ee)/nn
      sigma=sigma+float(ee)*float(ee)/n
  sigma = sigma - mean*mean
  if sigma > 0: 
     sigma = np.sqrt(sigma)
     return (mean,sigma,minimum,maximum)
  else:
     return (mean,0.0,minimum,maximum)

  
def partition_function(E):
 Z  = 0.0
 (mean,sigma,minimum,maximum)=distribution(E)
 p=minimum    
 #p=min([minimum+sigma,mean-3*sigma,minimum-(minimum/2)])    
 for m,ee in E.iteritems():
   e=(float(ee) - p)/(kB*T)
   try:
     if e > 1.0e+6: print "Out of limit Model %s is %f > 1.0e+6 \n"%(m,e)
     if e < 1.0e+6 :  Z = Z + np.exp(-e)
     if e >= 1.0e+6:  Z = Z 
     if e<0: print "Out of limit Model %s is  %f < 0\n"%(m,e)
   except  Exception as err:
      sys.stderr.write("ERROR: %s\n"%err)
      sys.stderr.write("ERROR: %s %s\n"%(m,str(e)))
      continue
 return (Z,p)

def average_ensemble (S,E):
  (Z,p)= partition_function( E )
  a    = 0.0
  for m,ee in E.iteritems():
      e=(float(ee) - p)/(kB*T)
      if S.has_key(m):
         try:
           if e < 1.0e+6:  a = a + float(S.get(m)) * np.exp( -e ) / Z
           if e >= 1.0e+6: a = a
         except  Exception as err:
           sys.stderr.write("ERROR: %s\n"%err)
           sys.stderr.write("ERROR: %s %s\n"%(m,str(e)))
           continue
  return a

def FEP (A,B):
  (ZA,p) = partition_function( A )
  x      = 0.0
  data   = [ (float(B.get(m)) - float(A.get(m))) for m,ee in A.iteritems() if B.has_key(m) ]
  if len(data) > 0: 
     d_min  = min(data)
  else:
     return 0.0
  for m,ee in A.iteritems():
    e=(float(ee) - p)/(kB*T)
    if B.has_key(m):
       try:
         d = float(B.get(m)) - float(A.get(m)) - d_min 
         if e < 1.0e+6 and d/(kB*T)< 1.0e+6 :  x = x + np.exp( -d/(kB*T) ) * np.exp(-e)/ZA
         if e >= 1.0e+6 or d/(kB*T)>=1.0e+6 :  x = x
       except  Exception as err:
         sys.stderr.write("ERROR: %s\n"%err)
         sys.stderr.write("ERROR: %s %s\n"%(m,str(e)))
         continue
  if x>0:
     fep = - kB * T * np.log(x) + d_min
  else:
     fep = d_min
  return fep
  

if len(sys.argv)==1:
 sys.stdout.write("   Execution is 'Exponential_Averaging_FEP.py <directory>    <variable>   <energy for Z>' \n \
 -- where directory is the output folder of modelist.py,  \n \
 -- 'energy for Z' is the score to calculate Z using Zwanzig equation (ddG_mean by default), and \n \
 -- 'variable' is the function to evaluate the comparison of exponential averages (ddG_mean by default)\n")
 exit(0)

input_dir=sys.argv[1]

if len(sys.argv)>3: 
  energy=sys.argv[3]
  score=sys.argv[2]
elif len(sys.argv)==3:
  energy="ddG_mean"
  score=sys.argv[2]
else: 
  energy="ddG_mean"
  score="ddG_mean"



files=os.listdir(input_dir)

ppi_sets=set()
for input_file in files:
   if input_file.endswith(".ppi"): 
      ppi_sets.add(input_file)

global_output = os.path.join(input_dir,"FEP_"+score+"_Z_with_"+energy+".out")
fg = open(global_output,"w")
fg.write("%5s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%s\n"%("#ppi","FEP direct","FEP invers","<FEP>","-kTln(Z'/Z)","score WT","score MUT","Diff.score","PPI models"))
for file_ppi in ppi_sets:
    output = os.path.join(input_dir,file_ppi+"_"+score+"_Z_with_"+energy+".out")
    print "Write Output in %s"%(output)
    fo = open(output,"w")
    fo.write("%5s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%s\n"%("#ppi","FEP direct","FEP invers","<FEP>","-kTln(Z'/Z)","score WT","score MUT","Diff.score","PPI models"))
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
    wt_file=data[1]+"::"+data[-1]
    if not os.path.exists(os.path.join(input_dir,"models",wt_file)): continue
    print "Check files in %s"%(wt_file)
    files_wt = os.listdir(os.path.join(input_dir,"models",wt_file))
    files_ddg_wt_path=set()
    files_scr_wt_path=set()
    number_of_clusters=set()
    for f in files_wt:
      if f.endswith("out"):
        cluster_group = f.split("cluster_")
        if len(cluster_group)>0:
           number_group=cluster_group[1].split(".list.")
           if len(number_group)>0:
              if number_group[0] != "0":
                 number_of_clusters.add(number_group[0])
    #print "Clusters %s"%number_of_clusters
    for number_group in number_of_clusters:
           ddg_wt_file=root_ppi+"_"+wt_file+"_cluster_"+number_group+".list."+energy+"."+root_ppi+".out"
           scr_wt_file=root_ppi+"_"+wt_file+"_cluster_"+number_group+".list."+score +"."+root_ppi+".out"
           files_ddg_wt_path.add(os.path.join(input_dir,"models",wt_file,ddg_wt_file))
           files_scr_wt_path.add(os.path.join(input_dir,"models",wt_file,scr_wt_file))
    E_wt={}
    S_wt={}
    for ddg_wt_path in files_ddg_wt_path:
       if not fileExist(ddg_wt_path): continue
       print "Open %s"%ddg_wt_path
       fd=open(ddg_wt_path,"r")
       for line in fd:
         model_file,ene = line.split()
         if model_file.endswith(".h") or model_file.endswith(".h_AB") or model_file.endswith(".h_BA"):
            model = os.path.basename(model_file)
            E_wt.setdefault(model,float(ene))
       fd.close()
    for scr_wt_path in files_scr_wt_path:
       if not fileExist(scr_wt_path): continue
       print "Open %s"%scr_wt_path
       fd=open(scr_wt_path,"r")
       for line in fd:
         model_file,scr = line.split()
         if model_file.endswith(".h") or model_file.endswith(".h_AB") or model_file.endswith(".h_BA"):
            model = os.path.basename(model_file)
            S_wt.setdefault(model,float(scr))
       fd.close()
    if len(files_ddg_wt_path)>0 and len(files_scr_wt_path)>0:
       n = 0;
       for ppi in ppis:
           type_file     = ppi[0]+"::"+ppi[1]
           if not os.path.exists(os.path.join(input_dir,"models",type_file)): continue
           print "Check files in %s"%(type_file)
           files_type = os.listdir(os.path.join(input_dir,"models",type_file))
           files_ddg_type_path=set()
           files_scr_type_path=set()
           number_of_clusters=set()
           for f in files_type:
               if f.endswith("out"):
                 cluster_group = f.split("cluster_")
                 if len(cluster_group)>0:
                    number_group=cluster_group[1].split(".list.")
                    if len(number_group)>0:
                       if number_group[0] != "0":
                          number_of_clusters.add(number_group[0])
           #print "Clusters %s"%number_of_clusters
           for number_group in number_of_clusters:
                    ddg_type_file=root_ppi+"_"+type_file+"_cluster_"+number_group+".list."+energy+"."+root_ppi+".out"
                    scr_type_file=root_ppi+"_"+type_file+"_cluster_"+number_group+".list."+score +"."+root_ppi+".out"
                    files_ddg_type_path.add(os.path.join(input_dir,"models",type_file,ddg_type_file))
                    files_scr_type_path.add(os.path.join(input_dir,"models",type_file,scr_type_file))
           E_type={}
           S_type={}
           for ddg_type_path in files_ddg_type_path:
               if not fileExist(ddg_type_path): continue
               print "Open %s"%ddg_type_path
               fd=open(ddg_type_path,"r")
               for line in fd:
                 model_file,ene = line.split()
                 if model_file.endswith(".h") or model_file.endswith(".h_AB") or model_file.endswith(".h_BA"):
                    model = os.path.basename(model_file)
                    E_type.setdefault(model,float(ene))
               fd.close()
           for scr_type_path in files_scr_type_path:
               if not fileExist(scr_type_path):continue
               print "Open %s"%scr_type_path
               fd=open(scr_type_path,"r")
               for line in fd:
                 model_file,scr = line.split()
                 if model_file.endswith(".h") or model_file.endswith(".h_AB") or model_file.endswith(".h_BA") :
                    model = os.path.basename(model_file)
                    S_type.setdefault(model,float(scr))
               fd.close()
           if len(files_ddg_type_path)>0 and len(files_scr_type_path)>0:
               n = n + 1
               Scr_wt       = average_ensemble( S_wt,E_wt )
               Scr_type     = average_ensemble( S_type,E_type )
               delta_scr    = Scr_type - Scr_wt
               FEP_direct   = FEP(E_wt,E_type)
               FEP_invers   = FEP(E_type,E_wt)
               FEP_average  = (FEP_direct - FEP_invers)/2
               (ZA,pA)      = partition_function( E_wt )
               (ZB,pB)      = partition_function( E_type )
               if ZA>0 and ZB>0: partition = - kB * T * np.log(ZB/ZA)
               fo.write("%5d\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%s\n"%(n,FEP_direct,FEP_invers,FEP_average,partition,Scr_wt,Scr_type,delta_scr,type_file))
               fg.write("%5d\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%s\n"%(n,FEP_direct,FEP_invers,FEP_average,partition,Scr_wt,Scr_type,delta_scr,type_file))
               print "Add %5d\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%10.5e\t%s\n"%(n,FEP_direct,FEP_invers,FEP_average,partition,ZA,ZB,Scr_wt,Scr_type,delta_scr,type_file)
    fo.close()
fg.close()         
print "Done"

    

