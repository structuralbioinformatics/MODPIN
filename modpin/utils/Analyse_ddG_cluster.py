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
import scipy.stats as st
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


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

def fileExist (file):
 return os.path.exists(file) and os.path.isfile(file)


def reduce_outliers_lstsq(y,x,m,c,t):
    yy=[]
    xx=[]
    d={}
    n = int(len(x)*t/100)
    out=0
    for i in xrange(len(x)):
        b = m*x[i] + c
        d.setdefault(i,(abs(y[i] - b)))
    for key, value in sorted(d.iteritems(), key=lambda (k,v): (v,k),reverse=True):
        out = out + 1
        if out > n:
           yy.append(y[key])
           xx.append(x[key])
        else:
           sys.stdout.write("Outlier: %10d \tY= %f \tX= %f \tRMSD= %f\n"%(key,y[key],x[key],value))
    return (yy,xx)



def parse_user_arguments(*args, **kwds):
    parser = argparse.ArgumentParser(
        description = 'Selection and energy comparison of  WT versus mutant forms of clusters of protein-protein interactions',
        epilog      = '@Oliva\'s lab 2018')
    parser.add_argument('-d', '--folder', dest = 'folder', action = 'store',
                        help = 'Input directory with MODPPI results')
    parser.add_argument('-g','--ddG',dest='ddg',action = 'store',default=None,
                        help = 'SKEMPI parsed file of affinities')
    parser.add_argument('-e', '--Energy_tested', dest = 'score', action = 'store',default="ddG_mean",
                        help = 'Energy to be used in the test')
    parser.add_argument('-o', '--output_directory', dest = 'outdir', action = 'store', default = 'selected_clusters',
                        help = 'Output directory (default is ModPPI_models)')
    parser.add_argument('-l','--label',dest='label',action = 'store',default='',
                        help = 'Label to identify the output "Compare_ddG_with_average_"')
#    parser.add_argument('-n', '--number_of_models', dest = 'nmodels', action = 'store', default = 1, type=int,
#                        help = 'Minimum number of models select the cluster for analysis (default is 1)')
    parser.add_argument('-i', '--interface_percentil_mutants', dest = 'ipercentil', action = 'store', default = 1.0, type=float,
                        help = 'Minimum percentage of models with mutation in the interface to select the cluster for analysis (default is 0.0)')
    parser.add_argument('-p', '--p_value_threshold', dest = 'p_value', action = 'store', default = 0.05, type=float,
                        help = 'P-value threshold on the comparison of WT-mutant to select the cluster for analysis (default is 0.05)')
    parser.add_argument('-r', '--cluster_rank', dest = 'rank', action = 'store', default = 1, type=int,
                        help = 'Cluster ranking by size (default is 1)')
    parser.add_argument('-lt','--outliers_threshold',dest='outlier_threshold',action = 'store',default=None,
                        help = 'Free Energy Perturbation Outliers energy (limits max and min of FEP)')
    parser.add_argument('-lp','--outliers_deviation',dest='outlier_percentil',action = 'store',default='0.0',
                        help = 'Percentage of Outliers (neglects a percentage of data, default is 0 and uses all)')
    parser.add_argument('-s','--show_canvas',dest='canvas',action = 'store_true',
                        help = 'Show canvas figure after execution (not shown by default)')
    parser.add_argument('-fi','--format_image',dest='format',action = 'store',default='png',
                        help = 'Format of the image in the output file: png(default), pdf, ps, eps and svg ')
    parser.add_argument('-k','--factor_lambda',dest='factor_lambda',action = 'store',default='1.0',type=float,
                        help = 'Factor factor_lambda to multiply the energies (default is 1) ')
    parser.add_argument('-c','--clean_null',dest='clean_null',action = 'store',default=None,
                        help = 'Clean ddG values near 0.0 because they are caused by errors (default is None, otherwise select a threshold near 0, i.e. 0.1)')
    parser.add_argument('-v', '--verbose', dest = 'show', action = 'store_true',
                        help = 'Flag for verbose mode (default is False)')
    options = parser.parse_args()
    return options



def main():

 options=parse_user_arguments()

 files       = os.listdir(options.folder)
 score       = options.score
# nmodels     = options.nmodels
 p_threshold = options.p_value
 percentil   = options.ipercentil
 outdir      = options.outdir
 input_dir   = options.folder
 rank        = options.rank
 verbose     = options.show
 factor      = options.factor_lambda
 label_out   = options.label

 if not os.path.exists(outdir): os.makedirs(outdir)

 if verbose: sys.stdout.write("Parsing %s\n"%(options.ddg))
 ddg_dict={}
 if not fileExist(options.ddg):
    sys.stderr.write("File %s not found\n"%options.ddg)
 fa=open(options.ddg,"r")
 for line in fa:
    if line.startswith("#"): continue
    form,ddg=line.split()
    ddg_dict.setdefault(form,float(ddg))
 fa.close()


 ppi_sets=set()
 for input_file in files:
   if input_file.endswith(".ppi"): 
      ppi_sets.add(input_file)

 ddg_real_dict={}
 ddg_pred_dict={}
 comparison={}
 comparison_all={}
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
    if verbose: print "Check file  %s"%os.path.join(input_dir,ddg_file)
    if fileExist(os.path.join(input_dir,ddg_file)):
      if verbose: print "Open %s"%os.path.join(input_dir,ddg_file)
      fd=open(os.path.join(input_dir,ddg_file),'r')
      for line in fd:
       pair=line.strip().split("versus")
       statistic=line.strip().split()
       if len(statistic)>1:
         if statistic[-1] == "Statistic": statistic_type=statistic[0]
       if len(pair)>1:
        ppia=pair[0].split()[0]
        ppib=pair[1].split()[0]
        pval_averages=pair[1].split("p-value")[1]
        pvalue=pval_averages.split()[0]
        average_preform_a=pval_averages.split()[1]
        average_inform_a=average_preform_a.split("(")[1].split("#")[0]
        average_form_a="::".join(average_inform_a.split("-"))
        average_form_a_r=average_inform_a.split("-")[1]+"::"+average_inform_a.split("-")[0]
        average_value_a=pval_averages.split()[2]
        average_preform_b=pval_averages.split()[3]
        average_inform_b=average_preform_b.split("(")[1].split("#")[0]
        average_form_b="::".join(average_inform_b.split("-"))
        average_form_b_r=average_inform_b.split("-")[1]+"::"+average_inform_b.split("-")[0]
        average_value_b=pval_averages.split()[4]
        p,z=ppia.split("-")
        q,conformer_a= z.split("#")
        cluster_a,interface_a= conformer_a.split("_")
        overlap_a=int(interface_a.lstrip("i"))
        x,z=ppib.split("-")
        y,conformer_b = z.split("#")
        cluster_b,interface_b= conformer_b.split("_")
        overlap_b=int(interface_b.lstrip("i"))
        if int(cluster_a) == int(cluster_b):
          diff = float(average_value_b)-float(average_value_a)
          ddg_a = 0.0
          ddg_b = 0.0
          root_a=average_form_a.split("_")[0]
          mut_a ="_".join([mx for mx in average_form_a.split("_")[1:]])
          root_b=average_form_b.split("_")[0]
          mut_b ="_".join([mx for mx in average_form_b.split("_")[1:]])
          interactors_all=wt[0]+"::"+wt[1]
          if ddg_dict.has_key(average_form_a): ddg_a = ddg_dict.get(average_form_a)
          if ddg_dict.has_key(average_form_a_r): ddg_a = ddg_dict.get(average_form_a_r)
          if ddg_dict.has_key(average_form_b): ddg_b = ddg_dict.get(average_form_b)
          if ddg_dict.has_key(average_form_b_r): ddg_b = ddg_dict.get(average_form_b_r)
          ddg = ddg_b - ddg_a
          comparison_all.setdefault((interactors_all,statistic_type),[]).append((cluster_b,overlap_b,overlap_a,average_form_b,average_form_a,line.strip(),float(pvalue),float(average_value_b),float(average_value_a),diff,ddg))
          if (p,q) == wt or (q,p) == wt :
           if overlap_b > 0:
            if x>y:  interactors= x+"::"+y
            else:    interactors= y+"::"+x
            diff=float(average_value_b)-float(average_value_a)
            if ddg_dict.has_key(average_form_b): ddg=ddg_dict.get(average_form_b)
            if ddg_dict.has_key(average_form_b_r): ddg=ddg_dict.get(average_form_b_r)
            comparison.setdefault((interactors,statistic_type),[]).append((overlap_b,cluster_b,line.strip(),float(pvalue),float(average_value_b),float(average_value_a),diff,ddg))
          if (x,y) == wt or (y,x) == wt  :
           if overlap_a >0:
            if p>q:  interactors= p+"::"+q
            else:    interactors= q+"::"+p 
            diff=float(average_value_a)-float(average_value_b)
            if ddg_dict.has_key(average_form_a): ddg=ddg_dict.get(average_form_a)
            if ddg_dict.has_key(average_form_a_r): ddg=ddg_dict.get(average_form_a_r)
            comparison.setdefault((interactors,statistic_type),[]).append((overlap_a,cluster_a, line.strip(),float(pvalue),float(average_value_a),float(average_value_b),diff,ddg))
      fd.close()

 for (file_name,analysis),result in comparison.iteritems(): 
  fo=open(os.path.join(outdir,"Select_Cluster_with_WT_"+score+"_"+file_name+"_"+analysis+".out"),"w") 
  if verbose: print "Write results in %s"%(os.path.join(outdir,"Select_Cluster_with_WT_"+score+"_"+file_name+"_"+analysis+".out"))     
  for overlap,cluster,data,pvalue,ave_b,ave_a,diff,ddg in sorted(result,key=lambda x: x[0],reverse=True):
    #print "Percentage mutation %s Cluster %s  Data %s \n"%(overlap,cluster,data)
    fo.write("Percentage of mutation in interface around %d %s   Cluster %s  P-value: %f Average-Mutant: %f Average-WT: %f Difference: %f ddG: %f Data %s \n"%(10*int(overlap),"%",cluster,pvalue,ave_b,ave_a,diff,ddg,data))
  fo.close()
         

 for (file_name,analysis),result in comparison_all.iteritems(): 
  fo=open(os.path.join(outdir,"Select_Cluster_all_mutants_"+score+"_"+file_name+"_"+analysis+".out"),"w") 
  if verbose: print "Write results in %s"%(os.path.join(outdir,"Select_Cluster_all_mutants_"+score+"_"+file_name+"_"+analysis+".out"))     
  for cluster,overlap_b,overlap_a,form_b,form_a,data,pvalue,ave_b,ave_a,diff,ddg in sorted(result,key=lambda x: int(x[0]),reverse=True):
    #print "Percentage mutation %s Cluster %s  Data %s \n"%(overlap,cluster,data)
    fo.write("Percentage of mutation in interface around %d %s and %d %s  Cluster %s  P-value: %f Average %s: %f Average %s: %f Difference: %f ddG: %f Data %s \n"%(10*int(overlap_b),"%",10*int(overlap_a),"%",cluster,pvalue,form_b,ave_b,form_a,ave_a,diff,ddg,data))
    if pvalue < p_threshold and int(cluster) <=rank:
        if (10*int(overlap_b) >= percentil or len(form_b.split("_"))==1) and (10*int(overlap_a) >= percentil or len(form_a.split("_"))==1):
            if form_b < form_a:
                ddg_real_dict.setdefault((form_b,form_a),ddg)
                ddg_pred_dict.setdefault((form_b,form_a),diff)
            else:
                ddg_real_dict.setdefault((form_a,form_b),ddg)
                ddg_pred_dict.setdefault((form_a,form_b),diff)
  fo.close()


 ddg_list=[]
 ave_list=[]
 for forms,ene in ddg_real_dict.iteritems():
   if ddg_pred_dict.has_key(forms):
      ddg_list.append(ene)
      ave_list.append(factor*ddg_pred_dict.get(forms))
      

 if len(ave_list)<=0:
     sys.stderr.write("Exit: this set is empty. Please change your conditions\n")
     exit(0)

 ddg_real=np.array(ddg_list)
 ddg_pred=np.array(ave_list)


 (pearson,prob)=st.pearsonr(ddg_real,ddg_pred)

 A = np.vstack([ddg_pred, np.ones(len(ddg_pred))]).T
 m, c = np.linalg.lstsq(A, ddg_real,rcond=-1 )[0]

 output=outdir+"/Compare_ddG_with_average_"+label_out+"_"+score+".dat"
 graph=outdir+"/Compare_ddG_with_average_"+label_out+"_"+score+"."+options.format
 datagraph=outdir+"/Compare_ddG_with_average_"+label+"_"+score+"."+options.format+".dat"

 
 
 if float(options.outlier_percentil) > 0:
   (ddg_list_reduced,ave_list_reduced)=reduce_outliers_lstsq(ddg_list,ave_list,m,c,float(options.outlier_percentil))
   ddg_real=np.array(ddg_list_reduced)
   ddg_pred=np.array(ave_list_reduced)
   (pearson,prob)=st.pearsonr(ddg_real,ddg_pred)
   A = np.vstack([ddg_pred, np.ones(len(ddg_pred))]).T
   m, c = np.linalg.lstsq(A, ddg_real,rcond=-1 )[0]


 fo=open(output,"w")
 fo.write("Pearson Correlation:\t%f\n"%pearson)
 fo.write("Two-tailed p-value: \t%f\n"%prob)
 fo.write("Slope of fitting:   \t%f\n"%m)
 fo.write("Cut Y-axis:         \t%f\n"%c)
 fo.close()

 fo=open(datagraph,"w")
 fo.write("%15s\t%15s\n"%("ddG predict","ddG real"))
 for i in xrange(0,len(ddg_pred)):
  fo.write("%15.5f\t%15.5f\n"%(ddg_pred[i],ddg_real[i]))
 fo.close()



 image, fig = plt.subplots()

 fig.plot(ddg_pred, ddg_real, 'o', label='SKEMPI ddG', markersize=2)
 fig.plot(ddg_pred, m*ddg_pred + c, 'r', label='Fitted line')
 fig.legend()

 fig.set(xlabel='ddG predicted (Kcal/mol)', ylabel='ddG experimental (Kcal/mol)',
       title='Compare ddG of SKEMPI and ddG predicted by MODPIN')
 image.savefig( graph)
 if options.canvas:
    plt.show()


  

# MAIN ###########################################################################################################

if __name__ == '__main__':
    main()   

