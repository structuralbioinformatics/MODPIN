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

def sum_of_tuples(a,b):
    x=[]
    size=max(len(a),len(b))
    for i in xrange(size):
        if i <  len(a) and i <  len(b): x.append(a[i]+b[i])
        if i >= len(a) and i <  len(b): x.append(b[i])
        if i <  len(a) and i >= len(b): x.append(a[i])
    return tuple(x)

def plotROC(a,out,title):
  p=[]
  q=[]
  r=[]
  fo=open(out+".dat","w")
  fo.write("#%9s\t%10s\t%10s\n"%("FPR","TPR","TNR"))
  for fpr in sorted(a.iterkeys()):
      p.append(fpr)
      tpr,tnr=a.get(fpr)
      fo.write("%10.5f\t%10.5f\t%10.5f\n"%(fpr,tpr,tnr))
      q.append(tpr)
      r.append(tnr)
  fo.close()
  x = np.array(p)
  y = np.array(q)
  z = np.array(r)
  plt.plot(x,y, ms=5, lw=2, alpha=0.7, mfc='cyan',label="TPR")
  plt.plot(x,z, ms=5, lw=2, alpha=0.7, mfc='firebrick',label="Specificity")
  plt.xlabel('false positive ratio')
  plt.ylabel('true positive ratio')
  plt.title("%s"%title)
  plt.grid(True)
  p=[]
  q=[]
  p.append(0.0)
  q.append(1.0)
  p.append(1.0)
  q.append(0.0)
  x = np.array(p)
  y = np.array(q)
  plt.plot(x,y, ms=5, lw=2, alpha=0.7, mfc='green',label="random")
  plt.legend(bbox_to_anchor=(0.85, 0.9), loc='upper left', borderaxespad=0.)
  plt.savefig(out)
  plt.close()

def plotSuccess(a,out,title):
  p=[]
  t=[]
  r=[]
  w=[]
  fo=open(out+".dat","w")
  fo.write("#%9s\t%10s\t%10s\t%10s\n"%("correlation","coverage","success","failure"))
  for correlation in sorted([float(c) for c in a.iterkeys() if c!= "neutral"]):
      p.append(float(correlation))
      coverage, success_ratio,fail_ratio = a.get(correlation)
      fo.write("%10s\t%10.5f\t%10.5f\t%10.5f\n"%(str(correlation),coverage,success_ratio,fail_ratio))
      t.append(coverage)
      r.append(success_ratio)
      w.append(fail_ratio)
  fo.close()
  x = np.array(p)
  y = np.array(t)
  z = np.array(r)
  q = np.array(w)
  plt.plot(x,y, ms=5, lw=2, alpha=0.7, mfc='cyan',label="coverage")
  plt.plot(x,z, ms=5, lw=2, alpha=0.7, mfc='firebrick',label="success")
  plt.legend(bbox_to_anchor=(0.15, 0.95), loc='upper left', borderaxespad=0.)
  plt.xlabel('Pearson correlation')
  plt.ylabel('ratio')
  plt.title("%s"%title)
  plt.grid(True)
  plt.savefig(out)
  plt.close()


def parse_user_arguments(*args, **kwds):
    parser = argparse.ArgumentParser(
        description = 'Prediction of gain/loss of interaction for  WT versus mutant forms using clusters of protein-protein interactions',
        epilog      = '@Oliva\'s lab 2018')
    parser.add_argument('-d', '--folder', dest = 'folder', action = 'store',
                        help = 'Input directory with MODPPI results')
    parser.add_argument('-s', '--skip', dest = 'skip', action = 'store',
                        help = 'file with PPIs that has to be skipped')
    parser.add_argument('-g','--ddG',dest='ddg',action = 'store',default=None,
                        help = 'SKEMPI parsed file of affinities to check TPR and FPR (if empty no analysis is done)')
    parser.add_argument('-e', '--Energy_tested', dest = 'score', action = 'store',default="ddG_mean",
                        help = 'Energy to be used in the test')
    parser.add_argument('-i', '--input_ranks', dest = 'file_ranks', action = 'store',default=None,
                        help = 'Input file with the ranking conditions of the best correlation between predicted and real ddG')
    parser.add_argument('-o', '--output_directory', dest = 'outdir', action = 'store', default = 'selected_clusters',
                        help = 'Output directory (default is ModPPI_models)')
    parser.add_argument('-l','--label',dest='label',action = 'store',default='',
                        help = 'Label to identify the output "Compare_ddG_with_average_"')
    parser.add_argument('-k','--factor_lambda',dest='factor_lambda',action = 'store',default='1.0',type=float,
                        help = 'Factor factor_lambda to multiply the energies (default is 1) ')
    parser.add_argument('-x', '--cross_ddG', dest = 'cross', action = 'store_true',
                        help = 'Flag to produce all crossed pairs of ddG (default is False)')
    parser.add_argument('-v', '--verbose', dest = 'show', action = 'store_true',
                        help = 'Flag for verbose mode (default is False)')
    options = parser.parse_args()
    return options

def predictor(ddg_prediction,ddg,epsilon):
  if   ddg_prediction >  epsilon :    state_prediction="loss"
  if   ddg_prediction < -epsilon :    state_prediction="gain"
  if  abs(ddg_prediction) < epsilon : state_prediction="neutral"
  ok   = 0
  fail = 0
  loss = (0, 0, 0, 0)
  gain = (0, 0, 0, 0)
  same = (0, 0, 0, 0)
  if ddg is not None:
     #      tp,fp,tn,fn
     loss = (0, 0, 0, 0)
     gain = (0, 0, 0, 0)
     same = (0, 0, 0, 0)
     if   ddg_prediction >  epsilon and ddg >  epsilon:      
          ok   = 1
          fail = 0
          #      tp,fp,tn,fn
          loss = (1, 0, 0, 0)
          gain = (0, 0, 1, 0)
          same = (0, 0, 1, 0)
     if   ddg_prediction < -epsilon and  ddg < -epsilon:      
          ok   = 1
          fail = 0
          #      tp,fp,tn,fn
          loss = (0, 0, 1, 0)
          gain = (1, 0, 0, 0)
          same = (0, 0, 1, 0)
     if  abs(ddg)<=epsilon and abs(ddg_prediction)<= epsilon: 
          ok   = 1
          fail = 0
          #      tp,fp,tn,fn
          loss = (0, 0, 1, 0)
          gain = (0, 0, 1, 0)
          same = (1, 0, 0, 0)
     if   ddg_prediction > epsilon and  ddg < epsilon:      
          ok   = 0
          fail = 1
          #      tp,fp,tn,fn
          loss = (0, 1, 0, 0)
          if ddg < -epsilon:
             gain = (0, 0, 0, 1)
          else:
             same = (0, 0, 0, 1)
     if   ddg_prediction < -epsilon and  ddg > -epsilon:      
          ok   = 0
          fail = 1
          #      tp,fp,tn,fn
          gain = (0, 1, 0, 0)
          if ddg > epsilon:
            loss = (0, 0, 0, 1)
          else:
            same = (0, 0, 0, 1)
     if  abs(ddg)>epsilon and abs(ddg_prediction)<= epsilon: 
          ok   = 0
          fail = 1
          #      tp,fp,tn,fn
          same = (0, 1, 0, 0)
          if ddg >  epsilon:
            loss = (0, 0, 0, 1)
          else:
            gain = (0, 0, 0, 1)
  print state_prediction
  print ddg_prediction,ddg
  return state_prediction,ok,fail, gain, loss, same

def test_equal(a,b):

    aa,ab=a.split("::")
    ba,bb=b.split("::")
    protein_A1=aa.split()[0]
    protein_A2=ab.split()[0]
    protein_B1=ba.split()[0]
    protein_B2=bb.split()[0]

    check=False

    if str(protein_A1) == str(protein_B1) and str(protein_A2) == str(protein_B2): check =True
    if str(protein_A1) == str(protein_B2) and str(protein_A2) == str(protein_B1): check =True
   
    return check
 

def tpr_fpr(tp,fp,tn,fn):
    tpr= fpr= tnr = 0.0
    if (tp+fp)>0: tpr = float(tp)/float(tp+fp)
    if (tn+fn)>0: fpr = float(fp)/float(tn+fn)
    if (tn+fn)>0: tnr = float(tn)/float(tn+fn)
    return tpr,tnr,fpr

def main():

 options=parse_user_arguments()

 files       = os.listdir(options.folder)
 skip_file   = os.listdir(options.skip)
 score       = options.score
 outdir      = options.outdir
 input_dir   = options.folder
 input_rank  = options.file_ranks
 verbose     = options.show
 factor      = options.factor_lambda
 label_out   = options.label
 cross       = options.cross

 if not os.path.exists(outdir): os.makedirs(outdir)
 if not  fileExist(input_rank):
   sys.stderr.write("Missing File with ranking conditions\n")
   exit(0) 
 if options.ddg is not None:
   if verbose: sys.stdout.write("Parsing %s\n"%(options.ddg))
   ddg_dict={}
   if not fileExist(options.ddg):
    sys.stderr.write("File %s not found \n"%options.ddg)
   fa=open(options.ddg,"r")
   for line in fa:
    if line.startswith("#"): continue
    form,ddg=line.split()
    ddg_dict.setdefault(form,float(ddg))
   fa.close()

 ppi_skip=set()
 fs=open(skip_file,"r")
 for line in fs:
   a,b=line.strip().split()
   ppi_skip.add((a,b))
   ppi_skip.add((b,a))
 fs.close()

 ppi_sets=set()
 for input_file in files:
   if input_file.endswith(".ppi"): 
      data=input_file.split("_")
      test_set=set()
      for data_x in data:
       for data_y in data:
           test_set.add((data_x,data_y))
      if len(test_set.intersection(ppi_skip))>0:continue      
      ppi_sets.add(input_file)

 conditions={}
 fa=open(input_rank,"r")
 for line in fa:
   if line.startswith("#"): continue
   rank, pmi, pv , cluster, correlation, nmodels, sign, slope, y_axis = line.split()
   conditions.setdefault(float(correlation),[]).append((float(sign),float(pmi),float(pv),int(cluster),float(slope),float(y_axis)))
 fa.close()


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
          if cross:
            if (ddg_dict.has_key(average_form_a) or ddg_dict.has_key(average_form_a_r)) and (ddg_dict.has_key(average_form_b) or ddg_dict.has_key(average_form_b_r)):
              if average_form_a == average_form_b or average_form_a_r == average_form_b or average_form_a == average_form_b_r or average_form_a_r == average_form_b_r:
                  ddg = None
              else:
                  ddg = ddg_b - ddg_a
              comparison_all.setdefault((interactors_all,statistic_type),[]).append((int(cluster_b),float(overlap_b),float(overlap_a),average_form_b,average_form_a,float(pvalue),float(average_value_b),float(average_value_a),diff,ddg))
          if (p,q) == wt or (q,p) == wt :
              diff=float(average_value_b)-float(average_value_a)
              if average_form_a == average_form_b or average_form_a_r == average_form_b or average_form_a == average_form_b_r or average_form_a_r == average_form_b_r:
                  ddg = None
              else:
                  if ddg_dict.has_key(average_form_b): ddg=ddg_dict.get(average_form_b)
                  if ddg_dict.has_key(average_form_b_r): ddg=ddg_dict.get(average_form_b_r)
              comparison_all.setdefault((interactors_all,statistic_type),[]).append((int(cluster_b),float(overlap_b),float(overlap_a),average_form_b,average_form_a,float(pvalue),float(average_value_b),float(average_value_a),diff,ddg))
          if (x,y) == wt or (y,x) == wt  :
              diff=float(average_value_a)-float(average_value_b)
              if average_form_a == average_form_b or average_form_a_r == average_form_b or average_form_a == average_form_b_r or average_form_a_r == average_form_b_r:
                  ddg = None
              else:
                  if ddg_dict.has_key(average_form_a): ddg=ddg_dict.get(average_form_a)
                  if ddg_dict.has_key(average_form_a_r): ddg=ddg_dict.get(average_form_a_r)
              comparison_all.setdefault((interactors_all,statistic_type),[]).append((int(cluster_b),float(overlap_a),float(overlap_b),average_form_a,average_form_b,float(pvalue),float(average_value_a),float(average_value_b),diff,ddg))
      fd.close()

 
 output =outdir+"/Predict_ddG_with_average_"+label_out+"_"+score+".dat"
 graph  =outdir+"/ROC_ddG_with_average_"+label_out+"_"+score
 outdata=outdir+"/ROC_ddG_with_average_"+label_out+"_"+score+".dat"

 epsilon= 0.5
 prediction={}
 tp_and_fp ={}
 okey  = 0
 wrong = 0
 #       tp,fp,tn,fn
 gain  = (0, 0, 0, 0)
 loss  = (0, 0, 0, 0)
 same  = (0, 0, 0, 0)
 for correlation in sorted(conditions.iterkeys(),reverse=True):
     if not conditions.has_key(correlation): continue
     parameters= conditions.get(correlation)
     step  = 0
     if verbose: sys.stdout.write("Checking correlation %s Gain= %s Loss= %s Neutral= %s\n"%(str(correlation),str(gain),str(loss),str(same)))
     for (sign,pmi,pv,cluster,slope,y_axis) in parameters:
       step = step + 1
       number_of_models = 0
       for (file_name,analysis),result in comparison_all.iteritems():
         if verbose and step==1 and number_of_models < 100: 
            sys.stdout.write("Checking Interaction %s with condition PMI %10.5f P-value %10.1e Cluster %d \n"%(file_name,float(pmi),float(pv),int(cluster)))
            number_of_models = number_of_models + 1
         if verbose and step==1 and number_of_models == 100: 
            sys.stdout.write("...continue with the rest of interactions\n")
            number_of_models = number_of_models + 1
         for check_cluster,overlap_b,overlap_a,form_b,form_a,pvalue,ave_b,ave_a,diff,ddg in sorted(result,key=lambda x: int(x[0]),reverse=True):
             if test_equal(form_b,form_a): continue
             if verbose and step==1 and number_of_models < 10:
                sys.stdout.write("\t-- Check features %s %s Cluster %s PMI_A %10s PMI_B %10s P-value %10s Predicted ddG: %10s Known ddG: %10s\n"%(form_a,form_b,str(check_cluster),str(overlap_a),str(overlap_b),str(pvalue),str(diff),str(ddg)))
             if prediction.has_key((form_a,form_b)) or prediction.has_key((form_b,form_a)): continue
             if (10*int(overlap_b) >= pmi or  len(form_b.split("_"))==1) and (10*int(overlap_a) >= pmi or len(form_a.split("_"))==1) and check_cluster<=cluster and pvalue < pv:
               ddg_prediction = factor*diff*slope + y_axis
               if verbose: sys.stdout.write("\t-- Found Interaction %s %s with condition PMI %10.5f P-value %10.1e Cluster %d with correlation %10s and significance %10.1e Gain= %s Loss= %s Neutral= %s\n"%(form_a,form_b,float(pmi),float(pv),int(cluster),str(correlation),float(sign),str(gain),str(loss),str(same)))
               if verbose: sys.stdout.write("\t\t--Diff %f\n"%diff)
               if verbose: sys.stdout.write("\t\t--ddg %f = %f * %f * %f + %f\n"%(ddg_prediction,factor,diff,slope,y_axis))
               state_prediction=None
               if   ddg_prediction >  epsilon :    state_prediction="loss"
               if   ddg_prediction < -epsilon :    state_prediction="gain"
               if  abs(ddg_prediction) < epsilon : state_prediction="neutral"
               if state_prediction is None: continue
               if ddg is not None and state_prediction is not None:
                 state_prediction, ok,fail, g, l, s = predictor(ddg_prediction,ddg,epsilon)
                 wrong= wrong + fail
                 okey = okey  + ok
                 gain=sum_of_tuples(gain,g)
                 loss=sum_of_tuples(loss,l)
                 same=sum_of_tuples(same,s)
               prediction.setdefault((form_b,form_a),(check_cluster,overlap_b,overlap_a,pvalue,ave_b,ave_a,correlation,sign,slope,y_axis,state_prediction,ddg_prediction,str(ddg),ok,fail,g,l,s))
     tp_and_fp.setdefault(correlation,(okey,wrong,gain,loss,same))

 size_data=float(okey+wrong)
 okey  = 0
 wrong = 0
 #       tp,fp,tn,fn
 gain  = (0, 0, 0, 0)
 loss  = (0, 0, 0, 0)
 same  = (0, 0, 0, 0)
 for (file_name,analysis),result in comparison_all.iteritems():
  for check_cluster,overlap_b,overlap_a,form_b,form_a,pvalue,ave_b,ave_a,diff,ddg in sorted(result,key=lambda x: int(x[0]),reverse=True):
    if prediction.has_key((form_a,form_b)) or prediction.has_key((form_b,form_a)): continue
    if test_equal(form_b,form_a): continue
    ddg_prediction = 0
    if ddg is not None:
      state_prediction, ok,fail, g, l, s = predictor(ddg_prediction,ddg,epsilon)
      wrong= wrong + fail
      okey = okey  + ok
      gain=sum_of_tuples(gain,g)
      loss=sum_of_tuples(loss,l)
      same=sum_of_tuples(same,s)
    prediction.setdefault((form_b,form_a),(check_cluster,overlap_b,overlap_a,pvalue,ave_b,ave_a,correlation,sign,slope,y_axis,state_prediction,ddg_prediction,str(ddg),ok,fail,g,l,s))
 tp_and_fp.setdefault("neutral",(okey,wrong,gain,loss,same))
               
 size_data_neutral=size_data+float(okey+wrong)

          
 fo=open(output,"w")
 fo.write("#%14s\t%15s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%20s\t%20s\t%20s\n"%("Protein_B","Protein_A","PMI_B","PMI_A","Pvalue","cluster","correlation","two-tail-prob","slope","Y-axis","dG_B","dG_A","Prediction","pred.ddG","ddG","success","fail","gain","loss","unaffected"))
 for pair,results in  prediction.iteritems():
     protein_b,protein_a = pair
     check_cluster,pmi_b,pmi_a,pvalue,dg_b,dg_a,correlation,sign,slope,y_axis,state_prediction,ddg_prediction,ddg,ok,f,g,l,s =  results
     fo.write("%15s\t%15s\t%10.5f\t%10.5f\t%10.1e\t%10d\t%10s\t%10.1e\t%10.5f\t%10.5f\t%10s\t%10s\t%10s\t%10.5f\t%10s\t%10d\t%10d\t%20s\t%20s\t%20s\n"%(protein_b,protein_a,pmi_b,pmi_a,pvalue,check_cluster,str(correlation),sign,slope,y_axis,dg_b,dg_a,state_prediction,ddg_prediction,ddg,ok,f,g,l,s))
 fo.close()

 
 tpr_tnr_fpr_gain={}
 tpr_tnr_fpr_loss={}
 tpr_tnr_fpr_same={}
 coverage_success={}

 fo=open(outdata,"w")
 fo.write("#%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n"%("Correlation","Success","Fails","Succ.Ratio","Fail Ratio","Coverage","TPR_gain","TNR_gain","FPR_gain","TP_gain","FP_gain","TN_gain","FN_gain","TPR_loss","TNR_loss","FPR_loss","TP_loss","FP_loss","TN_loss","FN_loss","TPR_neutral","TNR_neutral","FPR_neutral","TP_neutral","FP_neutral","TN_neutral","FN_neutral"))
 for correlation in sorted(tp_and_fp.iterkeys()):
     correct,fails,gain,loss,same = tp_and_fp.get(correlation)
     if correct+fails > 0: success_ratio = correct/float(correct+fails)
     if correct+fails > 0: fail_ratio    = fails/float(correct+fails)
     if correlation=="neutral": 
        if size_data_neutral<=0:continue
        coverage = float(correct+fails)/(size_data_neutral)
     else:
        if size_data<=0:continue
        coverage = float(correct+fails)/(size_data)
     tp_gain,fp_gain,tn_gain,fn_gain  = gain
     tp_loss,fp_loss,tn_loss,fn_loss  = loss
     tp_same,fp_same,tn_same,fn_same  = same
     tpr_gain,tnr_gain,fpr_gain = tpr_fpr(tp_gain,fp_gain,tn_gain,fn_gain)
     tpr_loss,tnr_loss,fpr_loss = tpr_fpr(tp_loss,fp_loss,tn_loss,fn_loss)
     tpr_same,tnr_same,fpr_same = tpr_fpr(tp_same,fp_same,tn_same,fn_same)
     tpr_tnr_fpr_gain.setdefault(fpr_gain,(tpr_gain,tnr_gain))
     tpr_tnr_fpr_loss.setdefault(fpr_loss,(tpr_loss,tnr_loss))
     tpr_tnr_fpr_same.setdefault(fpr_same,(tpr_same,tnr_same))
     coverage_success.setdefault(correlation,(coverage,success_ratio,fail_ratio))
     fo.write("%10s\t%10d\t%10d\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10d\t%10d\t%10d\t%10d\t%10.5f\t%10.5f\t%10.5f\t%10d\t%10d\t%10d\t%10d\t%10.5f\t%10.5f\t%10.5f\t%10d\t%10d\t%10d\t%10d\n"%(str(correlation),correct,fails,success_ratio,fail_ratio,coverage,tpr_gain,tnr_gain,fpr_gain,tp_gain,fp_gain,tn_gain,fn_gain,tpr_loss,tnr_loss,fpr_loss,tp_loss,fp_loss,tn_loss,fn_loss,tpr_same,tnr_same,fpr_same,tp_same,fp_same,tn_same,fn_same))
 fo.close()

 outroc=graph+"_gain_ROC.png"
 title="ROC gain interaction"
 plotROC(tpr_tnr_fpr_gain,outroc,title)

 outroc=graph+"_loss_ROC.png"
 title="ROC loss interaction"
 plotROC(tpr_tnr_fpr_loss,outroc,title)
 
 outroc=graph+"_neutral_ROC.png"
 title="ROC unaffected interaction"
 plotROC(tpr_tnr_fpr_same,outroc,title)
  
 outgraph=graph+"_success.png"
 title="Success,failure and coverage" 
 plotSuccess(coverage_success,outgraph,title)

  

# MAIN ###########################################################################################################

if __name__ == '__main__':
    main()   

