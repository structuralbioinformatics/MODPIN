import os,sys,re
import argparse
import ConfigParser
import shutil
import string
import numpy as np
import scipy.stats as st
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def fileExist (file):
 if file is not None:
  return os.path.exists(file) and os.path.isfile(file)
 else:
  return False

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
        description = "Parse SKEMPI and create MODPIN files",
        epilog      = "@oliva's lab 2018")
    parser.add_argument('-a','--ddG',dest='ddg',action = 'store',default=None,
                        help = 'SKEMPI parsed file of affinities')
    parser.add_argument('-b','--FEP',dest='fep',action = 'store',default=None,
                        help = 'Free Energy Perturbation file obtained with MODPPI models')
    parser.add_argument('-lt','--outliers_threshold',dest='outlier_threshold',action = 'store',default=None,
                        help = 'Free Energy Perturbation Outliers energy (limits max and min of FEP)')
    parser.add_argument('-lp','--outliers_deviation',dest='outlier_percentil',action = 'store',default='0.0',
                        help = 'Percentage of Outliers (neglects a percentage of data, default is 0 and uses all)')
    parser.add_argument('-u','--use_score',dest='score',action='store_true',
                        help = 'Flag to use average difference of scores instead of FEP')
    parser.add_argument('-z','--use_partition_function',dest='partition',action='store_true',
                        help = 'Flag to use ratio of partition functions: -kT log(Zmut/Zwt) instead of FEP')
    parser.add_argument('-s','--show_canvas',dest='canvas',action = 'store_true',
                        help = 'Show canvas figure after execution (not shown by default)')
    parser.add_argument('-o','--output_file',dest='out',action = 'store', default='output',
                        help = 'Output rootname for files (default is output)')
    parser.add_argument('-f','--format_image',dest='format',action = 'store',default='png',
                        help = 'Format of the image in the output file: png(default), pdf, ps, eps and svg ')
    parser.add_argument('-k','--factor_lambda',dest='factor_lambda',action = 'store',default='1.0',
                        help = 'Factor factor_lambda to multiply the energies (default is 1) ')
    parser.add_argument('-v','--verbose',dest='verbose',action = 'store_true',
                        help = 'Verbose execution ')
    parser.add_argument('-c','--clean_null',dest='clean_null',action = 'store',default=None,
                        help = 'Clean ddG values near 0.0 because they are caused by errors (default is None, otherwise select a threshold near 0, i.e. 0.1)')




    options=parser.parse_args()

    return options

  

def main():

 options   = parse_user_arguments()

 if options.ddg is None or options.fep is None  :
    sys.stderr.write("Missing arguments, Check help\n")
    exit(0)

 if options.verbose: sys.stdout.write("Parsing %s\n"%(options.ddg))
 ddg_dict={}
 if not fileExist(options.ddg):
    sys.stderr.write("File %s not found\n"%options.ddg)
 fa=open(options.ddg,"r")
 for line in fa:
    if line.startswith("#"): continue
    form,ddg=line.split()
    if options.clean_null is not None:
      if ddg < -float(options.clean_null) or ddg > float(options.clean_null):
        ddg_dict.setdefault(form,float(ddg))
    else:
      ddg_dict.setdefault(form,float(ddg))
 fa.close()


 if options.verbose: sys.stdout.write("Parsing %s\n"%(options.fep))
 fep_dict={}
 if not fileExist(options.fep):
    sys.stderr.write("File %s not found\n"%options.fep)
 fb=open(options.fep,"r")
 for line in fb:
    data=line.split()
    if data[0].startswith("#"): continue
    data=line.split()
    form=data[-1]
    if options.score:
       if data[-2] is not "nan" and  data[-2] is not "inf": fep=float(data[6])*float(options.factor_lambda)
    elif options.partition:
       if data[4] is not "nan" and  data[4] is not "inf": fep=float(data[4])*float(options.factor_lambda)
    else:
       if data[3] is not "nan" and  data[3] is not "inf": fep=float(data[3])*float(options.factor_lambda)
    if options.outlier_threshold is not None:
     if options.clean_null is not None:
      if fep < -float(options.clean_null) or fep > float(options.clean_null):
        if fep < float(options.outlier_threshold) and fep > -float(options.outlier_threshold): fep_dict.setdefault(form,fep)
     else:
      if fep < float(options.outlier_threshold) and fep > -float(options.outlier_threshold): fep_dict.setdefault(form,fep)
    else:
     if options.clean_null is not None:
      if fep < -float(options.clean_null) or fep > float(options.clean_null):
        fep_dict.setdefault(form,fep)
     else:
      fep_dict.setdefault(form,fep)
      
 fb.close()


 ddg_list=[]
 fep_list=[]
 for form,ene in ddg_dict.iteritems():
   if fep_dict.has_key(form):
      ddg_list.append(ene)
      fep_list.append(fep_dict.get(form))
      


 
 ddg_real=np.array(ddg_list)
 ddg_pred=np.array(fep_list)


 (pearson,p_value)=st.pearsonr(ddg_real,ddg_pred)

 A = np.vstack([ddg_pred, np.ones(len(ddg_pred))]).T
 m, c = np.linalg.lstsq(A, ddg_real,rcond=-1 )[0]

 output= options.out+".dat"
 graph = options.out+"."+options.format

 if float(options.outlier_percentil) > 0:
   (ddg_list_reduced,fep_list_reduced)=reduce_outliers_lstsq(ddg_list,fep_list,m,c,float(options.outlier_percentil))
   ddg_real=np.array(ddg_list_reduced)
   ddg_pred=np.array(fep_list_reduced)
   (pearson,p_value)=st.pearsonr(ddg_real,ddg_pred)
   A = np.vstack([ddg_pred, np.ones(len(ddg_pred))]).T
   m, c = np.linalg.lstsq(A, ddg_real,rcond=-1 )[0]


 fo=open(output,"w")
 fo.write("Pearson Correlation:\t%f\n"%pearson)
 fo.write("Two-tailed p-value: \t%f\n"%p_value)
 fo.write("Slope of fitting:   \t%f\n"%m)
 fo.write("Cut Y-axis:         \t%f\n"%c)
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



if __name__=="__main__":
 main()

 

