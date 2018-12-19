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
import functions
from SeqIO import *
from functions import *
from SBI.external.blast import blast_parser
from SBI.structure.contacts import Complex
from SBI.sequence import Sequence
from SBI.structure import PDB

class ModelException(Exception):
    pass


def main():

    options    = parse_user_arguments()
    verbose    = options.show
    dummy_dir  = options.dummy_dir
    if dummy_dir.startswith("/"):
     cluster_dummy_dir=dummy_dir
    else:
     cluster_dummy_dir=os.path.abspath(dummy_dir)
    outdir     = options.outdir
    parallel   = options.parallel
    query_list = options.query_list
    nmodels    = int(options.nmodels)
    optimize   = options.optimize
    did        = options.did
    fasta      = options.sequence_file
    analysis   = options.analysis
    label      = options.label
    cont       = options.cont
    hydrogens  = options.hbplus
    renumerate = options.renumerate
    zrank      = options.zrank
    foldx      = options.foldx
    rosetta    = options.rosetta
    ssp        = options.ssp

    if label is None:
       label=query_list.split("/")[-1]

    if not os.path.exists(outdir): os.makedirs(outdir)
    make_subdirs(outdir, subdirs =['models'])
    modeldir=os.path.join(outdir,'models')

    flags = " -skip -n %d "%(nmodels)
    energy = " "
    if optimize:  flags = flags + " -opt "
    if did:       flags = flags + " -3did "
    if verbose:   flags = flags + " -v "
    if hydrogens: flags = flags + " --hydrogens "
    if renumerate:flags = flags + " --renumerate "
    if zrank:     energy = energy + " --zrank"
    if foldx:     energy = energy + " --foldx"
    if ssp:       energy = energy + " --split-potentials"
    if rosetta:   energy = energy + " --rosetta"   
    translate={}
    if fileExist(fasta):
     for protein in FASTA_iterator(fasta):
        name = protein.get_identifier()
        if re.search('[a-z][|A-Z0-9][|A-Z0-9]', name):
           accession=name.split("|")[1] 
           entry    =name.split("|")[2] 
        else:
           accession=name
           entry=name
        translate.setdefault(entry.split()[0],accession.split()[0])
    else:
        sys.stderr.write('EXIT: Missing FASTA sequence list file %s\n' %(fasta))
        exit(0)


    fd=open(query_list,"r")
    pair={}
    for line in fd:
       if line.startswith("#"): continue
       try:
        p,q=line.strip().split()
       except:
        continue
       try:
        a=translate[p]
        b=translate[q]
       except:
        a=p
        b=q
       x=a.split("_")[0]
       y=b.split("_")[0]
       if pair.has_key((x,y)) or pair.has_key((y,x)):
        if pair.has_key((x,y)):
         pair.setdefault((x,y),set()).add((a,b))
        else:
         pair.setdefault((y,x),set()).add((b,a))
       else:
         pair.setdefault((x,y),set()).add((a,b))
         pair.setdefault((x,y),set()).add((x,y))
    fd.close()
  
    
    for wt,forms in pair.iteritems():
     output=open(os.path.join(outdir,label+"_"+wt[0]+"_"+wt[1]+".ppi"),"w")
     for p,q in forms:
      output.write("%s\t%s\n"%(p,q))
     output.close()


    src_path    = config.get('Paths','modppi_path')   
    scripts_path= os.path.join(src_path,'scripts')
    python_path = config.get('Paths', 'python_path')

    for wt,forms in pair.iteritems():
      ppi=os.path.join(outdir,label+"_"+wt[0]+"_"+wt[1]+".ppi")
      if parallel:
       if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
       else: cluster_queue=config.get("Cluster", "cluster_queue")
       functions.submit_command_to_queue("%s %s %s -seq %s -ppi %s -o %s -d %s" % ( os.path.join(python_path, "python"),  os.path.join(scripts_path,'modppi.py'), flags, os.path.abspath(fasta), os.path.abspath(ppi), os.path.abspath(modeldir), cluster_dummy_dir),cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),cluster_dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
      else:
       os.system("%s %s %s -seq %s -ppi %s -o %s -d %s" % ( os.path.join(python_path, "python"),  os.path.join(scripts_path,'modppi.py'), flags, os.path.abspath(fasta), os.path.abspath(ppi), os.path.abspath(modeldir),dummy_dir))
 
    if parallel and analysis and not cont:
      sys.stderr.write("Wait until all submissions have finished, then run again with flag '--continue'\n")
      exit(0)  
      
    if analysis:
     pair_model={}
     interactions_done=open(os.path.join(modeldir,'interactions_done.list'),"r")
     for line in interactions_done:
       word=line.strip().split()
       has_model=word[2]
       if has_model=="FAILED":continue
       a,b,has_model,path_list=line.strip().split()
       if pair_model.has_key((a,b)):
          pair_model.pop((a,b))
          pair_model.setdefault((a,b),path_list)
       else:
          pair_model.setdefault((a,b),path_list)
     interactions_done.close()
     for wt,forms in pair.iteritems():
      done_list =[(p,q) for p,q in forms if pair_model.has_key((p,q)) or  pair_model.has_key((q,p)) ]
      study_list=open(os.path.join(outdir,label+"_"+wt[0]+"_"+wt[1]+".list"),"w")
      for (p,q) in done_list:
        if pair_model.has_key((p,q)):
          study_list.write('%s\t%s\tDONE\t%s\n'%(p,q,pair_model[(p,q)]))
        else:
          study_list.write('%s\t%s\tDONE\t%s\n'%(q,p,pair_model[(q,p)]))
      study_list.close()
      if parallel:
       if  config.get("Cluster", "cluster_queue") == "None": cluster_queue=None
       else: cluster_queue=config.get("Cluster", "cluster_queue")
       functions.submit_command_to_queue("%s %s -l %s -ppi %s --hydrogens --renumerate  -boxplot -v -d %s -o %s -seq %s %s " % ( os.path.join(python_path, "python"),  os.path.join(scripts_path,'analysis.py'),   label+"_"+wt[0]+"_"+wt[1], os.path.abspath(os.path.join(outdir,label+"_"+wt[0]+"_"+wt[1]+".list")), cluster_dummy_dir, os.path.abspath(outdir),os.path.abspath(fasta),energy ), cluster_queue, int(config.get("Cluster", "max_jobs_in_queue")),os.path.join(scripts_path,config.get("Cluster","command_queue")),cluster_dummy_dir,config.get("Cluster","cluster_submit"),config.get("Cluster","cluster_qstat"))
      else:  
       os.system("%s %s  -l %s -ppi %s --hydrogens --renumerate  -boxplot -v -d %s -o %s -seq %s %s " % ( os.path.join(python_path, "python"),  os.path.join(scripts_path,'analysis.py'),   label+"_"+wt[0]+"_"+wt[1], os.path.abspath(os.path.join(outdir,label+"_"+wt[0]+"_"+wt[1]+".list")), dummy_dir,os.path.abspath(outdir),os.path.abspath(fasta),energy)) 
 

      
    sys.stdout.write("All requested submissions are done.\n")
      
        
       
     
def parse_user_arguments(*args, **kwds):
    parser = argparse.ArgumentParser(
        description = 'Automatic modelling of protein-protein interactions',
        epilog      = '@Oliva\'s lab 2016')
    parser.add_argument('-i', '--query_list', dest = 'query_list', action = 'store',
                        help = 'Input file with a list of pairs of proteins to group and test')
    parser.add_argument('-l', '--label_name', dest = 'label', action = 'store', default = None,
                        help = 'Label to store the analyses (default is the name of the input query)')
    parser.add_argument('-seq', '--sequences', dest = 'sequence_file', action = 'store',
                        help = 'Input file with a list of sequences in FASTA format')
    parser.add_argument('-o', '--output_directory', dest = 'outdir', action = 'store', default = 'ModPPI',
                        help = 'Output directory (default is ModPPI)')
    parser.add_argument('-n', '--number_of_models', dest = 'nmodels', action = 'store', default = 1, type=int,
                        help = 'Number of models for each template (default is 1)')
    parser.add_argument('-d', '--dummy_dir', dest = 'dummy_dir', action = 'store', default = '/tmp/modppi_dummy',
                        help = 'Specifies the location of the dummy folder (default is /tmp/modppi_dummy)')
    parser.add_argument('-opt', '--optimize', dest = 'optimize', action = 'store_true',
                        help = 'Flag to allow model optimization (default is False)')
    parser.add_argument('-3did', '--use_domain_interactions', dest = 'did', action = 'store_true',
                        help = 'Flag to include domain-domain interactions from 3DiD (default is False)')
    parser.add_argument('-a','--analysis', dest = 'analysis', action = 'store_true',
                        help = 'Flag to include the analysis of ddG in the runs (default is False)')
    parser.add_argument('-c','--continue', dest = 'cont', action = 'store_true',
                        help = 'Flag to continue with analysis after all models are done in a cluster (default is False)')
    parser.add_argument('-v', '--verbose', dest = 'show', action = 'store_true',
                        help = 'Flag for verbose mode (default is False)')
    parser.add_argument("-j","--parallel", default=False, action="store_true", dest="parallel", 
                        help="Submit JOBS to Queues in parallel (default = False)")
    parser.add_argument('-hydro','--hydrogens', dest = 'hbplus', action = 'store_true',
                        help = 'Flag to include hydrogens in modelling (default is False, and always True for analyses)')
    parser.add_argument('-r','--renumerate', dest = 'renumerate' , action = 'store_true',
                        help = 'Flag to renumber the sequences as in the original FastA (default is False, and always True for analyses)')
    parser.add_argument('-zrank', '--zrank', dest = 'zrank', action = 'store_true',
                        help = 'Flag to calculate energies with ZRANK (default is False)')
    parser.add_argument('-foldx', '--foldx', dest = 'foldx', action = 'store_true',
                        help = 'Flag to calculate energies with FoldX (default is False)')
    parser.add_argument('-rosetta', '--rosetta', dest = 'rosetta', action = 'store_true',
                        help = 'Flag to calculate energies with ROSETTA (default is False)')
    parser.add_argument('-ssp', '--split-potentials', dest = 'ssp', action = 'store_true',
                        help = 'Flag to calculate energies with Split-Statistic Potentials (default is False unless no other method is used)')
    options = parser.parse_args()
    return options


if __name__ == '__main__':
    main()
