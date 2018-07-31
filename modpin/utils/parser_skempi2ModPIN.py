import os,sys,re
import argparse
import ConfigParser
import shutil
import string
import numpy

# Add "." to sys.path #
src_path =  os.path.abspath("/home/boliva/PROGRAMS/Interactome/scripts/")
sys.path.append(src_path)

# Read configuration file #
config = ConfigParser.ConfigParser()
config_file = os.path.join(src_path, "config.ini")
config.read(config_file)

# Read SBI path
sbi_path = os.path.join(src_path, config.get("Paths", "sbi_path"))
sys.path.append(sbi_path)
from SBI.structure import *


def fileExist (file):
 if file is not None:
  return os.path.exists(file) and os.path.isfile(file)
 else:
  return False



def printverbose(f,flag,message):
 """Define verbose to print in file 'f', if 'flag'=True, message given as 'message'"""
 if flag: f.write("%s"%(message))


def split(pdb_file,out):
 """
  Split PDB in chains: it creates the PDB and the FastA sequences of each chain and returns the set of pairs (name,sequence)
 """
 if fileExist(pdb_file):
   pdb=PDB(pdb_file)
 else:
  sys.stderr.write("Missing PDB file %s \n"%pdb_file)
  return None
 sequences_split=set()
 for code in pdb.chain_identifiers:
   chain=pdb.get_chain_by_id(code)
   pdb_new=PDB()
   pdb_new.add_chain(chain)
   new_name=out+code+".pdb"
   new_fasta=out+code+".fa"
   name=out.split("/")[-1]+code
   sequence=None
   if chain.chaintype == "N": sequence=chain.nucleotide_sequence()
   if chain.chaintype == "P": sequence=chain.protein_sequence
   #pdb_new.clean
   #pdb_new.write(new_name,force=True)
   #fasta=open(new_fasta,"w")
   #fasta.write(">%s\n%s\n"%(name,sequence))
   #fasta.close()
   if sequence is not None: sequences_split.add((name,sequence))

 return sequences_split
 
def parse_skempi(skempi_file):

   if fileExist(skempi_file):
      fd=open(skempi_file)
   else:
      sys.stderr.write("Missing SKEMPI CSV file %s \n"%skempi_file)
      return None
   ddG_skempi  = {}
   for line in fd:
     if line.startswith('"#'): continue
     data      = line.split('\t')
     try:
      pdb_chain = data[0].lstrip('"').rstrip('"')
      pdb_name  = pdb_chain.split("_")[0]
      chain_a   = pdb_chain.split("_")[1]
      chain_b   = pdb_chain.split("_")[2]
      mutation  = data[2].lstrip('"').rstrip('"')
      Kd_mut    = float(data[7].lstrip('"').rstrip('"'))
      Kd_wt     = float(data[9].lstrip('"').rstrip('"'))
      temp      = float("".join(data[13].lstrip('"').rstrip('"').split()[0][0:2]))
      ddG       =  (8.314/4184.0) * (273.15 + temp) * numpy.log(Kd_mut) - (8.314/4184.0) * (273.15 + temp) * numpy.log(Kd_wt)
      set_a    = set()
      set_b    = set()
      for mutant in mutation.split(","):
           aa_wt    = mutant[0]
           chain    = mutant[1]
           position = mutant[2:-1]
           aa_mut   = mutant[-1]
           if chain == chain_a:
              set_a.add(aa_wt + position + aa_mut)
           if chain == chain_b:
              set_b.add(aa_wt + position + aa_mut)
      form_a="+".join([x for x in set_a])
      form_b="+".join([x for x in set_b])
      pair_wt = pdb_name + "+" + chain_a + "::" + pdb_name + "+" + chain_b
      if len(set_a)>0:
         pair_mut= pdb_name + "+" + chain_a + "_" + form_a + "::"
      else:
         pair_mut= pdb_name + "+" + chain_a  + "::"
      if len(set_b)>0:
         pair_mut= pair_mut + pdb_name + "+" + chain_b + "_" + form_b
      else:
         pair_mut= pair_mut + pdb_name + "+" + chain_b
      ddG_skempi.setdefault(pair_mut,ddG)  
     except Exception as e:
      sys.stderr.write("Error %s. Skip %s %s\n"%(e,pdb_name,mutation))
   return ddG_skempi


def parse_user_arguments(*args, **kwds):

    parser = argparse.ArgumentParser(
        description = "Parse SKEMPI and create MODPIN files",
        epilog      = "@oliva's lab 2018")
    parser.add_argument('-i','--skempi_file',dest='skempi',action = 'store',
                        help = 'SKEMPI Input CSV table')
    parser.add_argument('-p','--pdb_folder',dest='pdb',action = 'store',
                        help = 'PDB directory (i.e. downloaded from SKEMPI)')
    parser.add_argument('-o','--output_file',dest='out',action = 'store', default='output',
                        help = 'Output rootname for files (default is output)')
    parser.add_argument('-m','--max_mutations',dest='maxmut',action='store', default=None,
                        help = 'Maximum number of mutations per chain')
    parser.add_argument('-v','--verbose',dest='verbose',action = 'store_true',
                        help = 'Verbose execution ')


    options=parser.parse_args()

    return options

def main():

 options   = parse_user_arguments()
 if options.verbose: sys.stdout.write("Parsing Skempi\n")
 ddG       = parse_skempi(options.skempi)
 ppi_file  = options.out + ".ppi"
 seq_file  = options.out + ".fa"
 ddG_file  = options.out + ".ddG"
 if options.maxmut is not None:
     maxmut=int(options.maxmut)

 interactions=set()
 fasta={}
 fg=open(ddG_file,"w")
 fg.write("#%39s\t%10s\n"%("Mutant","ddG"))

 for pair,binding in ddG.iteritems():
     if options.verbose: sys.stdout.write("Parse %s\n"%(pair))
     a,b  = pair.split("::")
     wt_a = a.split("_")[0]
     wt_b = b.split("_")[0]
     form_a=None
     form_b=None
     if len(a.split("_"))>1: form_a=a.split("_")[1]
     if len(b.split("_"))>1: form_b=b.split("_")[1]
     pdbA = wt_a.split("+")[0]
     pdbB = wt_b.split("+")[0]
     if len(wt_a.split("+"))>1: chain_a = wt_a.split("+")[1]
     if len(wt_b.split("+"))>1: chain_b = wt_b.split("+")[1]
     if len(chain_a)>1 or len(chain_b)>1:
        if options.verbose: sys.stdout.write("\t-- Skip this interaction\n")
        continue
     if form_a is not None and options.maxmut is not None:
       if len(form_a.split("+"))>maxmut:
        if options.verbose: sys.stdout.write("\t-- Skip this interaction\n")
        continue
     if form_b is not None and options.maxmut is not None:
       if len(form_b.split("+"))>maxmut:
        if options.verbose: sys.stdout.write("\t-- Skip this interaction\n")
        continue
     if pdbA == pdbB:
        pdb_file=options.pdb+"/"+pdbA+".pdb"
        out_file=options.pdb+"/"+pdbA+"+"
        if fileExist(pdb_file):
           sequences=split(pdb_file,out_file)
        else:
           sys.stderr.write("Missing PDB file %s \n"%pdb_file)
     else:
        sys.stderr.write("Interaction in diferent PDB files %s vs %s\n"%(pdbA,pdbB))
     interactions.add(wt_a+"::"+wt_b)
     interactions.add(pair)
     fg.write("%40s\t%10.5f\n"%(pair,binding))
     for (name,seq) in  sequences:
        if name == wt_a:
           fasta.setdefault(name,seq)
           if form_a is not None:
              for mutant in form_a.split("+"):
                  aa_wt    = mutant[0]
                  position = int(mutant[1:-1])
                  aa_mut   = mutant[-1]
                  if seq[position-1] == aa_wt:
                     seq_new=seq[0:position-1] + aa_mut +seq[position:]
                     fasta.setdefault(a,seq_new)
                  else:
                     if options.verbose: sys.stderr.write("Wrong substitution %s in %s\n"%(mutant,a))

        if name == wt_b:
           fasta.setdefault(name,seq)
           if form_b is not None:
              for mutant in form_b.split("+"):
                  aa_wt    = mutant[0]
                  position = int(mutant[1:-1])
                  aa_mut   = mutant[-1]
                  if seq[position-1] == aa_wt:
                     seq_new=seq[0:position-1] + aa_mut +seq[position:]
                     fasta.setdefault(b,seq_new)
                  else:
                     if options.verbose: sys.stderr.write("Wrong substitution %s in %s\n"%(mutant,b))
     
 fg.close()

 fa=open(seq_file,"w")
 for name,seq in fasta.iteritems():
   fa.write(">%s\n%s\n"%(name,seq))
 fa.close()
 fi=open(ppi_file,"w")
 for ppi in interactions:
   a,b=ppi.split("::")
   fi.write("%s\t%s\n"%(a,b))
 fi.close()


if  __name__ == "__main__":
    main()


