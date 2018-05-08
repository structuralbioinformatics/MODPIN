import sys
import os
import re


if len(sys.argv)>3:
  fd=open(sys.argv[1],"r")
  ff=open(sys.argv[2],"r")
  outdir=sys.argv[3]
elif len(sys.argv)==3:
  fd=open(sys.argv[1],"r")
  ff=open(sys.argv[2],"r")
  outdir="output.fa"
else:
  print "Run as parse_seq.py 'input_sequences' 'uniprot_database' [FOLDER NAME]"
  print "Missing input list"
  exit(0)

translate={}
for line in ff:
  if line.startswith(">"):
    accession = line.lstrip(">").strip().split("|")[1]
    entry     = line.lstrip(">").strip().split("|")[2].split()[0]
    translate.setdefault(entry,accession)


ff.close()
fo=open(outdir,"w")
for line in fd:
  if line.startswith(">"):
    name_sequence=line.lstrip(">").strip().split("_")
    if len(name_sequence)>2:
     name="_".join(name_sequence[0:-1])
     extension=name_sequence[-1]
    else:
     name="_".join(name_sequence)
     extension=None
    if translate.has_key(name):
     if extension is None: 
       accession=translate[name]
       entry=name
     else:
       accession=translate[name]+"_"+extension
       entry=name+"_"+extension
     fo.write(">sp|%s|%s\n"%(accession,entry))
    else:
     if extension is None:
      entry=name
     else:
      entry=name+"_"+extension
     fo.write(">%s\n"%(entry))
  else:
   fo.write("%s\n"%(line.strip()))

fd.close()
fo.close()
 

    

