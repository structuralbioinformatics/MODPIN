import sys
import os
import re


if len(sys.argv)>2:
  fd=open(sys.argv[1],"r")
  outdir=sys.argv[2]
elif len(sys.argv)==2:
  fd=open(sys.argv[1],"r")
  outdir="./"
else:
  print "Run as parse_list.py 'input_list' [FOLDER NAME]"
  print "Missing input list"
  exit(0)

pair={}
for line in fd:
   a,b=line.strip().split()
   x=a.split("_")[0]
   y=b.split("_")[0]
   if pair.has_key((x,y)) or pair.has_key((y,x)):
    if pair.has_key((x,y)):
     pair.setdefault((x,y),set()).add((a,b))
    else:
     pair.setdefault((y,x),set()).add((b,a))
   else:
     pair.setdefault((x,y),set()).add((a,b))
fd.close()

if not os.path.exists(outdir): os.makedirs(outdir)
  
for wt,forms in pair.iteritems():
   output=open(os.path.join(outdir,wt[0]+"_"+wt[1]+".ppi"),"w")
   for p,q in forms:
    output.write("%s\t%s\n"%(p,q))
   
  

