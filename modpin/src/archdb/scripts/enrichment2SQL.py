import sys, os, argparse, gzip

from SBI.beans           import Path

from archdb.src import start_transaction, end_transaction

# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'enrichment2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',          dest     = "inputfile", 
                        action = "store",         required = True,  help = "File of enrichments")
    parser.add_argument('-c', '--classification', dest     = "clust", choices = ['class','subclass'],
                        action = "store",         required = True,  help = "Class or Subclass")
    parser.add_argument('-r', '--relation',       dest     = "rel",   choices = ['scop','go','enzyme','drugBank'],
                        action = "store",         required = True,  help = "Class or Subclass")
    parser.add_argument('-v', '--verbose',        dest     = "verbose", 
                        action = "store_true",    default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    options.__dict__['table'] = "_".join([options.rel,options.clust,'enrichment'])
    return options

def enrichment2SQL(infile, table, rel, verbose):
   
    if verbose: sys.stderr.write("Retrieving data from {0} ...\n".format(infile))
    if verbose: sys.stderr.write("\tTo table {0} ...\n".format(table))
    start_command = "INSERT INTO {0} VALUES ".format(table)
    sys.stdout.write(start_transaction()+"\n")
    fd = open(infile)
    for line in fd:
        d = line.strip().split()
        d[1] = d[1] if rel in ['scop','go'] else "'"+d[1]+"'"
        sys.stdout.write("{0} ({1[0]},{1[1]},'{1[7]}',{1[2]},'{1[8]}','{1[9]}','{1[10]}');\n".format(start_command,d))
    sys.stdout.write(end_transaction())
 
    if verbose: sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    enrichment2SQL(options.inputfile, options.table, options.rel, options.verbose)