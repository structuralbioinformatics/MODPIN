import sys, os, argparse, gzip

from SBI.beans           import Path

from archdb.src import CDhit
from archdb.src import start_transaction, end_transaction

# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'CDhit2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database", 
                        action = "store",       required = True,  help = "CDhit output file.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile", 
                        action = "store",       required = True,  help = "Output SQL file")
    parser.add_argument('-v', '--verbose',      dest     = "verbose", 
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options

def CDhit2SQL(database, sqlfile, verbose):
   
    if verbose: sys.stderr.write("Retrieving data from {0} ...\n".format(database))
    cdhit = CDhit(database)

    if verbose: sys.stderr.write("Writing {0} ....\n".format(sqlfile))
    Path.mkdir(os.path.split(os.path.abspath(sqlfile))[0])
    sql_fd = gzip.open(sqlfile, 'wb')
    sql_fd.write(start_transaction())

    sql_fd.write(cdhit.toSQL())

    sql_fd.write(end_transaction())
    sql_fd.close()
    if verbose: sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    CDhit2SQL(options.database, options.sqlfile, options.verbose)

    