import sys, os, argparse, gzip

from SBI.beans           import Path
from SBI.databases import DrugBanklink

from archdb.src import Source, Drug
from archdb.src import start_transaction, end_transaction

# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'DrugBank2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database", 
                        action = "store",       required = True,  help = "Destination Directory for drugBank database.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile", 
                        action = "store",       required = True,  help = "Output SQL file")
    parser.add_argument('-k', '--skipdownload', dest     = "skip_download", 
                        action = "store_true",  default  = False,  help = "Skips downloading the database.")
    parser.add_argument('-v', '--verbose',      dest     = "verbose", 
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options

def DrugBank2SQL(database, sqlfile, skip_download, verbose):
    drugbank_connect = DrugBanklink(local = database)
    newsource        = None
    if not skip_download:
        if verbose: sys.stderr.write("Downloading drugBank database to {0} ...\n".format(database))
        # drugbank_connect.download()
        newsource = Source(name = 'DrugBank', source = drugbank_connect.source)
        if verbose: sys.stderr.write("Download Finished.\n")
    else:
        if verbose: sys.stderr.write("Using previously downloaded database.\n")

    if verbose: sys.stderr.write("Parsing drugBank.\n")
    if verbose: sys.stderr.write("Writing {0} ....\n".format(sqlfile))
    Path.mkdir(os.path.split(os.path.abspath(sqlfile))[0])
    sql_fd = gzip.open(sqlfile, 'wb')
    sql_fd.write(start_transaction())
    sql_fd.write(Drug.preuniprotdeleted())
    if newsource is not None:
        sql_fd.write(newsource.toSQL())

    for drg_line in drugbank_connect.localDrugs:
        newdrg = Drug(inline = drg_line)
        sql_fd.write(newdrg.toSQL())

    sql_fd.write(Drug.afteruniprotdeleted())
    sql_fd.write(end_transaction())
    sql_fd.close()
    if verbose: sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    DrugBank2SQL(options.database, options.sqlfile, options.skip_download, options.verbose)