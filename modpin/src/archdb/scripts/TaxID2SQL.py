import sys, os, argparse, gzip

from SBI.beans           import Path
from SBI.databases import TaxIDlink

from archdb.src import Source, TaxID
from archdb.src import start_transaction, end_transaction

# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'TaxID2SQL', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database", 
                        action = "store",       required = True,  help = "Destination Directory for TaxID database.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile", 
                        action = "store",       required = True,  help = "Output SQL file")
    parser.add_argument('-k', '--skipdownload', dest     = "skip_download", 
                        action = "store_true",  default  = False,  help = "Skips downloading the database.")
    parser.add_argument('-v', '--verbose',      dest     = "verbose", 
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options

def TaxID2SQL(database, sqlfile, skip_download, verbose):
    taxid_connect = TaxIDlink(local = database)
    newsource        = None
    if not skip_download:
        if verbose: sys.stderr.write("Downloading TaxID database to {0} ...\n".format(database))
        taxid_connect.download()
        newsource = Source(name = 'taxid', source = taxid_connect.source)
        if verbose: sys.stderr.write("Download Finished.\n")
    else:
        if verbose: sys.stderr.write("Using previously downloaded database.\n")

    has_new = []
    if verbose: sys.stderr.write("Parsing TaxID.\n")
    if verbose: sys.stderr.write("Writing {0} ....\n".format(sqlfile))
    Path.mkdir(os.path.split(os.path.abspath(sqlfile))[0])
    sql_fd = gzip.open(sqlfile, 'wb')
    sql_fd.write(start_transaction())
    if newsource is not None:
        sql_fd.write(newsource.toSQL())

    for tax_line in taxid_connect.localTaxIDs:
        newtax = TaxID(inline = tax_line)
        if newtax.has_new: has_new.append(newtax.toSQL())
        else:              sql_fd.write(newtax.toSQL() + "\n")

    sql_fd.write("\n".join(has_new) + "\n")
    sql_fd.write(end_transaction())
    sql_fd.close()
    if verbose: sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    TaxID2SQL(options.database, options.sqlfile, options.skip_download, options.verbose)
