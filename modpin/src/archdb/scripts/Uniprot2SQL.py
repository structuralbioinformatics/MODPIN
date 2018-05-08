import sys, os, argparse, gzip

from SBI.beans           import Path
from SBI.databases import Uniprotlink

from archdb.src import Source, Uniprot
from archdb.src import start_transaction, end_transaction

# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'Uniprot2SQL',  formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database',     dest     = "database",
                        action = "store",       required = True,  help = "Destination Directory for Uniprot database.")
    parser.add_argument('-s', '--sql',          dest     = "sqlfile",
                        action = "store",       required = True,  help = "Output SQL file pattern (_)")
    parser.add_argument('-k', '--skipdownload', dest     = "skip_download",
                        action = "store_true",  default  = False,  help = "Skips downloading the database.")
    parser.add_argument('-v', '--verbose',      dest     = "verbose",
                        action = "store_true",  default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options

def Uniprot2SQL(database, sqlfile, skip_download, verbose):
    uniprot_connect = Uniprotlink(local = database)
    newsource       = None

    if not options.skip_download:
        if verbose: sys.stderr.write("Downloading Uniprot database to {0} ...\n".format(database))
        uniprot_connect.download()
        newsource = Source(name = 'uniprot', source = uniprot_connect.source)
        if verbose: sys.stderr.write("Download Finished.\n")
    else:
        if verbose: sys.stderr.write("Using previously downloaded database.\n")

    file_counter  = 1
    file_sequence = 0
    file_sql_name = sqlfile.replace('_','{0:03}')

    if verbose: sys.stderr.write("Parsing Uniprot.\n")
    if verbose: sys.stderr.write("Writing {0} ....\n".format(file_sql_name.format(file_counter)))
    Path.mkdir(os.path.split(os.path.abspath(sqlfile))[0])
    sql_fd = gzip.open(file_sql_name.format(file_counter), 'wb')
    sql_fd.write(start_transaction())
    if newsource is not None:
        sql_fd.write(newsource.toSQL())

    for uni_line in uniprot_connect.localUniprots:
        newuni = Uniprot(inline = uni_line)
        if file_sequence > 500000:
            sql_fd.write(end_transaction())
            sql_fd.close()
            file_sequence = 0
            file_counter += 1
            if verbose: sys.stderr.write("Writing {0} ....\n".format(file_sql_name.format(file_counter)))
            sql_fd = gzip.open(file_sql_name.format(file_counter), 'wb')
            sql_fd.write(start_transaction())

        sql_fd.write(newuni.toSQL())
        file_sequence += 1
    if verbose: sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    Uniprot2SQL(options.database, options.sqlfile, options.skip_download, options.verbose)
