import sys
import os
import argparse
import gzip

from SBI.beans     import Path
from SBI.data      import element_dic
from SBI.databases import PDBeChemlink

from archdb.src import Source, PDBeChem, Element
from archdb.src import start_transaction, end_transaction


# User Options
def options(*args, **kwds):
    parser = argparse.ArgumentParser(prog = 'PDBeChem2SQL', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database', dest     = "database",
                        action = "store", required = True, help = "Destination Directory for PDBeChem database.")
    parser.add_argument('-s', '--sql', dest     = "sqlfile",
                        action = "store", required = True, help = "Output SQL file")
    parser.add_argument('-k', '--skipdownload', dest     = "skip_download",
                        action = "store_true", default  = False, help = "Skips downloading the database.")
    parser.add_argument('-v', '--verbose', dest     = "verbose",
                        action = "store_true", default  = False, help = "Verbose Mode.")
    options = parser.parse_args()
    return options


def PDBeChem2SQL(database, sqlfile, skip_download, verbose):
    pdbechem_connect = PDBeChemlink(local = database)
    newsource        = None
    if not skip_download:
        if verbose: sys.stderr.write("Downloading PDBeChem database to {0} ...\n".format(database))
        pdbechem_connect.download()
        newsource = Source(name   = 'PDBeChem', source = pdbechem_connect.source)
        if verbose: sys.stderr.write("Download Finished.\n")
    else:
        if verbose: sys.stderr.write("Using previously downloaded database.\n")

    noparent_chems = []
    parent_chems   = []
    if verbose: sys.stderr.write("Parsing PDBeChem.\n")
    for chem_file in pdbechem_connect.localPDBeChems:
        if verbose: sys.stderr.write("\tReading {0} ....\n".format(chem_file))
        newchem = PDBeChem(chem_file)
        if newchem.parent is None: noparent_chems.append(newchem.toSQL())
        else: parent_chems.append(newchem.toSQL())

    if verbose: sys.stderr.write("Writing {0} ....\n".format(sqlfile))
    Path.mkdir(os.path.split(os.path.abspath(sqlfile))[0])
    sql_fd = gzip.open(sqlfile, 'wb')
    sql_fd.write(start_transaction())
    if newsource is not None:
        sql_fd.write(newsource.toSQL())
    for e in element_dic.values():
        newelement = Element(e.number, e.symbol, e.name)
        sql_fd.write(newelement.toSQL() + "\n")
    sql_fd.write("\n".join(noparent_chems) + "\n")
    sql_fd.write("\n".join(parent_chems) + "\n")
    sql_fd.write(end_transaction())
    sql_fd.close()
    if verbose: sys.stderr.write("End execution.\n")

if __name__ == "__main__":
    options = options()
    PDBeChem2SQL(options.database, options.sqlfile, options.skip_download, options.verbose)
