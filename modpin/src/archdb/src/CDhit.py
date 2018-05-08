from SBI.external.CDhit import CDhitList

tables = {'chain' : 'chain',
          'homo'  : 'chain_homology'}

class CDhit(CDhitList):
    def toSQL(self):
        command = ''
        for hit in self.clusters:
            mpdb, mchain = hit.master.name.split('_')[0], hit.master.name.split('_')[1]
            command += "SET @master = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], mpdb,mchain)
            command += "UPDATE {0} SET clustered = 'M' WHERE nid=@master;\n".format(tables['chain'])
            for s in hit.sequences:
                seq = hit.sequences[s]
                spdb, schain = seq.name.split('_')[0], seq.name.split('_')[1]
                command += "SET @homo = (SELECT nid FROM {0} WHERE pdb='{1}' AND chain='{2}');\n".format(tables['chain'], spdb,schain)
                command += "UPDATE {0} SET clustered = 'H' WHERE nid=@homo;\n".format(tables['chain'])
                command += "INSERT INTO {0} VALUES (@master, @homo, {1.length}, {2.length}, {2.homology});\n".format(tables['homo'], hit.master, seq)
        return command
