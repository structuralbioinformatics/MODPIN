from SBI.databases import TM
import re 

tables = {'regions'  :'PDBTM_regions',
          'pdbtm'    :'PDBTM',
          'tm2chain' :'PDBTM2chain',
          'chain'    :'chain',
          'main'     :'PDB',
          'old'      :'oldPDB'}

class TM(TM):
    def toSQL(self):
        command = ''
        command += "INSERT INTO {0} (tmres, type, kwres) VALUES ({1.tmres},'{1.tmtype}',{2});\n".format(tables['pdbtm'], self,0 if self.kwres=='no' else 1)
        command += "SET @tm = LAST_INSERT_ID();\n";
        for chain in self.chains:
            chachain = chain[0] if chain[0] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyz" else 'A'
            command += "SET @chain = NULL;\n"
            command += "SET @chain = (SELECT nid FROM {0} WHERE pdb=COALESCE((SELECT pdb from {3} WHERE oldid = '{1}'),'{1}') AND chain='{2}');\n".format(tables['chain'], self.pdb,chachain, tables['old'])
            command += "SET @chain = (SELECT IF (@chain IS NULL, (SELECT nid FROM {0} WHERE pdb='0000' AND chain='A'), @chain));\n".format(tables['chain'])
            for section in self.chains[chain]:
                ms, me = re.search('(\-*\d+)(\w*)',section[0]), re.search('(\-*\d+)(\w*)',section[1])
                start, idxs, end, idxe = ms.group(1), ms.group(2) if ms.group(2) != '' else ' ', me.group(1), me.group(2) if me.group(2) != '' else ' '
                region = self._switch_region(str(section[-1]))
                command += "INSERT INTO {0} (PDBTM,chain,start,idxs,end,idxe,region) VALUES ".format(tables['tm2chain'])
                command += "(@tm,@chain,{0},'{1}',{2},'{3}','{4}');\n".format(start, idxs, end, idxe, region)
        return command

    @staticmethod
    def prepdbdeleted():
        command = "INSERT INTO {0} (pdb,date,header,method) VALUES ('0000','000-00-00','FAKE_PROTEIN','X-RAY DIFFRACTION');\n".format(tables['main'])
        chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyz"
        for c in chains:
            command += "INSERT INTO {0} (pdb,chain,start,end) VALUES ('0000','{1}',1,2);\n".format(tables['chain'], c)
        return command
    @staticmethod
    def afterpdbdeleted():
        command = "DELETE FROM {0} WHERE chain IN (SELECT nid FROM {1} WHERE pdb='0000');\n".format(tables['tm2chain'], tables['chain']) 
        command += "DELETE FROM {0} WHERE pdb = '0000';\n".format(tables['chain'])
        command += "DELETE FROM {0} WHERE pdb = '0000';\n".format(tables['main'])
        return command
    def _switch_region(self, region):
        if (region != '1' and region != '2') or self.side == '': return region
        oposites = {'Inside':'Outside','Outside':'Inside'}
        translate = {'Inside':'i','Outside':'o'}
        if region == '1':   return translate[self.side]
        else:               return translate[oposites[self.side]]

    @staticmethod
    def regions2SQL():
        command = ''
        for r in TM.section_types:
            command += "INSERT INTO {0} (id, definition) VALUES ('{1}','{2}');\n".format(tables['regions'], r, TM.section_types[r])
        return command