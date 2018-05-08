import re

tables = {'scop'     : 'scop',
          'main'     : 'PDB',
          'scop_d'   : 'scop_description',
          'scop_pdb' : 'chain2scop',
          'chain'    : 'chain',
          'oldPDB'   : 'oldPDB'}

class SCOP(object):
    accepted = ['cl','cf','sf','fa','dm']
    def __init__(self):
        self._descriptions = []
        self._relations    = []
        self._norepel      = set()
        self._maps         = []

    def add_description(self, line):
        line = line.strip().split()
        # if line[1] in self.accepted and line[2] != 'unassigned-sccs':
        if line[1] in self.accepted:
            if line[1] == 'cl': line[-1] = ''
            self._descriptions.append([line[0], line[1], line[2], " ".join(line[4:]).replace("'","\\'")])

    def add_relation(self, line):
        line = line.strip().split()
        # if line[0] == 'unassigned-sccs':return

        data = tuple([tuple(x.split('=')) for x in line[5].split(',')][0:5])

        if not data[4][1] in self._norepel:
            self._relations.append(data)
            self._norepel.add(data[4][1])

        info = []
        if ',' in line[2]:
            multi = line[2].split(',')
            for m in multi:
                info.append(list(line))
                info[-1][2] = m
        else:
            info.append(line)
        for l in info:
            assignation = [l[1],
                           'A' if l[2] == '-' else l[2] if not ':' in l[2] else l[2].split(':')[0],
                           '' if (l[2] == '-' or not ':' in l[2]) else l[2].split(':')[1],
                           data[4][1]]

            self._maps.append(assignation)

    @staticmethod
    def prepdbdeleted():
        command = "INSERT INTO {0} (pdb,date,header,method) VALUES ('0000','000-00-00','FAKE_PROTEIN','X-RAY DIFFRACTION');\n".format(tables['main'])
        chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyz"
        for c in chains:
            command += "INSERT INTO {0} (pdb,chain,start,end) VALUES ('0000','{1}',1,2);\n".format(tables['chain'], c)
        return command
    @staticmethod
    def afterpdbdeleted():
        command = "DELETE FROM {0} WHERE chain IN (SELECT nid FROM {1} WHERE pdb='0000');\n".format(tables['scop_pdb'], tables['chain']) 
        command += "DELETE FROM {0} WHERE pdb = '0000';\n".format(tables['chain'])
        command += "DELETE FROM {0} WHERE pdb = '0000';\n".format(tables['main'])
        return command
    def toSQL(self):
        command = ''
        for d in self._descriptions:
            command += "INSERT INTO {0} (id,type,code,description) VALUES ({1[0]},'{1[1]}','{1[2]}','{1[3]}');\n".format(tables['scop_d'], d)
        for r in self._relations:
            command += "INSERT INTO {0} (domain,family,superfamily,fold,class) VALUES ({1[4][1]},{1[3][1]},{1[2][1]},{1[1][1]},{1[0][1]});\n".format(tables['scop'], r)
        for a in self._maps:
            command += "SET @nid = NULL;\n";
            command += "SELECT nid, start, idxs, end, idxe INTO @nid, @p1, @i1, @p2, @i2 FROM {0} WHERE pdb=COALESCE((SELECT pdb from {3} WHERE oldid = '{1}'),'{1}') AND chain='{2}';\n".format(tables['chain'],a[0].upper(),a[1], tables['oldPDB'])
            command += "SET @nid = (SELECT IF (@nid IS NULL, (SELECT nid FROM {0} WHERE pdb='0000' AND chain='A'), @nid));\n".format(tables['chain'])
            if a[2] != '':
                m = re.search('(-*\d+)(\w*)-{1}(-*\d+)(\w*)',a[2])
                command += "SET @p1 = {0};\n".format(  m.group(1))
                command += "SET @i1 = '{0}';\n".format(' ' if m.group(2) == '' else m.group(2))
                command += "SET @p2 = {0};\n".format(  m.group(3))
                command += "SET @i2 = '{0}';\n".format(' ' if m.group(4) == '' else m.group(4))
            command += "INSERT INTO {0} (domain,start,idxs,end,idxe,chain) VALUES ({1},@p1,@i1,@p2,@i2,@nid);\n".format(tables['scop_pdb'],a[3])
        return command