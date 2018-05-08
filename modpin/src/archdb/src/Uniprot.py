from SBI.databases.uniprot import Uniprot

tables = {'uniprot' : 'uniprot',
          'ent2acc' : 'uniprot_entry2accession',
          'uni2go'  : 'uniprot2GO',
          'uni2tax' : 'uniprot2taxid',
          'unienz'  : 'uniprot2enzyme',

          'taxid'   : 'taxid',
          'taxold'  : 'taxid_old',

          'go'      : 'GO',
          'GO_alt'  : 'GO_alternative'}

class Uniprot(Uniprot):
    def __init__(self, entry = None, source = None, inline = None):
        super(Uniprot, self).__init__(entry, source, inline)

    def toSQL(self):
        command = []
        command.append('INSERT INTO {0} VALUES ("{1.entry}","{1.source}");\n'.format(tables['uniprot'],self))
        for a in self.accession:
            command.append('INSERT INTO {0} VALUES ("{1.entry}","{2}");\n'.format(tables['ent2acc'], self, a))
        if self.taxid is not None:
            command.append('SET @taxid = (SELECT COALESCE((SELECT id FROM {0} WHERE id={1}),(SELECT taxid FROM {2} WHERE oldid={1})));\n'.format(tables['taxid'], self.taxid, tables['taxold']))
            command.append('INSERT INTO {0} VALUES ("{1.entry}",@taxid,"is");\n'.format(tables['uni2tax'], self))
        for h in self.hosts:
            command.append('SET @taxid = (SELECT COALESCE((SELECT id FROM {0} WHERE id={1}),(SELECT taxid FROM {2} WHERE oldid={1})));\n'.format(tables['taxid'], h, tables['taxold']))
            command.append('INSERT INTO {0} VALUES ("{1.entry}",@taxid,"host");\n'.format(tables['uni2tax'], self))
        for db in self.databases:
            if db[0] == 'GO':
                command.append('SET @goid = (SELECT COALESCE((SELECT nid FROM {0} WHERE id="{1}"),(SELECT GO_nid FROM {2} WHERE alternative="{1}")));\n'.format(tables['go'], db[1], tables['GO_alt']))
                command.append('INSERT IGNORE INTO {0} VALUES ("{1}",@goid);\n'.format(tables['uni2go'], self.entry))
        return "".join(command)
