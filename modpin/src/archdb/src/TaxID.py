from SBI.databases import TaxID

tables = {'taxid'  : 'taxid',
          'taxold' : 'taxid_old'}

class TaxID(TaxID):
    def __init__(self, taxid = None, inline = None):
        super(TaxID, self).__init__(taxid = taxid, inline = inline)

    def toSQL(self):
        if not self.has_old:
            return 'INSERT INTO {0} VALUES ({1.taxid},"{1.name}",{1.parent},"{1.rank}");'.format(tables['taxid'], self)
        else:
            if not self.has_new:
                return 'INSERT INTO {0} (oldid) VALUES ({1.taxid});'.format(tables['taxold'], self)
            else:
                return 'INSERT INTO {0} VALUES ({1.taxid},{1.new});'.format(tables['taxold'], self)