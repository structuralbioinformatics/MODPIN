from SBI.databases import GOterm

tables = {'go'     : 'GO', 
          'GO_alt' : 'GO_alternative',
          'GO_ext' : 'GO_extended',
          'GO_rel' : 'GO_relationship'}

class GOterm(GOterm):
    def __init__(self, id = None, inline = None):
        super(GOterm, self).__init__(id = id, inline = inline)

    """"""
    def toSQL(self):
        data = []
        data.append('INSERT INTO {0} (id, name, namespace, obsolete) VALUES ("{1.id}","{1.name}","{1.namespace_mini}",{1.obsolete_binary});'.format(tables['go'], self))
        data.append('SET @goid = (SELECT nid FROM {0} WHERE id="{1.id}");'.format(tables['go'], self))
        for alt in self.alt_id:
            data.append('INSERT INTO {0} VALUES (@goid, "{1}");'.format(tables['GO_alt'], alt))
        return "\n".join(data)

    def relations2SQL(self):
        data = []
        data.append('SET @goid1 = (SELECT nid FROM {0} WHERE id="{1.id}");'.format(tables['go'], self))
        for rel in self.relations:
            data.append('SET @goid2 = (SELECT nid FROM {0} WHERE id="{1}");'.format(tables['go'], rel[1]))
            data.append('INSERT INTO {0} VALUES (@goid1, @goid2, "{1}");'.format(tables['GO_rel'], rel[0]))
        return "\n".join(data)

    def parents2SQL(self):
        data = []
        data.append('SET @goid1 = (SELECT nid FROM {0} WHERE id="{1.id}");'.format(tables['go'], self))
        for rel in self.parents:
            data.append('SET @goid2 = (SELECT nid FROM {0} WHERE id="{1}");'.format(tables['go'], rel))
            data.append('INSERT INTO {0} VALUES (@goid2, @goid1);'.format(tables['GO_ext']))
        return "\n".join(data)