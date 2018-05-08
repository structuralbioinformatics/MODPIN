from datetime import date

tables = {'sources':'source_databases'}

class Source(object):
    def __init__(self, name, source):
        self.name   = name
        self.source = source
        self.date   = date.today()

    @property
    def date2string(self):
        return '{0}_{1}_{2}'.format(self.date.strftime('%d'), self.date.strftime('%b'), self.date.strftime('%y'))

    @property
    def date2SQL(self):
        return '{0}-{1}-{2}'.format(self.date.strftime('%Y'), self.date.strftime('%m'), self.date.strftime('%d'))

    def toSQL(self):
        command  = ""
        command += "INSERT INTO {0} VALUES ('{1.name}','{1.source}','{1.date2SQL}')\n".format(tables['sources'], self)
        command += "ON DUPLICATE KEY UPDATE "
        command += "{0}.date='{1.date}';\n".format(tables['sources'], self)
        return command