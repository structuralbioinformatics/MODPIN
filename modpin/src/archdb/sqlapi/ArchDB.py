import mysql.connector

class ArchDB(object):
    def __init__(self, dbname=None, dbhost=None, dbuser=None, dbpassword=None):
        self._database  = dbname
        self._host      = dbhost
        self._user      = dbuser
        self._password  = dbpassword

        self._db        = None
        self._cursor    = None

        self._connect()


    #
    # PRIVATE FUNCTIONS
    #

    def _connect(self):

        # Connection to MySQL server 
        if not self.dbuser is None and not self.dbpassword is None:
            self._db = mysql.connector.connect(host= self.dbhost, user=self.dbuser, passwd=self.dbpassword)
        elif not self.dbuser is None:
            self._db = mysql.connector.connect(host= self.dbhost, user=self.dbuser)
        else:
            self._db = mysql.connector.connect(host= self.dbhost)

        # Cursor (from where operations are done and information is retrieved)
        self._cursor = self.db.cursor()

        # Specify Database 
        if( self.dbname is not None ):
            self.db.database=self.dbname