import sys
import os
import inspect
import traceback

class parameters(object):
    """parameters controls some IO minor properties

        To use, include the SBIglobals variable in your script and modify.
    """
    def __init__(self):
        self._verbose   = False      # Promps extra information during different processes performed by the library
        self._debug     = False      # Promps even more information than verbose (use when fixing or modifying the library)
        self._ddebug    = False      # As debug, but for things that are going to repeat itsef a lot (like in Residue or Contact)
        self._overwrite = False      # General set to force overwrite existing files in all the processes.
                                     # Local overwrite (available in different functions) MUST take precedence
        self._warning   = False      # Promps advice if needed
        self._error     = True       # Explains errors
        self._fd        = sys.stderr # Where to print the stuff
    #
    # ATTRIBUTES
    #
    @property
    def verbose(self):          return self._verbose
    @verbose.setter
    def verbose(self, value):   self._verbose = value

    @property
    def debug(self):            return self._debug
    @debug.setter
    def debug(self, value):
        if value:  self._verbose, self._debug = value, value
        else:      self._debug                = value

    @property
    def deepdebug(self):        return self._debug
    @debug.setter
    def deepdebug(self, value):
        if value:  self._verbose, self._debug, self._ddebug = value, value, value
        else:      self._ddebug                             = value

    @property
    def warning(self):          return self.warning
    @warning.setter
    def warning(self, value):   self._warning = value

    @property
    def overwrite(self):        return self._overwrite
    @overwrite.setter
    def overwrite(self, value): self._overwrite = value

    #
    # METHODS
    #
    def output_file(self, debug_file):
        self._fd = open(debug_file,'w')

    def set_file(self):
        execname = ".".join([os.path.splitext(inspect.getfile(inspect.currentframe()))[0], os.getpid(), 'debug'])
        self._fd = open(execname,'w')

    def alert(self, level, source_object, message, error = None, killit = True):
        if self._active_level(level):
            name = getattr(source_object, '__class__', source_object)
            tag  = ''
            if level   == 'error':   tag = '[[E!]]'
            elif level == 'warning': tag = '[[W!]]'
            self._fd.write('{0}[{1}]: {2}\n'.format(tag, name, message))
            if error is not None:
                self._fd.write(traceback.format_exc())
            if level == 'error':
                self._fd.flush()
                if killit: sys.exit()

    def error(self, source_object, message, error = None, killit = True):
        self.alert('error', source_object, message, error, killit)

    def decide_overwrite(self, local_overwrite):
        return self.overwrite if local_overwrite is None else local_overwrite
    #
    # PRIVATE METHODS
    #
    def _active_level(self, level):
        return (level == 'error'                         or (level == 'warning' and self.warning) or \
               (level == 'verbose'   and self.verbose)   or (level == 'debug'   and self.debug)   or \
               (level == 'deepdebug' and self.deepdebug))

SBIglobals = parameters()
