# -*-
#
# @author: jaumebonet
# @email:  jaume.bonet@gmail.com
# @url:    jaumebonet.github.io
#
# @date:   2015-05-08 16:16:50
# @lab:    LPDI/EPFL
#
# @last modified by:   jaumebonet
# @last modified time: 2015-05-11 13:33:39
#
# -*-
'''
Given a PDB file and a minimal set of parameters generates models
through homology modeling with MODELLER.

[!] The script does not contain any input check.
[!] The simple_model function redirects STDOUT and STDERR to file logs.
    The redirection is terminated before the function ends.

'''
import sys
import os
import warnings
from argparse import ArgumentParser

from modpyHelper import basic_parser
from modpyHelper import identify_pir

from modeller import *
from modeller.automodel import *


def set_options(*args, **kargs):
    '''
    Set the specific options for this script

    @return: parsed ArgumentParser
    '''
    parser = ArgumentParser(parents = [basic_parser()],
                            conflict_handler='resolve')

    parser.add_argument("--models", dest="numMod",     action="store",
                        type=int,   metavar="MOD_NUM", default=5,
                        help="Number of models to be done (min: 1) [Default: 5]")

    parser.add_argument("--optimize",  dest="optimize", action="store_true",
                        default=False, help="Activate optimization")
    return parser.parse_args()


def simple_model(alignment, instances,
                 output=None, optimize=False, verbose=False):

    contents = identify_pir(alignment)

    # We redirect STDOUT and STDERR to output files
    if output is None: output = contents['seq']
    warnings.warn('STDOUT & STDERR will be redirected to log files.')
    sys.stdout = open(output + '.log', 'w')
    sys.stderr = open(output + '.err', 'w')

    if verbose: log.verbose()  # Commands MODELLER to display all log output

    env = environ()  # Initializes the 'environment' for this modeling run

    env.io.hetatm = True  # Allow the presence of hetatoms

    a = automodel(env,  # Loading the environment
                  alnfile=alignment,         # Assigning the PIR alignment
                  knowns=contents['str'],    # Listing the known structures
                  sequence=contents['seq'],  # Identify the Query Sequence
                  assess_methods=(assess.DOPE,
                                  assess.GA341,
                                  assess.DOPEHR))  # Energy evaluation methods

    # Setting Starting and Ending Model number
    a.starting_model, a.ending_model = 1, instances

    if optimize:
        # Very thorough VTFM optimization:
        a.library_schedule   = autosched.slow
        a.max_var_iterations = 300

        # Thorough MD optimization:
        a.md_level = refine.slow

        # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
        a.repeat_optimization = 2
        a.max_molpdf          = 1e6

    a.initial_malign3d = True

    a.make()  # Create the Model

    sys.stdout = sys.__stdout__  # restore stdout back to normal
    sys.stderr = sys.__stderr__  # restore stdout back to normal
    warnings.warn('STDOUT & STDERR have been restored.')

    # We delete the error file if no error has occurred
    ERRsize = os.path.getsize(output + '.err')
    if int(ERRsize) is 0:
        os.remove(output + '.err')

if __name__ == '__main__':
    options = set_options()
    simple_model(options.alignment, options.numMod, options.out,
                 options.optimize, options.verbose)
