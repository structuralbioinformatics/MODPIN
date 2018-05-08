# -*-
#
# @author: jaumebonet
# @email:  jaume.bonet@gmail.com
# @url:    jaumebonet.github.io
#
# @date:   2015-05-08 16:15:32
# @lab:    LPDI/EPFL
#
# @last modified by:   jaumebonet
# @last modified time: 2015-05-10 15:07:37
#
# -*-
import re

from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


def basic_parser(*args, **kargs):
    '''
    Include those options that most likely are going to be shared
    by all (or most of) the modpy scripts.

    Advantages:
        1) Write less
        2) Homogenization of parameters

    @return: argparse.ArgumentParser object
    '''
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('--pir', dest='alignment', type=str, action='store',
                        help="PIR formated alignment", metavar="PIR_FILE")
    parser.add_argument("--out", dest="out", type=str, action="store",
                        default=None, metavar="OUT_PREFIX",
                        help="Prefix for the log and error file")
    parser.add_argument("-v", dest="verbose", action="store_true",
                        default=False, help="Verbose Mode")
    return parser


def identify_pir(filename):
    '''
    Identify known structures and query sequence in the PIR alignment.

    @return: {'str': [str(), ...], 'seq': str()}
    '''
    _idr = re.compile('^>\w+\;(\S+)')
    data = {'str': [], 'seq': None}

    with open(filename) as fd:
        nextline = False
        mcontent = None
        for line in fd:
            m = _idr.match(line)
            if m:
                nextline = True
                mcontent = m.group(1).strip()
                continue
            if nextline and line.startswith('structureX'):
                    data['str'].append(mcontent)
                    nextline, mcontent = False, None
            elif nextline and line.startswith('sequence'):
                    data['seq'] = mcontent
                    nextline, mcontent = False, None
    return data
