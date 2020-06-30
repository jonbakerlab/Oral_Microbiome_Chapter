#!/usr/bin/env python
# encoding: utf-8
# -*- mode: python -*-
from __future__ import division, absolute_import, print_function, unicode_literals  # as python3 as possible

import sys
import errno
import argparse
import logging
from gzip import GzipFile
from bz2 import BZ2File
from enum import Enum
import re


log = logging.getLogger(__name__)


# utility functions
def open_any(filename):
    """Opens uncompressed or compressed files, returning the file object,
    or, if passed a filename of '-', returns STDIN"""
    if filename == '-':
        fh = sys.stdin
    elif filename[-3:] == '.gz':
        fh = GzipFile(filename, 'r')
    elif filename[-4:] == '.bz2':
        fh = BZ2File(filename, 'r')
    else:
        fh = open(filename, 'r', encoding='latin1')

    return fh


class Rec(object):
    """holder for records"""
    def __init__(self, id=None, source=None, bitscore=None, pct_ids=None, evalues=None):
        if pct_ids is None:
            pct_ids = []
        if evalues is None:
            evalues = []
        self.id = id
        self.source = source
        self.bitscore = bitscore
        self.pct_ids = pct_ids
        self.evalues = evalues

    def to_tab(self):
        return '{}\t{}\t{}\t{}\t{}'.format(self.id, self.source, self.bitscore, ','.join((str(i) for i in self.pct_ids)), ','.join((str(i) for i in self.evalues)))

States = Enum('States', 'start in_rec have_id have_source have_bitscore in_hits')


class App(object):
    """
        Main script wrapper: handles argument parsing and file opening
    """

    def __init__(self, argv=None):
        if not argv:
            argv = sys.argv
        # parse arguments
        self.exe_name = argv[0]
        self.options = None
        self.parse_options(argv[1:])

        # set up the logger
        log_level = logging.WARNING
        if self.options.verbose:
            log_level = logging.DEBUG
        log.level = log_level

        log.debug('%s %s', self.exe_name, ' '.join(argv[1:]))

    def parse_options(self, args):
        p = argparse.ArgumentParser(description='<description of the program function>')
        p.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose logging')
        # TODO add options here
        p.add_argument('files', metavar='file', nargs='*', help='File to read. If none, read from STDIN')
        self.options = p.parse_args(args)

    re_id = re.compile(r'^1. (\S+)')
    re_source = re.compile(r'^Source: (.+)')
    re_bitscore = re.compile(r'^Cumulative Blast bit score: (\d+)')
    re_start_hits = re.compile(r'^Table of Blast hits')
    re_pct_id_line = re.compile(r'^\S+\t\S+\t(\d+)\t\d+\t\d+(?:\.\d+)?(?:e-?\d+)?\t(\d+(?:\.\d+)?(?:e-?\d+)?)')

    def main(self):

        files = self.options.files
        if not files:
            files = '-'

        for f in files:
            with open_any(f) as fh:
                # TODO process files
                state = States.start
                cur = Rec()

                for line in (l.rstrip('\n') for l in fh):
                    if state == States.start and line == '>>':  # first record
                        state = States.in_rec
                        if cur.id:
                            print(cur.to_tab())
                            cur = Rec()
                    elif state == States.in_rec:
                        m = App.re_id.match(line)
                        if m:
                            state = States.have_id
                            cur.id = m.group(1)
                    elif state == States.have_id:
                        m = App.re_source.match(line)
                        if m:
                            state = States.have_source
                            cur.source = m.group(1)
                    elif state == States.have_source:
                        m = App.re_bitscore.match(line)
                        if m:
                            state = States.have_bitscore
                            cur.bitscore = m.group(1)
                    elif state == States.have_bitscore:
                        m = App.re_start_hits.match(line)
                        if m:
                            state = States.in_hits
                    elif state == States.in_hits:
                        if line == '>>':
                            state = States.in_rec
                            print(cur.to_tab())
                            cur = Rec()
                        else:
                            m = App.re_pct_id_line.match(line)
                            if m:
                                cur.pct_ids.append(int(m.group(1)))
                                cur.evalues.append(float(m.group(2)))
                # last record
                if cur.id is not None:
                    print(cur.to_tab())

                    
# -----------------------------------------
if __name__ == '__main__':
    logging.basicConfig()
    app = App()
    # trap Ctrl-C and sigpipe to block tracebacks
    try:
        retval = app.main()
    except KeyboardInterrupt:
        sys.exit(1)
    except IOError as e:
        if e.errno == errno.EPIPE:
            sys.exit(1)
        else:
            raise
    else:
        sys.exit(retval)