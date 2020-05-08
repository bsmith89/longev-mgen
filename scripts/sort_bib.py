#!/usr/bin/env python
"""Sort all of the entries in the files found in the arguments.

Prints the sorted entries to stdout.

Entry recognition is based on the @ARTICLE command syntax, so @PREAMBLE,
@COMMENT, and @STRING are all treated as bibleography entries.
with the string leading up to the first comma used as the key.

TODO: Fix the sorting to deal correctly with non-entry items in the file.

"""


import sys
import itertools
from biblib import bib

def entries(path):
    with open(path) as f:
        db = bib.Parser().parse(f, log_fp=sys.stderr).get_entries()
    for entry in db.values():
        yield entry

def get_key(entry):
    return entry.key

if __name__ == '__main__':
    entries_chain = itertools.chain(*[entries(path) for path in sys.argv[1:]])
    sorted_entries = sorted(entries_chain, key=lambda e: get_key(e))
    for entry in sorted_entries:
        print(entry.to_bib())
