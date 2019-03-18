#!/usr/bin/env python3

import sys

if __name__ == "__main__":
    with open(sys.argv[1]) as handle:
        name = ''
        description = []
        for i, line in enumerate(handle):
            line = line.strip()
            if line.startswith('//'):
                print(name, '; '.join(description), sep='\t')
                name = ''
                description = []
            elif line.startswith('NAME'):
                assert name == '', f'Name for record at line {i} is already set to "{name}"'
                name = line[6:]
            elif line.startswith('DESC'):
                description += [line[6:]]
        if name:
            print(name, '; '.join(description), sep='\t')

