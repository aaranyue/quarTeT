#!/usr/bin/env python3

import subprocess
import sys

usage = '''quarTeT: Telomere-to-telomere Toolkit
version 1.2.0

Usage: python3 quartet.py <module> <parameters>

Modules:
AssemblyMapper\t| am\tAssemble draft genome.
GapFiller\t| gf\tFill gaps in draft genome.
TeloExplorer\t| te\tIdentify telomeres.
CentroMiner\t| cm\tIdentify centromere candidates.

Use <module> -h for module usage.'''

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(usage)
        sys.exit(0)
    elif len(sys.argv) > 1:
        module = sys.argv[1]
        if module not in ['AssemblyMapper', 'GapFiller', 'CentroMiner', 'TeloExplorer', 'am', 'gf', 'cm', 'te', '-h', '--help']:
            print('Unexpected parameters. Use -h for help.')
            sys.exit(0)
        parameter = '' if len(sys.argv) == 2 else ' '.join(sys.argv[2:])
    
    if module == 'AssemblyMapper' or module == 'am':
        subprocess.run(f'python3 {sys.path[0]}/quartet_assemblymapper.py {parameter}', shell=True)
    elif module == 'GapFiller' or module == 'gf':
        subprocess.run(f'python3 {sys.path[0]}/quartet_gapfiller.py {parameter}', shell=True)
    elif module == 'CentroMiner' or module == 'cm':
        subprocess.run(f'python3 {sys.path[0]}/quartet_centrominer.py {parameter}', shell=True)
    elif module == 'TeloExplorer' or module == 'te':
        subprocess.run(f'python3 {sys.path[0]}/quartet_teloexplorer.py {parameter}', shell=True)
    elif module == '-h' or module == '--help':
        print(usage)
