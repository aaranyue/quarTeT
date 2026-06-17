#!/usr/bin/env python3
# Last modified: V1.3.0

import sys
from logging.config import dictConfig
from quartet import assemblymapper, gapfiller, centrominer, teloexplorer

USAGE = '''quarTeT: Telomere-to-telomere Toolkit
version 1.3.0

Usage: quartet <module> <parameters>

Modules:
  AssemblyMapper  | am    Assemble draft genome.
  GapFiller       | gf    Fill gaps in draft genome.
  TeloExplorer    | te    Identify telomeres.
  CentroMiner     | cm    Identify centromere candidates.

Use <module> -h for module usage.
'''

LOGGING_CONFIG = {
    'version': 1,
    'formatters': {
        'default': {
            'format': '[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s',
        },
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'default',
        },
    },
    'root': {
        'level': 'INFO',
        'handlers': ['console'],
    },
}

dictConfig(LOGGING_CONFIG)

def main():
    if len(sys.argv) == 1:
        print(USAGE)
        sys.exit(0)
    
    module = sys.argv[1]
    parameters = sys.argv[2:] if len(sys.argv) > 2 else []
    
    valid_modules = {
        'AssemblyMapper': assemblymapper.main,
        'am': assemblymapper.main,
        'GapFiller': gapfiller.main,
        'gf': gapfiller.main,
        'CentroMiner': centrominer.main,
        'cm': centrominer.main,
        'TeloExplorer': teloexplorer.main,
        'te': teloexplorer.main,
        '-h': None,
        '--help': None,
    }
    
    if module not in valid_modules:
        print('Unexpected module name.')
        print(USAGE)
        sys.exit(1)
    
    if module in ['-h', '--help']:
        print(USAGE)
        sys.exit(0)
    
    valid_modules[module](parameters)

if __name__ == '__main__':
    main()