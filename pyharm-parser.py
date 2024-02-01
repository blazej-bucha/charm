# Extracts symbolic constants from the CHarm's header files that are inside
# "PATH" and saves the constants to a Python "FILE" to be used with PyHarm.

import os
import ctypes as ct
import pyclibrary

# List of paths to read the header files from
PATHS = ['./charm/', './src/shs']

# Name of the output Python file to save the symbolic constants to
FILE = './wrap/pyharm/_constants.py'

# Dictionary specifying the names of the CHARM symbolic constants and their new
# names after parsing all header files from "PATH".
REPLACE_DICT = {'ULONG_MAX': '%s' % (1 - ct.c_ulong(-1).value),
                'CHARM_SHC_NMAX_MODEL': 'CHARM_SHC__NMAX_MODEL',
                'CHARM_SHC_NMAX_ERROR': 'CHARM_SHC__NMAX_ERROR'}

VALID_SUFFICES = ['CHARM_', 'GRAD_']

# Create the output file to write the symbolic constants to
fid = open(FILE, 'w')

fid.write('\"\"\"\nGenerated automatically.  Do not modify!\n\"\"\"\n\n'
          '_globals = {')

# Loop over all header files in "path" and in its subfolders
counter = 0
for path in PATHS:

    print(f'Parsing C-header files in {PATHS}...')

    for root, dirs, files in os.walk(path):

        for f in files:

            print(f'    Found \"{f}\", ', end='')

            # Skip all files that are not a header file
            if os.path.splitext(f)[1] != '.h':
                print('skipping...')
                continue
            else:
                print('parsing...')

            # Parse the header file "f" in "path"
            parser = pyclibrary.CParser(os.path.join(path, f),
                                        replace=REPLACE_DICT)

            # List of symbolic constants found
            symbols = list(parser.defs['values'].keys())

            if len(symbols) < 0:
                continue

            for symbol in symbols:

                # Skip any constants that do not start with any of the keywords
                # specified in "VALID_SUFFICES"
                stop = [''] * len(VALID_SUFFICES)
                for i, sufix in enumerate(VALID_SUFFICES):
                    if not symbol.startswith(sufix):
                        stop[i] = True
                    else:
                        stop[i] = False

                # Continue with the next "symbol" if the suffix of "symbol"
                # does not match any of "VALID_SUFFICES"
                if sum(stop) == len(VALID_SUFFICES):
                    continue

                # Value of the symbolic constant
                val = parser.defs['values'][symbol]

                out = f'\n    \'{symbol}\': '
                if isinstance(val, str):
                    out += f'\'{val}\''
                else:
                    out += f'{val}'
                out += ','
                fid.write(out)

            counter += 1

fid.write('\n'
          '}\n\n')
fid.close()

print(f'Number of C-header files parsed: {counter}.')
print('Done.')

