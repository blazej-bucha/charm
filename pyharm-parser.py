# Extracts symbolic constants from the CHarm's header files that are inside
# "PATH" and saves the constants to a Python "FILE" to be used with PyHarm.

import os
import pyclibrary

# Path to read the header files from
PATH = './charm/'

# Name of the output Python file to save the symbolic constants to
FILE = './wrap/pyharm/_constants.py'

# Create the output file to write the symbolic constants to
fid = open(FILE, 'w')

fid.write('\"\"\"\nGenerated automatically.  Do not modify!\n\"\"\"\n\n'
          '_globals = {')

print(f'Parsing C-header files in {PATH}...')

# Loop over all header files in "path" and in its subfolders
counter = 0
for root, dirs, files in os.walk(PATH):

    for f in files:

        print(f'    Found \"{f}\", ', end='')

        # Skip all files that are not a header file
        if os.path.splitext(f)[1] != '.h' or (not f.startswith('charm_')  \
            and not f == 'charm.h'):
            print('skipping...')
            continue
        else:
            print('parsing...')

        # Parse the header file "f" in "PATH"
        parser = pyclibrary.CParser(os.path.join(PATH, f))

        # List of symbolic constants found
        symbols = list(parser.defs['values'].keys())

        if len(symbols) < 0:
            continue

        for symbol in symbols:

            # Skip any constants that do not start with the "CHARM_" prefix
            if not symbol.startswith('CHARM_'):
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

