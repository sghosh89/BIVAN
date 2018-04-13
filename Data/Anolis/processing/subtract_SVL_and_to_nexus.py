#!/usr/bin/env python
from __future__ import print_function
import csv
import sys
from dendropy.dataio.nexusprocessing import escape_nexus_token


def read_csv_and_convert_non_first_to_floats(fn):
    rows = []
    with open(fn, 'rU') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            rows.append(row)
    header = rows.pop(0)
    for i in range(len(rows)):
        row = rows[i]
        for j in range(1, len(row)):
            row[j] = float(row[j])
    return header, rows

out = sys.stdout
fn = sys.argv[1]
header, rows = read_csv_and_convert_non_first_to_floats(fn)
svl_inds = [i for i, h in enumerate(header) if 'SVL' in h]
first_svl = svl_inds[0]
avoid = set(svl_inds)
nc = len(header) - len(avoid) - 1 # first slot is species name, not a character
nt = len(rows)

# Write preamble...
out.write('''#NEXUS
BEGIN DATA;
    DIMENSIONS  NTAX={} NCHAR={};
    FORMAT DATATYPE = CONTINUOUS GAP = - MISSING = ?;
'''.format(nt, nc))

# Write CharLabels...
out.write('CharLabels')
sys.stderr.write('contrast,path length')
for col_ind, label in enumerate(header):
    if col_ind == 0 or col_ind in avoid:
        continue
    corrected_lab = label + ' minus SVL'
    out.write('\n    {}'.format(escape_nexus_token(corrected_lab)))
    sys.stderr.write(',{}'.format(corrected_lab))
sys.stderr.write('\n')
# Write Matrix...
out.write(''';
MATRIX
''')
## compose a format string to right-pad with spaces all shorter OTU names
max_tax_label = max([len(row[0]) for row in rows])
nfs = '{:<' + str(max_tax_label + 1) + '} '
for row in rows:
    out.write(nfs.format(row[0]))
    for col_ind, cell in enumerate(row):
        if col_ind == 0 or col_ind in avoid:
            continue
        d = cell - row[first_svl]
        esd = str(d)
        out.write('{}\t'.format(esd))
    out.write('\n')
out.write(' ;')

# terminate file
out.write('\nEND;\n')
