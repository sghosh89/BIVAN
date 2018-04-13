import csv
import sys
from dendropy.dataio.nexusprocessing import escape_nexus_token


def read_csv(fn):
    rows = []
    with open(fn, 'rU') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            rows.append(row)
    header = rows.pop(0)
    return header, rows

out = sys.stdout
fn = sys.argv[1]
header, rows = read_csv(fn)  
nc = len(header) - 1 # first slot is species name, not a character
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
for label in header [1:]:
    out.write('\n    {}'.format(escape_nexus_token(label)))
    sys.stderr.write(',{}'.format(label))
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
    for cell in row[1:]:
        out.write('{}\t'.format(escape_nexus_token(cell)))
    out.write('\n')
out.write(' ;')

# terminate file
out.write('\nEND;\n')
