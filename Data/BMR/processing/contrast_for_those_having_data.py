#!/usr/bin/env python
from __future__ import print_function
import dendropy
from dendropy.model import continuous
import sys, csv, math
from dendropy.dataio.nexusprocessing import escape_nexus_token


# Taken from file of same name in DendroPy source/examples and slightly

def print_to_root(nd, s=0, i=0):
    blen = nd.edge.length
    s += blen
    msg = '{}{} (sum={}, i={})\n'
    msg = msg.format('  '*i, blen, s, i)
    sys.stderr.write(msg)
    if nd.parent_node is not None:
        print_to_root(nd.parent_node, s, 1 + i)

def main(char_csv, tree_fn):
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(
            path=tree_fn,
            schema="nexus",
            taxon_namespace=taxa)
    headers, rows = read_csv(char_csv)
    assert(headers[1].lower() == 'mass')
    assert(headers[2].lower() == 'bmr')
   
    # tnd = tree.find_node_with_taxon_label('Neomys fodiens')
    # print_to_root(tnd, 0.0)
    
    data_dict = {}

    for row in rows:
        label = row[0]
        mass = float(row[1])
        bmr = float(row[2])

        if label in data_dict:
            sys.exit('label "{}" is repeated\n'.format(label))
        data_dict[label] = (mass, bmr, math.log(mass))
    to_prune = [i.taxon for i in tree.leaf_node_iter() if i.taxon.label not in data_dict]
    # print('to_prune = {}'.format(len(to_prune)))
    # print('to_keep = {}'.format(len(taxa) - len(to_prune)))
    # print('data_dict = {}'.format(len(data_dict)))
    tree.prune_taxa(to_prune)
    for t in to_prune:
        taxa.remove_taxon(t)
    
    # tnd = tree.find_node_with_taxon_label('Neomys fodiens')
    # tree.print_plot()
    # print_to_root(tnd, 0.0)
    
    cm = dendropy.ContinuousCharacterMatrix(taxon_namespace=taxa)
    for taxon in taxa:
        cm.new_sequence(taxon, data_dict[taxon.label])

    crows = []
    for cind, nd in enumerate(tree.postorder_internal_node_iter()):
        if len(crows) <= cind:
            crows.append([])
        crow = crows[cind]
        left_child, right_child = nd.child_nodes()
        crow.append(left_child.edge.length + right_child.edge.length)
    pic = continuous.PhylogeneticIndependentConstrasts(
            tree=tree,
            char_matrix=cm)

    for cidx in range(cm.vector_size):
        ctree1 = pic.contrasts_tree(character_index=cidx,
                                    annotate_pic_statistics=True,
                                    state_values_as_node_labels=True,
                                    corrected_edge_lengths=False)
        for cind, nd in enumerate(ctree1.postorder_internal_node_iter()):
            crow = crows[cind]
            crow.append(nd.pic_contrast_standardized)

    out = sys.stdout
    out.write('contrast,path length,mass,BMR,ln mass\n')
    for cind, crow in enumerate(crows):
        out.write('c{},'.format(cind))
        out.write(','.join(['{:8.6f}'.format(i) for i in crow]))
        out.write('\n')

    return tree
    
def read_csv(fn):
    rows = []
    with open(fn, 'rU') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            rows.append(row)
    header = rows.pop(0)
    return header, rows

if __name__ == '__main__':
    assert len(sys.argv) == 3
    tree = main(sys.argv[1], sys.argv[2])

"""
chars = dendropy.ContinuousCharacterMatrix.get_from_path(
        char_fn,
        "nexus",
        taxon_namespace=taxa)




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
"""
