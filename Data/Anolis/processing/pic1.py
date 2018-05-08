import dendropy
from dendropy.model import continuous
import sys

# Taken from file of same name in DendroPy source/examples and slightly

char_fn, tree_fn = sys.argv[1:]

taxa = dendropy.TaxonNamespace()
tree = dendropy.Tree.get(
        path=tree_fn,
        schema="nexus",
        taxon_namespace=taxa)


chars = dendropy.ContinuousCharacterMatrix.get_from_path(
        char_fn,
        "nexus",
        taxon_namespace=taxa)

crows = []
for cind, nd in enumerate(tree.postorder_internal_node_iter()):
    if len(crows) <= cind:
        crows.append([])
    crow = crows[cind]
    left_child, right_child = nd.child_nodes()
    crow.append(left_child.edge.length + right_child.edge.length)

pic = continuous.PhylogeneticIndependentConstrasts(
        tree=tree,
        char_matrix=chars)
        

for cidx in range(chars.vector_size):
    ctree1 = pic.contrasts_tree(character_index=cidx,
            annotate_pic_statistics=True,
            state_values_as_node_labels=True,
            corrected_edge_lengths=False)
    for cind, nd in enumerate(ctree1.postorder_internal_node_iter()):
        crow = crows[cind]
        crow.append(nd.pic_contrast_standardized)

out = sys.stdout
for cind, crow in enumerate(crows):
    out.write('c{},'.format(cind))
    out.write(','.join(['{:8.6f}'.format(i) for i in crow]))
    out.write('\n')
