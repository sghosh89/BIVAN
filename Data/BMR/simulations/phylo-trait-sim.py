#!/usr/bin/env python
from __future__ import print_function
import random
import dendropy
import csv
import sys
import os

SCRIPT_NAME = os.path.split(sys.argv[0])[-1]
RNG = random.Random()

def warn(msg):
    sys.stderr.write('{} WARNING: {}\n'.format(SCRIPT_NAME, msg))

def debug(msg):
    sys.stderr.write('{} DEBUG: {}\n'.format(SCRIPT_NAME, msg))


def get_first_tree(tree_filepath):
    tree_list = dendropy.TreeList.get(path=tree_filepath,
                                      schema="nexus")
    if len(tree_list) != 1:
        m = 'expecting 1 tree in "{}" found {}'.format(tree_filepath, len(tree_list))
        if len(tree_list) > 1:
            warn(m)
        else:
            raise RuntimeError(m) 
    return tree_list[0]

def read_trait_offsets(csv_filepath):
    trait_pairs = []
    with open(csv_filepath, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='\"')
        ri = iter(reader)
        header = next(ri)
        for n, row in enumerate(ri):
            p = tuple([float(row[-2]), float(row[-1])])
            trait_pairs.append(p)
    return tuple(trait_pairs)

def do_sim(tree, trait_offsets):
    ntv = len(trait_offsets)
    assert ntv > 0
    nt = len(trait_offsets[0])
    out = sys.stdout
    tl = ["T{}".format(1 + i) for i in range(nt)]
    header_values = [''] + tl
    out.write('"{}"\n'.format('","'.join(header_values)))
    leaf_num = 1
    for nd in tree.preorder_node_iter():
        offset_idx = RNG.randrange(ntv)
        offset_pair = trait_offsets[offset_idx]
        par = nd.parent_node
        if par is None:
            nd.trait_values = offset_pair
        else:
            x = [par.trait_values[i] + offset_pair[i] for i in range(nt)]
            nd.trait_values = tuple(x)
        if nd.is_leaf():
            row = ['"{}"'.format(leaf_num)] + [str(v) for v in nd.trait_values]
            out.write('{}\n'.format(','.join(row)))
            leaf_num += 1

            


def main(tree_filepath, bivariate_traits_filepath):
    tree = get_first_tree(tree_filepath)
    trait_offsets = read_trait_offsets(bivariate_traits_filepath)
    do_sim(tree, trait_offsets)

if __name__ == '__main__':
    try:
        assert os.path.isfile(sys.argv[1])
        assert os.path.isfile(sys.argv[2])
        seed = int(sys.argv[3]) if len(sys.argv) > 3 else None
        if seed is not None:
            assert seed > 1
    except Exception as x:
        sys.exit('Expecting 2 mandatory arguments: filepaths to NEXUS tree and 3-column csv files with trait change values in the last 2 columns. A random number seed (integer > 1) can be sent in as an optional 3rd argument')
    try:
        if seed is None:
            seed = 1
            while seed < 2:
                seed = abs(hash(os.urandom(4)))
        debug('setting RNG seed={}'.format(seed))
        RNG.seed(seed)
        main(sys.argv[1], sys.argv[2])
    except Exception as x:
        sys.exit('Exception: {}\n{}\n'.format(type(x).__name__, str(x)))