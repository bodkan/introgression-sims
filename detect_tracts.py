import argparse
from pathlib import Path

import random
from itertools import chain
from collections import defaultdict

import tskit
import pyslim
import pandas
import numpy
import tqdm

def get_introgressed_leaves(tree, neand_id):
    """Get SLiM IDs of nodes coalescing with the Neanderthal
    population root node in a given tree.
    """
    start, end = tree.interval

    admixed_list = []
    for root in tree.roots:
        # is the current root rooted in a Neanderthal population?
        if tree.population(root) == neand_id:
            # if it is, collect all non-Neanderthal nodes under this root
            # (these are Eurasians carrying Neanderthal introgression
            # at this locus)
            admixed = [node for node in tree.leaves(root)
                       if tree.population(node) != neand_id]
            admixed_list.extend(admixed)

    return admixed_list

def detect_introgression(ts, neand_id):
    """Iterate over all trees in the tree sequence and look for
    Eurasian node tracts which have been introgressed from
    Neanderthals.
    """
    result = []

    for tree in tqdm.tqdm(ts.trees()):
        leaves = get_introgressed_leaves(tree, neand_id)
        start, end = tree.interval
        result.append((int(start), int(end), leaves))

    tracts = pandas.DataFrame(result, columns=["start", "end", "node"]) \
        .explode("node") \
        .dropna() \
        .sort_values(by=["node", "start"]) \
        .reset_index(drop=True) \
    
    return tracts

def merge_adjacent(tracts):
    """Merge all detected introgressed tracts for each node
    (i.e. chromosome).
    """
    all_merged = []

    for node in tqdm.tqdm(tracts.node.unique()):
        node_tracts = []

        for i, seg in enumerate(tracts.query(f"node == {node}").itertuples()):
            if i == 0: # initialize detection of continuous tracts with the first tract
                tract_start, prev_end = seg.start, seg.end
            elif seg.start == prev_end: # the current tract starts where
                prev_end = seg.end      # the previous ended -- continue
            else: # end the continuous tract, start a new one
                node_tracts.append((tract_start, prev_end, node))
                tract_start, prev_end = seg.start, seg.end

        # close the last remaining tract
        node_tracts.append((tract_start, seg.end, node))

        all_merged.append(node_tracts)

    all_merged = list(chain.from_iterable(all_merged))

    tracts = pandas.DataFrame(all_merged, columns=["start", "end", "node"]) \
        .sort_values(by=["node", "start"]) \
        .reset_index(drop=True)

    return tracts

#
# parse command line arguments
#
parser = argparse.ArgumentParser(
    "Script for extracting Neanderthal tracts from simulated "
    "tree sequence data"
)

parser.add_argument("--slendr", required=True, metavar="DIRECTORY",
                    help="Location of the slendr model directory")
parser.add_argument("--trees", required=True, metavar="FILE",
                    help="Location of the .trees tree sequence file")
parser.add_argument("--output", required=True, metavar="FILE",
                    help="Where to save .tsv.gz table with coordinates of introgressed tracts")
args = parser.parse_args()

model_dir = Path(args.slendr)

#
# load and process all required data
#

# read information about populations, their names, and SLiM IDs
populations = pandas.read_table(model_dir / "populations.tsv")

print("Loading tree sequence input data...")

# load the complete tree sequence output
ts = tskit.load(args.trees)

# get nodes (chromosomes) of remembered individuals and simplify
# the tree sequence only to those nodes
remembered_nodes = list(chain.from_iterable(
    ind.nodes for ind in ts.individuals()
    if ind.flags & pyslim.INDIVIDUAL_REMEMBERED
))
ts = ts.simplify(remembered_nodes)

#
# detect introgressed tracts
#
neand_id = populations.query("pop == 'NEA'").pop_id.values[0]

print("Detecting introgressed tracts...")
tracts = detect_introgression(ts, neand_id=neand_id)

print("Merging adjacent tracts in each node...")
merged_tracts = merge_adjacent(tracts)

#
# save the results
#

# before we can save the table of Neanderthal tracts, we need to
# assign to each node/chromosome its unique numeric SLiM id to be
# able to trace those tracts back to the easy-to-read individual
# names we have in slendr
merged_tracts["slim_id"] = pandas.Series(
    ts.node(i).metadata["slim_id"]
    for i in merged_tracts.node
)

print(f"Saving the output to {args.output}")
merged_tracts[["start", "end", "slim_id"]] \
    .to_csv(args.output, sep="\t", index=False)
