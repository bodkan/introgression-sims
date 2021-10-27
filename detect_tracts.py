import argparse
import pathlib
import sys
import os

import tskit
import pyslim
import pandas
import numpy
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--trees", required=True, metavar="FILE", help="Tree sequence file")
parser.add_argument("--model", required=True, metavar="DIRECTORY", help="slendr model directory")
parser.add_argument("--source-label", required=True, type=str,
                    help="Introgressing population label")
parser.add_argument("--output", required=True, metavar="FILE",
                    help="Output file with coordinates of introgressed tracts")
args = parser.parse_args("--trees results/output_ts.trees --model model/ --source-label NEA --output out.txt".split())

model_dir = os.path.expanduser(args.model)

if not os.path.exists(model_dir):
    sys.exit(f"Model directory {model_dir} does not exist")

populations_path = pathlib.Path(model_dir, "populations.tsv")
source_id = pandas.read_table(populations_path).query(f"pop == '{args.source_label}'").pop_id[0]

ts = pyslim.load("test.trees")

# Load the .trees file and assess true local ancestry
breaks = numpy.zeros(ts.num_trees + 1)
ancestry = numpy.zeros(ts.num_trees + 1)

tree = ts.first()
# iterate over trees and ...
for tree in ts.trees():
    subpop_sum, subpop_weights = 0, 0
    # ... for each tree, check in which population does its root originate
    for root in tree.roots:
        leaves_count = tree.num_samples(root) - 1  # subtract one for the root, which is a sample
        subpop_sum += tree.population(root) * leaves_count
        subpop_weights += leaves_count
    breaks[tree.index] = tree.interval[0]
    ancestry[tree.index] = subpop_sum / subpop_weights
breaks[-1] = ts.sequence_length
ancestry[-1] = ancestry[-2]

# Make a simple plot
plt.plot(breaks, ancestry)
plt.show()