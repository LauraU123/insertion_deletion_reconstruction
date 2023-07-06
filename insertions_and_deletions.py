from Bio import Phylo, SeqIO
from itertools import groupby, count
from string import digits
import collections
import pandas as pd
import re, os, json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import argparse


def feature_dictionary(refseq):

    """returns a dictionary of CDS and their locations"""

    ref = SeqIO.read(refseq, 'genbank')
    dict_of_genes = dict()
    for feature in ref.features:
        if feature.type == "CDS": dict_of_genes[feature.qualifiers["gene"][0]] = feature.location

    return(dict_of_genes)


def deletion_recursive(node, list_=None, dictionary_=None):

    """ this function returns a dictionary with node name as key and a list of mutations along that branch as the info"""

    if list_ is None:
        list_ = []
    if dictionary_ is None:
        dictionary_ = dict()

    if 'mutations' in node['branch_attrs']:
        if 'nuc' in node['branch_attrs']['mutations']:
            list_=(node['branch_attrs']['mutations']['nuc'])
            
    if 'name' in node:
            dictionary_[node['name']] = list(list_)

    if 'children' in node:
        for child in node['children']:
           deletion_recursive(child, list(list_), dictionary_)
           
    return(dictionary_)


def reconstruct_insertions(insertionfile, tree):

    """ This function reconstructs insertions along a newick tree file, provided an insertion CSV file for terminal nodes"""

    tree_ = Phylo.read(tree, "newick")
    reconstruct_insertions = dict()
    csv_ = pd.read_csv(insertionfile)
    tree_.root_at_midpoint()
    tree_.find_clades()
    
    for branch in tree_.get_nonterminals(order='postorder'):
        shared_insertions = set()
        branch_ = []

        for b in branch:
            if b.is_terminal():
                insertion = csv_.loc[csv_['seqName'] == b.name]['insertions'].values.item()
                if pd.isna(insertion):  reconstruct_insertions[b.name] = set()
                else: reconstruct_insertions[b.name] = set(insertion.split(';'))
            branch_.append(reconstruct_insertions[b.name])
        shared_insertions = set.intersection(*branch_)

        for b in branch:  
            reconstruct_insertions[b.name] = reconstruct_insertions[b.name].difference(shared_insertions)
        reconstruct_insertions[branch.name] = shared_insertions
    return(reconstruct_insertions)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Reconstructing Insertions and Deletions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--reference", required=True, help="genbank reference")
    parser.add_argument("--json", required=True, help="json file")
    parser.add_argument("--output", help="tsv output file")
    parser.add_argument("--insertions", help="insertions_file")
    parser.add_argument("--tree", help="tree nwk file")
    args = parser.parse_args()

    # Reconstructing insertions, deletions and finding relevant indices for G
    with open (args.json) as file_:
        f = json.load(file_)  
        del_by_node = deletion_recursive(f['tree'])
    ins_by_node = reconstruct_insertions(args.insertions, args.tree)
    G_gene = feature_dictionary(args.reference)["G"]


    #Finding all branches with insertions 
    insertions_dict, deletions_dict = (defaultdict(list) for i in range(2))
    for branch, insertions in ins_by_node.items():
        if insertions != set():
            for i in insertions:
                location = int(i.split(":")[0])
                if location in G_gene: insertions_dict[branch].append(i)

    #Finding all branches with deletions
    for branch, deletions in del_by_node.items():
        if deletions != set():
            for i in deletions:
                location = int(i[1:-1])
                if location in G_gene: deletions_dict[branch].append(i)

    # Finding deletions and insertion locations and lengths
    ins_length, del_length = (dict() for i in range(2))
    ins_location, del_loc = (defaultdict(list) for i in range(2))
    for branch, insertion in insertions_dict.items():
        ins_length[branch] = len(insertion)
        for i in insertion: ins_location[branch].append(int(i.split(":")[0]))

    for branch, deletion in deletions_dict.items():
        del_length[branch] = len(deletion)
        for d in deletion: del_loc[branch].append(int(d[1:-1]))

    #Constructing TSV output file for insertions and deletions
    dictionaries = [del_length, ins_length, del_loc, ins_location]
    df = pd.DataFrame(dictionaries)
    df.index =['deletion length', 'insertion length', 'deletion locations', 'insertion locations']
    df = df.T 
    df =df.reset_index()
    df.to_csv(args.output, sep='\t')