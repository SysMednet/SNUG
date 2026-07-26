# -*- coding: utf-8 -*-
#python3 3_individual_signature_selection.py -w weight.txt -o individual_signature.txt -hc ./path/hcSIN
import sys
import networkx as nx
import numpy
import numpy as np
import random
import math
from scipy.stats import hypergeom
from kneed import DataGenerator, KneeLocator
import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-w",
        "--weight",
        required=True,
        help="Sample weight file (from step 1-1)"
    )

    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output file (.txt)"
    )

    parser.add_argument(
        "-hc",
        "--hcSIN_folder",
        required=True,
        help="Path of hcSIN_folder"
    )

    return parser.parse_args()

args = parse_args()
file_w = args.weight
file_name = args.output_file
hcSIN_folder = args.hcSIN_folder

file_disease = open(file_w,"r")
disease_content = file_disease.read()
disease_content = disease_content.splitlines()
del disease_content[0]
DIC_Disease = {}
Disease_list = set()
DIC_Disease_EZ = {}
for i in disease_content:
  i = i.split("\t")
  lis = []
  Disease_list.add(i[0])

for p in ["0.01"]:
    new = open(file_name,"w")
    new.write("Sample\tElbow point of degree\tNumber of genes\tGene list\n")
    for cell_line in Disease_list:
        with open(hcSIN_folder+"/hcSIN_"+str(cell_line)+"_network.txt","r") as single_net:
            single_net = single_net.read().splitlines()
            edge_list = []
            all_node = set()
            for i in single_net:
                i = i.split("\t")
                edge_list.append(i[0]+"\t"+i[1])
                all_node.add(i[0])
                all_node.add(i[1])
            G = nx.parse_edgelist(edge_list,delimiter = '\t',nodetype = str)
            gene_component = {}
            list_node = []
            Node_biggest_sub = {}
            for node,degree in nx.degree(G):
                Node_biggest_sub[str(node)] = degree
        score_node_sort = sorted(list(Node_biggest_sub.values()),reverse =False)
        node_sort= sorted(Node_biggest_sub.items(), key = lambda x:x[1] ,reverse = True)
        x = list(n for n in range(1, len(score_node_sort)+1))
        y = score_node_sort
        kneedle = KneeLocator(x, y, S=1.0, curve="convex", direction="increasing")
        cutoff = int(kneedle.knee_y)
        lis = []
        for node in node_sort:
            if int(node[1]) >= cutoff:
                lis.append(node[0])
            else:
                break
        new.write(str(cell_line)+"\t"+str(cutoff)+"\t"+str(len(lis))+"\t"+";".join(lis)+"\n")
    new.close()