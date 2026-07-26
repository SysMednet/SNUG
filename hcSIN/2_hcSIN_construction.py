import networkx as nx
import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-n",
        "--network",
        required=True,
        help="Network file"
    )

    parser.add_argument(
        "-w",
        "--weight",
        required=True,
        help="Sample weight file (from step 1-1)"
    )

    parser.add_argument(
        "-s",
        "--SIN_folder",
        required=True,
        help="Path of SIN folder"
    )

    parser.add_argument(
        "-hc",
        "--hcSIN_folder",
        required=True,
        help="Path of hcSIN folder"
    )

    return parser.parse_args()

args = parse_args()
network = args.network
file_w = args.weight
SIN_folder = args.SIN_folder
hcSIN_folder = args.hcSIN_folder

with open(network,"r") as network_file :
    network_file = network_file.read().splitlines()
    string_net = set()
    for line in network_file[1:]:
        line = line.split("\t")
        if int(line[0]) > int(line[1]):
            string_net.add(line[1]+"\t"+line[0])
        else:
            string_net.add(line[0]+"\t"+line[1])
            
file_disease = open(file_w,"r")
disease_content = file_disease.read()
disease_content = disease_content.splitlines()
del disease_content[0]
Disease_list = []
for i in disease_content:
  i = i.split("\t")
  Disease_list.append(i[0])

for cell_line in Disease_list:
    with open(SIN_folder+"/0.01_"+str(cell_line)+"_network.txt","r") as single_net:
        new = open(hcSIN_folder+"/hcSIN_"+str(cell_line)+"_network.txt","w")
        single_net = single_net.read().splitlines()
        single_edge = []
        map_string = []
        single_node = set()
        map_node = set()
        for line in single_net:
            line = line.split("\t")
            single_node.add(line[0])
            single_node.add(line[1])
            if int(line[0]) > int(line[1]):
                single_edge.append(line[1]+"\t"+line[0])
                if line[1]+"\t"+line[0] in string_net:
                    map_string.append(line[1]+"\t"+line[0])
                    map_node.add(line[0])
                    map_node.add(line[1])
            else:
                single_edge.append(line[0]+"\t"+line[1])
                if line[0]+"\t"+line[1] in string_net:
                    map_string.append(line[0]+"\t"+line[1])
                    map_node.add(line[0])
                    map_node.add(line[1])

        for edge in map_string:
            new.write(edge+"\n")
            edge = edge.split("\t")

        new.close()
