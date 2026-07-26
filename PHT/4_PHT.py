#python3 4_PHT.py -o SNUG.txt -hc ./path/hcSIN -i individual_signature.txt -d drug_target.txt
import networkx as nx
from scipy.stats import hypergeom
import argparse

def parse_args():
    parser = argparse.ArgumentParser()

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
        help="Path of hcSIN folder"
    )

    parser.add_argument(
        "-i",
        "--IS_file",
        required=True,
        help="Individual signature file (from step 3)"
    )

    parser.add_argument(
        "-d",
        "--drug_target",
        required=True,
        help="Drug target file (.txt)"
    )

    return parser.parse_args()

args = parse_args()
file_name = args.output_file
hcSIN_folder = args.hcSIN_folder
individual_signature_file = args.IS_file
drug_target_file = args.drug_target

#print(time_start)
file_disease = open(individual_signature_file,"r")
disease_content = file_disease.read()
disease_content = disease_content.splitlines()
del disease_content[0]
DIC_disease = {}
for i in disease_content:
  i = i.split("\t")
  lis = []
  DIC_disease[i[0]] = i[-1].split(";")

'''
record_pair = open("../../IC50/GDSC_combine.txt","r")
file = record_pair.read()
file = file.splitlines()
disease_set = set()
del file[0]
all_record_pair = set()
for line in file:
    line = line.split("\t")
    all_record_pair.add(line[0]+"\t"+line[1])
    disease_set.add(line[0])
'''

file_drug = open(drug_target_file,"r")
drug_content = file_drug.read()
drug_content = drug_content.splitlines()
DIC_Drug = {}
Drug_list = []
for i in drug_content:
  i = i.split("\t")
  Drug_list.append(i[0])
  DIC_Drug[i[0]] = i[-1].split(";")

file = open(file_name,"w")
file.write("Drug"+"\t"+"Sample"+"\t"+"m"+"\t"+"n"+"\t"+"bigM"+"\t"+"bigN"+"\t"+"p_value"+"\n")
for cell_line in list(DIC_disease.keys()):
    #if cell_line in disease_set:
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
        max_dis = int(nx.diameter(G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])))
        for a in range(0,len(sorted(nx.connected_components(G), key=len, reverse=True))):
            G_sub_big = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[a])
            for node,degree in nx.degree(G_sub_big):
                gene_component[node] = a
                if a == 0 :
                    list_node.append(node)
        distance_dic = {}
        all_dis = []
        for gene1 in all_node:
            for gene2 in all_node:
                lis = [int(gene1),int(gene2)]
                lis.sort()
                if str(lis[0])+"\t"+str(lis[1]) in distance_dic:
                    all_dis.append(int(distance_dic[str(lis[0])+"\t"+str(lis[1])]))
                    continue
                elif gene_component[gene1] == gene_component[gene2]:
                    distance_dic[str(lis[0])+"\t"+str(lis[1])] = int(nx.shortest_path_length(G, gene1,gene2))
                    all_dis.append(distance_dic[str(lis[0])+"\t"+str(lis[1])])
                elif  gene_component[gene1] != gene_component[gene2] :
                    distance_dic[str(lis[0])+"\t"+str(lis[1])] = max_dis 
                    all_dis.append(max_dis )
        bigN = int(sum(all_dis))
        gene_to_network = []
        for target in DIC_disease[cell_line]:
            for gene in all_node:
                lis = [int(target),int(gene)]
                lis.sort()
                gene_to_network.append(int(distance_dic[str(lis[0])+"\t"+str(lis[1])]))
        n = int(sum(gene_to_network))
        for drug in Drug_list:
            #if cell_line+"\t"+drug in all_record_pair: 
            drug_list = []
            for gene in DIC_Drug[drug]:
                if gene in all_node:
                    drug_list.append(gene)
                else:
                    continue
            if len(drug_list) > 0 :
                target_to_network = []
                target_to_disease =[]
                for target in drug_list:
                    for gene in all_node:
                        lis = [int(target),int(gene)]
                        lis.sort()
                        target_to_network.append(int(distance_dic[str(lis[0])+"\t"+str(lis[1])]))
                for target in drug_list:
                    for gene in DIC_disease[cell_line]:
                        lis = [int(target),int(gene)]
                        lis.sort()
                        target_to_disease.append(int(distance_dic[str(lis[0])+"\t"+str(lis[1])]))
                bigM = int(sum(target_to_network))
                m = int(sum(target_to_disease))
                p_value = float(hypergeom.cdf(m,bigN,bigM,n))
                file.write(str(drug)+"\t"+str(cell_line)+"\t"+str(m)+"\t"+str(n)+"\t"+str(bigM)+"\t"+str(bigN)+"\t"+str(p_value)+"\n")
file.close()