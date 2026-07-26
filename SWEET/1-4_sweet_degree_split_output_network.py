# -*- coding: utf-8 -*-
import os
import numpy as np
import math
from scipy import stats
import datetime
import argparse

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-g",
        "--gem",
        required=True,
        help="Gene expression matrix (GEM)"
    )

    parser.add_argument(
        "--k",
        type=float,
        default=0.1
    )

    parser.add_argument(
        "-w",
        "--weight",
        required=True,
        help="Sample weight file (from step 1-1)"
    )

    parser.add_argument(
        "-m",
        "--mean",
        required=True,
        help="Sample mean file (from step 1-2)"
    )

    parser.add_argument(
        "-s",
        "--std",
        required=True,
        help="Sample std file (from step 1-3)"
    )

    parser.add_argument(
        "-o",
        "--output_folder",
        required=True,
        help="Path of output folder"
    )

    return parser.parse_args()
'''
def timeuse(func) :
    def wrapper(*args,**kwargs) :
        import datetime
        start = datetime.datetime.now()
        res = func(*args,**kwargs)
        finish = datetime.datetime.now()
        print("finish time: ",str(finish - start))  
        return res
    return wrapper
'''
#@timeuse
def ssn_define_degree_matrix(file_gem,file_w,amp,file_m,file_s,save) :
    os.makedirs(save, exist_ok=True)

    pvalue = np.array([math.pow(10,-2),])
    zscore = stats.norm.isf(pvalue/2)
    pvaluelen = len(pvalue)
    
    gene , value = [] , []
    with open(file_gem,mode='r') as r_line :
        pat = r_line.readline().strip('\n').split('\t')[1:]
        n_line = r_line.readline().strip('\n')
        while n_line :
            g , *val = n_line.split('\t')
            v_t = sum(np.array(val,dtype=float))
            if v_t != 0 :
                gene.append(g)
                value += val
            n_line = r_line.readline().strip('\n')
    gene_l , pat_l = len(gene) , len(pat)
    value = np.array(value,dtype=float).reshape(gene_l,pat_l)

    #eat weight file
    pat_d = {}
    with open(file_w,mode='r') as t_line :
        tem = t_line.readline()
        for n_line in t_line :
            tem = n_line.strip('\n').split('\t')
            pat_d[tem[0]] = float(tem[1])
            
    #eat mean file
    with open(file_m,mode='r') as t_line :
        tem = t_line.readline().strip('\n').split('\t')
        tem = t_line.readline().strip('\n').split('\t')
        nor_m = float(tem[0])
     
    #eat std file
    with open(file_s,mode='r') as t_line :
        tem = t_line.readline().strip('\n').split('\t')
        tem = t_line.readline().strip('\n').split('\t')
        nor_s = float(tem[0])

    start_p = 0
    stop_p = pat_l
    cut_p = pat_l

    agg = np.round(np.corrcoef(value),6)
    tri = np.full(agg.shape,False)
    tri[np.triu_indices(gene_l,1)] = True
    agg = agg[tri]
    agg[np.isnan(agg)] = 0

    p1 = np.full((gene_l,gene_l),False)
    for l,pa,c in zip(range(start_p,stop_p),pat[start_p:stop_p],range(0,cut_p)) :
        value_s = np.c_[value,value[:,l]]
        value_s = np.round(np.corrcoef(value_s)[tri],6)
        value_s[np.isnan(value_s)] = 0
        value_s = pat_d[pa] * amp * pat_l * (value_s - agg) + agg
        value_s = np.abs((value_s-nor_m)/nor_s)
        for p , z , r in zip(pvalue,zscore,range(pvaluelen)) :
            temp = np.full((gene_l, gene_l), False, dtype=bool)
            temp[tri] = value_s >= z
            p1 = temp | temp.T
            edge_set = set()
            n=0
            new_file = open(save+"/"+str(p)+"_"+str(pa)+"_network.txt","w")
            for number in range(np.shape(p1)[1]):
                for loc in range(len(p1[number])):
                    if p1[number][loc] == True :
                        n+=1
                        if int(gene[number]) > int(gene[loc]):
                            edge_set.add(gene[loc]+"\t"+gene[number]+"\n")
                        else:
                            edge_set.add(gene[number]+"\t"+gene[loc]+"\n")
            print(len(edge_set),n/2)
            for pair in edge_set:
                new_file.write(pair)
            new_file.close()
  
if __name__ == '__main__' :
    args = parse_args()
    file_gem = args.gem
    file_w = args.weight
    amp = args.k
    file_m = args.mean
    file_s = args.std
    save = args.output_folder
    ssn_define_degree_matrix(file_gem,file_w,amp,file_m,file_s,save)