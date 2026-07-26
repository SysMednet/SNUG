# -*- coding: utf-8 -*-
import numpy as np
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
        "-c",
        "--corr",
        required=True,
        help="Sample to sample correlation matrix (.txt)"
    )

    parser.add_argument(
        "-w",
        "--weight",
        required=True,
        help="Sample weight file (.txt)"
    )

    return parser.parse_args()

def sample_corr(file1,save) :
    
    value = []
    gene = []
    
    with open(file1,mode='r') as t_line :
        pat = t_line.readline().strip('\n').split('\t')[1:]
        pat_l = len(pat)
        for n_line in t_line :
            tem = n_line.strip('\n').split('\t')
            value += tem[1:]
            gene.append(tem[0])
    gene_l = len(gene)
    value = np.array(value,dtype=float).reshape(gene_l,pat_l).T
    
    pat_c = np.corrcoef(value)
    
    with open(save,mode='w') as w_line :
        w_line.write('corr\t'+'\t'.join( n for n in pat)+'\n')
        for n1,p in zip(pat_c,pat) :
            tem = p + '\t' + '\t'.join( str(s) for s in n1 ) + '\n'
            w_line.write(tem)
            
'''
(x-min+0.01)/(max-min+0.01)
'''
def corr_maxmin(file1,save) :
    
    value = []
    
    with open(file1,mode='r') as t_line :
        pat = t_line.readline().strip('\n').split('\t')[1:]
        pat_l = len(pat)
        for n_line in t_line :
            tem = n_line.strip('\n').split('\t')
            value += tem[1:]

    value = np.array(value,dtype=float).reshape(pat_l,pat_l)
    value = (np.sum(value,axis=1)-1)/(pat_l-1)
    rmax , rmin = np.max(value) , np.min(value)
    dif = rmax - rmin + 0.01
    value = (value - rmin + 0.01)/dif
    
    with open(save,mode='w') as w_line :
        w_line.write('Sample\tcorr\n')
        for p,v in zip(pat,value) :
            tem = p + '\t' + str(v) + '\n'
            w_line.write(tem)

if __name__ == '__main__' :
    args = parse_args()
    file_gem = args.gem
    file_corr = args.corr
    file_w = args.weight
    sample_corr(file_gem,file_corr)
    corr_maxmin(file_corr,file_w)