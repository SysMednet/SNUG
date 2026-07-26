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
        "-o",
        "--output_file",
        required=True,
        help="Output file (.txt)"
    )

    return parser.parse_args()

def timeuse(func) :
    def wrapper(*args,**kwargs) :
        import datetime
        start = datetime.datetime.now()
        res = func(*args,**kwargs)
        finish = datetime.datetime.now()
        print("%s\tfinish time:%s"%(func.__name__+'\t'+'\t'.join(str(s) for s in args),str(finish - start)))
        return res
    return wrapper

@timeuse
def ssn_define_mean(file_gem,file_w,amp,save):

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
    
    pat_d = {}
    with open(file_w,mode='r') as t_line :
        tem = t_line.readline()
        for n_line in t_line :
            tem = n_line.strip('\n').split('\t')
            pat_d[tem[0]] = float(tem[1])
    
    agg = np.corrcoef(value)
    # tri = np.triu_indices( gene_l , 1 )
    tri = np.full(agg.shape,False)
    tri[np.triu_indices(gene_l,1)] = True
    agg = agg[tri]
    agg[np.isnan(agg)] = 0
    pair_l = len(agg)
    
    start_p = 0
    stop_p = pat_l

    nor_sum , pair_num = 0 , 0
    for l,p in zip(range(start_p,stop_p),pat[start_p:stop_p]) : 
        value_s = np.corrcoef(np.c_[value,value[:,l]])
        value_s = value_s[tri]
        value_s[np.isnan(value_s)] = 0
        value_s = pat_d[p] * amp * pat_l * (value_s - agg) + agg
        
        #nor_mean += np.sum(value_s)
        nor_sum += np.sum(value_s)
        pair_num += pair_l

    nor_mean = nor_sum / pair_num
        
    with open(save, "w") as w_line:
        w_line.write("nor_mean\n")
        w_line.write(f"{nor_mean}\n")
        
if __name__ == '__main__' :
    args = parse_args()
    file_gem = args.gem
    file_w = args.weight
    amp = args.k
    save = args.output_file
    ssn_define_mean(file_gem,file_w,amp,save)