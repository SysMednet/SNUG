# SNUG

SNUG is a single-sample network-based drug efficacy prediction framework that improves traditional network proximity analysis by integrating individual signatures (IS), high-confidence single-sample networks (hcSINs) and proximity-based hypergeometric test (PHT) scoring.  

The SNUG workflow consists of four main steps:
* Step 1: single-sample networks (SINs) construction
* Step 2: high-confidence single-sample networks (hcSINs) construction
* Step 3: individual signature (IS) selection
* Step 4: proximity-based hypergeometric test (PHT)

# Dependencies

The code was developed and tested with the following software:
- Python 3.14.8
- NumPy 2.4.3
- SciPy 1.17.1
- NetworkX 3.6.1
- kneed 0.8.5  

Install the required packages using:

```bash
pip install "package name"
```

# Input File Format

### Gene expression matrix (GEM) file

The input GEM file should be a tab-separated `.txt` file with the following structure:
| gene id | sample 1 | sample 2 | ... | sample n |
| ------- | --- | --- | --- | --- |
| gene 1 | ... | ... | ... | ... |
| gene 2 | ... | ... | ... | ... |
| ... | ... | ... | ... | ... |
| gene n | ... | ... | ... | ... |
| ... | ... | ... | ... | ... |

Note: Gene IDs should be converted to Entrez Gene IDs before running SNUG.

### Background network file

The input network file should be a tab-separated `.txt` file with the following columns:  

`Entrez1	Entrez2	protein1	protein2	combined_score`

### Drug target file

The input drug target file should be a `.txt` file in the following format (without a header). Each row should contain a drug name followed by its target genes.

```text
drug 1    target 1;target 2; ... ;target n; ...
drug 2    target 1;target 2; ... ;target n; ...
...
drug n    target 1;target 2; ... ;target n; ...
...
```

Note:
- Target genes should be converted to Entrez Gene IDs before running SNUG.
- The targets for each drug should be separated by `;`.

# Basic Usage

The demo dataset is a small example for testing the workflow and does not represent the complete dataset used in the manuscript.

### Step 1: single-sample networks (SINs) construction

**Step 1-1: calculate genome-wide sample weights:**

```bash
python3 ./SWEET/1-1_correlation_to_weight.py -g ./demo/demo_GEM.txt -c ./demo/demo_correlation.txt -w ./demo/demo_weight.txt
```
`-h`: Get help with the commands.  
`-g`: Gene expression matrix (GEM) file.  
`-c`: Output file containing sample to sample correlation matrix.  
`-w`: Output file containing weight of each sample.  

**Step 1-2: Calculate the mean value of all samples.**

```bash
python3 ./SWEET/1-2_sweet_mean.py -g ./demo/demo_GEM.txt -w ./demo/demo_weight.txt -o ./demo/demo_mean.txt
```
`-h`: Get help with the commands.  
`-g`: Gene expression matrix (GEM) file.  
`-w`: Sample weight file (i.e., the output file from step 1-1).  
`-o`: Output file containing mean value of all samples.  
`--k`: Balance parameter (default: 0.1).  

**Step 1-3: Calculate the standard deviation (std) of all samples.**

```bash
python3 ./SWEET/1-3_sweet_std.py -g ./demo/demo_GEM.txt -w ./demo/demo_weight.txt -m ./demo/demo_mean.txt -o ./demo/demo_std.txt
```
`-h`: Get help with the commands.  
`-g`: Gene expression matrix (GEM) file.  
`-w`: Sample weight file (i.e., the output file from step 1-1).  
`-m`: Sample mean value file (i.e., the output file from step 1-2).  
`-o`: Output file containing standard deviation of all samples.  
`--k`: Balance parameter (default: 0.1).  

**Step 1-4: Construct SINs of each sample.**

```bash
python3 ./SWEET/1-4_sweet_degree_split_output_network.py -g ./demo/demo_GEM.txt -w ./demo/demo_weight.txt -m ./demo/demo_mean.txt -s ./demo/demo_std.txt -o ./demo/demo_SIN
```
`-h`: Get help with the commands.  
`-g`: Gene expression matrix (GEM) file.  
`-w`: Sample weight file (i.e., the output file from step 1-1).  
`-m`: Sample mean value file (i.e., the output file from step 1-2).  
`-s`: Sample standard deviation file (i.e., the output file from step 1-3).  
`-o`: Path to the folder containing SIN files.  

### Step 2: high-confidence single-sample networks (hcSINs) construction

```bash
python3 ./hcSIN/2_hcSIN_construction.py -n ./demo/demo_network.txt -w ./demo/demo_weight.txt -s ./demo/demo_SIN -hc ./demo/demo_hcSIN
```
`-h`: Get help with the commands.  
`-n`: Background network file.  
`-w`: Sample weight file (i.e., the output file from step 1-1).  
`-s`: Path to the folder containing SIN files (i.e., the output file from step 1-4).  
`-hc`: Path to the folder containing hcSIN files.  

### Step 3: individual signature (IS) selection

```bash
python3 ./IS/3_individual_signature_selection.py -w ./demo/demo_weight.txt -hc ./demo/demo_hcSIN -o ./demo/demo_individual_signature.txt
```
`-h`: Get help with the commands.  
`-w`: Sample weight file (i.e., the output file from step 1-1).  
`-hc`: Path to the folder containing hcSIN files (i.e., the output file from step 2).  
`-o`: Output file containing individual signatures of each sample.  

### Step 4: proximity-based hypergeometric test (PHT)

```bash
python3 ./PHT/4_PHT.py -hc ./demo/demo_hcSIN -i ./demo/demo_individual_signature.txt -d ./demo/demo_drug_target.txt -o ./demo/demo_SNUG.txt
```
`-h`: Get help with the commands.  
`-hc`: Path to the folder containing hcSIN files (i.e., the output file from step 2).  
`-i`: Individual signature file (i.e., the output file from step 3).  
`-d`: Drug target file.  
`-o`: Output file containing the PHT results for each sample.  