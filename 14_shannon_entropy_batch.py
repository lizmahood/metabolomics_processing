#!/usr/bin/env python

from __future__ import division
import sys, math, glob, os, re


def get_groups(coln):
    '''
    coln is a list of strings
    This function is hard coded to look for groups at
    [32:-1] of this list
    RETURNS: dict of columns in each condition
    '''
    outd = {}
    condl = []
    newl = [i.split('_')[0] for i in coln[32:-1]]
    for i in newl:
        if i not in condl:
            condl.append(i)

    for con in condl:
        matching = []
        recon = re.compile(f'^{con}')
        for i in range(32, len(coln)):
            if re.match(recon, coln[i]):
                matching.append(i)
        
        outd[con] = matching 

    return outd

def parse_input_all(infil):
    '''
    infil is a string, path to intput file
    This file contains all metabolites, not just DA ones
    RETURNS: dict, with condition (string) as key, and nested 
    list as values. In each sublist is: peak rt, mz and average area
    If the average area is 0 for a condition, it is not included
    '''
    outd = {}
    with open(infil, 'r') as inf:
        inl = inf.readline()
        conditions = get_groups(inl.strip().split('\t'))
        print(conditions)
        while inl:
            if inl[0].isdigit():
                inlst = inl.split('\t')
                rt = inlst[1]; mz = inlst[2]
                for cond, cols in conditions.items():
                    pkareas = [float(inlst[col]) for col in cols]
                    avgarea = sum(pkareas) / len(pkareas)
                    if avgarea > 0:
                        if cond in outd.keys():
                            outd[cond].append([rt, mz, avgarea])
                        else:
                            outd[cond] = [[rt, mz, avgarea]]
            
            inl = inf.readline()
    
    return outd

def get_total_pk_area(ll):
    '''
    ind is a list of nested lists. One value of parse_input_all output.
    the third entry of each sublist is peak area
    RETURNS: float, sum of all peak area
    '''
    sm = 0
    for i in ll:
        sm += float(i[2])

    return sm

def get_shannon(peaks, total):
    '''
    peaks is a list, one value of parse_input_all output
    total is a float, output of get_total_pk_area
    RETURNS: float, the shannon entropy for a condition
    '''
    print(total)
    intlist = []
    for peak in peaks:
        relint = peak[2]/total ## relative intensity
        logg = math.log(relint,2)
        pij = relint*logg
        intlist.append(pij)

    return -sum(intlist) ## Shannon Entropy

def make_output(condd, odir, expn):
    '''
    condd is a dict, output (total) of parse_inpt_all
    odir is a string, path to where output file should go
    expn is a string
    OUTPUTS: Each line of output is a condition
    Outputs cond name, total # of peaks, shannon entropy
    '''
    
    with open(f'{odir}/{expn}_shannon_entropy_per_condition.tab', 'w') as ofil:
        ofil.write('Condition\tTotal_peak_num\tshannon\n')
        for k, v in condd.items():
            totalpk = get_total_pk_area(v)
            shannon = get_shannon(v, totalpk)
            ofil.write(f'{expn}_{k}\t{len(v)}\t{shannon}\n')

def main(infil, odir, expn):
    '''
    all inputs are strings
    '''
    condd = parse_input_all(infil)
    make_output(condd, odir, expn)
    print('Done!')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit('ARGS: 1) full path to input file 2) desired output directory 3) What\n'\
            'experiment type (tissue, sym)')

    main(sys.argv[1], sys.argv[2], sys.argv[3])
