import sys, glob
import pandas as pd
import numpy as np
from fisher import FisherExactTest

## args: canopus summary file, canopus performance file, 
## canopus conf score file, dir of DAM per stress

def filter_canopus_performance(perf):
    '''
    :param perf: string, path to canopus' performance file (w f1 per class)
    :returns: list of classes (or superclasses, etc) passing thresholds
    '''
    perfdf = pd.read_table(perf)
    passing = perfdf[(perfdf['f1'] >= 0.66) & (perfdf['count'] >= 50)]
    return passing['name'].tolist()

def parse_canopus(canopus, conf, class_surv):
    '''
    :param canopus: string, path to canopus summary file 
    :param conf: string, path to canopus conf score file
    :param class_surv: list of classes surviving, output of filter_canopus_perf
    :returns: df of the canopus file with usable alignment IDs and conf scores
    '''
    canopusdf = pd.read_table(canopus)
    confdf = pd.read_table(conf)

    ## parse canopus id
    canopusdf['AlignmentID'] = canopusdf['name'].str.split('scans').str[1]

    ## parse conf id
    confdf['AlignmentID'] = confdf['name'].str.split('scan').str[1]

    ## merging 
    out = pd.merge(canopusdf, confdf, how = 'left', on = 'AlignmentID')
    passing = out[(out['prob'] >= 0.5) & (out['class_y'].isin(class_surv))]
    return passing[['AlignmentID', 'class_y', 'superclass', 'prob']]

def make_fisher_input(canopus_surv, dam_id):
    '''
    :param canopus_surv: df, output of parse_canopus
    :param dam_id: list, of alignment IDs da in stress
    :returns: dict, with class name as key and list as value
    '''
    odict = {}

    sdam_id = [str(x) for x in dam_id]
    canopus_stress = canopus_surv[(canopus_surv['AlignmentID'].isin(sdam_id))]
    canopus_notstress = canopus_surv[(~canopus_surv['AlignmentID'].isin(sdam_id))]
    ## getting list and valuecounts of the DAM's classes
    uniq = canopus_stress['class_y'].unique().tolist()
    stress_counts = canopus_stress['class_y'].value_counts()
    notstress_counts = canopus_notstress['class_y'].value_counts()

    ## getting valuecounts per class in stress and not stressed metabolites
    for clas in uniq:
        cls_str = stress_counts[clas]
        nocls_str = canopus_stress.shape[0] - cls_str
        try:
            cls_nostr = notstress_counts[clas]
        except:
            cls_nostr = 0
        nocls_nostr = canopus_notstress.shape[0] - cls_nostr
        odict[clas] = [cls_str, nocls_str, cls_nostr, nocls_nostr]

    return odict

def fisher(tuple_t):# every 4 number tuple
    """
    :param tuple_t: one value in the output of make_fisher_input
    * k number of objects in the selection having a given property
    * n size of the selection
    * C number of objects in the population having this property
    * G size of the population
    """
    t1=int(tuple_t[0])
    t2=int(tuple_t[1])
    t3=int(tuple_t[2])
    t4=int(tuple_t[3])

    f = FisherExactTest()
    p=f.pvalue(t1,t2+t1,t1+t3,t1+t2+t4+t3)
    en=f.enrichment(t1,t2+t1,t1+t3,t1+t2+t4+t3)

    if en >= 1:
        en="+" #using right hand side Ha :u
    else:
        en="-"

    pv=p[2]
    return en,pv

def main(perf, canopus, conf, dam_dir, outn):
    '''
    all inputs strings
    '''
    passing_classes = filter_canopus_performance(perf)
    passing_metabs = parse_canopus(canopus, conf, passing_classes)

    dir_to_parse = f'{dam_dir}*FDR*'
    stress_dams = glob.glob(dir_to_parse)

    for dam in stress_dams:
        damn = dam.split('_')[-1]
        print(damn)
        damdf = pd.read_table(dam)
        ids = damdf['Alignment.ID'].tolist()
        stressd = make_fisher_input(passing_metabs, ids)
        print(ids)
        print(stressd)

        if stressd:
            for clas, tup in stressd.items():
                pv, en = fisher(tup)
                tup.append(pv)
                tup.append(en)

            ## make to df, sort based on p-value, write out
            out = pd.DataFrame.from_dict(stressd, orient = 'index')
            out.columns = ['NumClassDAM', 'NumNotClassDAM', 'NumClassNotDAM', 'NumNotClassNotDAM','En', 'pValue']
            out = out.sort_values('pValue')
            out.to_csv(f'{outn}_{damn}', sep = '\t')
        else:
            print(f'{damn} had no Canopus-classified metabolites!')

if __name__ == '__main__':

    if len(sys.argv) != 6:
        sys.exit('ARGS: 1) Canopus performance file (from Canopus authors) '\
            '2) Canopus summary file 3) Canopus conf score file 4) Directory of DAM files (with "FDR" in the name)'\
                '5) output path and full name')

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

    print('Done!')