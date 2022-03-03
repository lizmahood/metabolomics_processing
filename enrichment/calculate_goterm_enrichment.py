import numpy as np
import pandas as pd
from fisher import FisherExactTest
from statsmodels.stats.multitest import multipletests
import sys

def get_df_good(inpt):
    '''
    :param inpt: string, path to input file
    :returns: df
    '''
    df = pd.read_table(inpt, sep = ' ')
    if df.shape[1] != 5:
        df = pd.read_table(inpt, sep = '\t')
    print(df)
    return df

def get_enrichment(df):
    '''
    :param df: input df, with cols as NAME-t1-t2-t3-t4
    :returns: df as above but also with enrichment and pvalue columns
    '''
    fdict = df.to_dict(orient = 'index')
    pvall, enl = [], []
    f = FisherExactTest()
    for k, v in fdict.items():
        t1 = v['In_stress_with_GOterms']
        t2 = v['In_stress_without_GOterms']
        t3 = v['Exclusion_with_GOterms']
        t4 = v['Exclusion_without_GOterms']
        pval = f.pvalue(t1, t1 + t2, t1 + t3, t1 + t2 + t3 + t4)
        en = f.enrichment(t1, t1 + t2, t1 + t3, t1 + t2 + t3 + t4)
        if en >= 1:
            en="+" #using right hand side Ha :u
        else:
            en="-"
        enl.append(en)
        pv=pval[2]
        pvall.append(pv)

    padj = multipletests(pvall, method = 'fdr_bh')
    outdf = df
    print(outdf.shape)
    outdf['enrichment'] = enl    
    outdf['pvalue'] = pvall
    outdf['qvalue'] = padj[1]
    outdf = outdf.sort_values('qvalue')

    return outdf

def main(inpt):
    '''
    inpt is string
    '''

    indf = get_df_good(inpt)
    outdf = get_enrichment(indf)
    oname = '.'.join(inpt.split('.')[:-1])
    print(oname)
    outdf.to_csv(f'{oname}_enrichment.tsv', sep = '\t', index = False)

if __name__ == '__main__':

    if len(sys.argv) != 2:
        sys.exit('ARGS: 1) input file of: GOTERM-t1-t2-t3-t4')

    main(sys.argv[1])

    print('Done!')