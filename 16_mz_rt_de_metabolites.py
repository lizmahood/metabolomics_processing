import sys
import seaborn as sns
import matplotlib.pyplot as plt

def parse_alignment_file(infil, normtype):
    '''
    infil is a string
    normtype is a string, is OR vsn
    RETURNS: RTs and m/z of Differentially abundant metabolites
    as 2 lists
    '''
    mzl = []; rtl = []
    with open(infil, 'r') as inf:
        inl = inf.readline()
        inlst = inl.strip().split('\t')
        rtc = inlst.index('Average.Rt.min.')
        mzc = inlst.index('Average.Mz')
        pvalc = inlst.index('adj_pval')
        if normtype == 'vsn':
            fcc = inlst.index('ofc')
        elif normtype == 'is':
            fcc = inlst.index('nfold_change')
        
        while inl:
            if inl[0].isdigit():
                inlst = inl.strip().split('\t')
                try: 
                    float(inlst[fcc])
                    if float(inlst[pvalc]) <= 0.05 and abs(float(inlst[fcc])) >= 2:
                        mzl.append(float(inlst[mzc]))
                        rtl.append(float(inlst[rtc]))
                except: 
                    pass
            inl = inf.readline()

    return mzl, rtl

def make_hists(mzl, rtl, odir):
    '''
    mzl and rtl are lists, output of parse_alignment_file
    OUTPUTS: 2 histograms
    '''
    mzlen = len(mzl); rtlen = len(rtl)
    f, axes = plt.subplots(2, 1, figsize=(9,6), sharex=False)
    sns.distplot(mzl, color = 'red', kde = True,label = 'DAM, n=' + str(mzlen), ax = axes[0])

    axes[0].legend()
    axes[0].set_title('m/z')
    sns.distplot(rtl, color = 'red', kde = False, bins = 16,label = '', ax = axes[1])

    #axes[1].legend()
    axes[1].set_title('Retention Time')
    plt.savefig(odir + 'rt_mz_hist.pdf', bbox_inches = 'tight')

def main(infil, normtype, odir):
    '''
    all inputs strings
    '''
    mzl, rtl = parse_alignment_file(infil, normtype)
    make_hists(mzl, rtl, odir)
    print('Done!')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit('ARGS: 1) input file, 2) normalization type (is OR vsn)\n'\
            '3) desired output directory')

    main(sys.argv[1], sys.argv[2], sys.argv[3])
