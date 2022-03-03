import sys, os

def get_all_metabs(clustfil, impclust):
    '''
    clustfil is a string, path to cluster file input
    impclust is an int, id of importnat clusters
    OUTPUTS: list of metabolite ids to make msps from
    '''
    olist = []
    print(impclust)
    with open(clustfil, 'r') as fil:
        lin = fil.readline(); lin = fil.readline()
        while lin:
            if not lin.startswith('all'):
                ll = lin.strip().split('\t')
                if ll[1] == str(impclust):
                    olist.append(ll[0])
            lin = fil.readline()

    return olist

def parse_mgf(mgfin):
    '''
    mgfin is a string, path to input mgf
    RETURNS: dict with scan as key and rest as values
    '''
    odict = {}
    with open(mgfin, 'r') as mgff:
        lin = mgff.readline()
        while lin:
            if lin.startswith('SCAN'):
                name = lin.strip().split('=')[1]
                odict[name] = []
                peaks = []
            elif lin.startswith('PEP'):
                pm = float(lin.strip().split('=')[1])
                odict[name].append(pm)
            elif lin.startswith('RT'):
                odict[name].append(lin.strip().split('=')[1])
            elif lin.startswith('ION'):
                odict[name].append(lin.strip().split('=')[1])
            elif lin[0].isdigit():
                peaks.append(lin)
            elif lin.startswith('END'):
                odict[name].append(peaks)
            lin = mgff.readline()
    
    return odict

def write_out_msps(odir, mgfd, clusts, mode):
    '''
    odir is a string, path to output dir
    mgfd is a dict, output of above func
    clusts is a list of strings
    mode is a string, "Positive" or Negative
    OUTPUTS: a msp file for each alignment id in clusts
    '''
    for k, v in mgfd.items():
        if k in clusts:
            with open(f'{odir}/Cmp{k}.msp', 'w') as ofil:
                ofil.write(f'Name:{k}\n')
                ofil.write(f'PRECURSORMZ:{v[0]}\n')
                ofil.write(f'RETENTIONTIME:{v[1]}\n')
                ofil.write(f'IONMODE:{mode}\n')
                ofil.write(f'PRECURSORTYPE:{v[2]}\n')
                ofil.write(f'Num Peaks:{len(v[-1])}\n')
                for pk in v[-1]:
                    ofil.write(f'{pk}')

    print('Done!')

def main(clustfil, impclust, mgfin, mode, odir):
    '''
    all inputs are strings
    '''
    clustl = get_all_metabs(clustfil, impclust)
    mgfd = parse_mgf(mgfin)
    write_out_msps(odir, mgfd, clustl, mode)

if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.exit('ARGS: 1) input to cluster file (_node.txt file)'\
            '2) What cluster do you want to send to MSFINDER?'\
                '3) path to mgf input'\
                    '4) Mode: "Positive" or "Negative"'\
                        '5) Desired output directory to write msp files to')

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])