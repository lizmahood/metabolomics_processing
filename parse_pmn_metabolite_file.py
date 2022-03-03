import sys
import pandas as pd
pd.options.mode.chained_assignment = None 
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import itertools

def read_compounds(compdat):
    '''
    :param compdat: string, path to compounds.dat
    :returns: dict of mono-isotopic masses as key and ID as values
    '''
    odict = {}
    with open(compdat) as cd:
        inlist = cd.readlines()
        for inl in inlist:
            if inl.startswith('UNIQUE'):
                cid = inl.strip().split(' - ')[1]
                tmp = [cid]
            elif inl.startswith('TYPES'):
                typ = inl.strip().split(' - ')[1]
                tmp.append(typ)
            elif inl.startswith('COMMON'):
                name = inl.strip().split(' - ')[1]
                name.replace(',', '').replace('\'', '').replace(' ', '')
                tmp.append(name)
            elif inl.startswith('MONOISO'):
                mass = round(float(inl.strip().split(' - ')[1]), 4)
                tmp.append(mass)
            elif inl.startswith('//'):
                if mass not in odict.keys():
                    odict[mass] = [cid, typ, name]
                else:
                    odict[mass][0] += '|' + cid
                    odict[mass][1] += '|' + typ
                    odict[mass][2] += '|' + name

    return odict

def get_masses_of_compounds(compd, ion_mode):
    '''
    :param compd: dict, output of read_compounds
    :param ion_mode: string, pos OR neg
    :returns: dict of dicts, keys of outer are plantcyc
    masses, values of outer is dict of adduct masses
    '''
    if ion_mode == 'pos':
        adducts = ['[M+H]+', '[M+NH4]+', '[M+ACN+H]+']
        adduct_mass = [1.007276, 18.033823, 42.033823]
    elif ion_mode == 'neg':
        adducts = ['[M-H]-', '[M-H2O-H]-', '[M+FA-H]-']
        adduct_mass = [-1.007825, -19.0183897, 44.997654]

    odict = {}
    for k in compd.keys():
        odict[k] = {adducts[0]: k+adduct_mass[0], 
        adducts[1]: k+adduct_mass[1], adducts[2]: k+adduct_mass[2]}

    return odict

def read_mgf(mgfin):
    '''
    :param mgfin: string, path to our mgf
    :returns: dict, with masses as keys and IDs as values
    '''
    odict = {}
    with open(mgfin) as mgf:
        mgfl = mgf.readline()
        while mgfl:
            if mgfl.startswith('PEP'):
                mass = float(mgfl.strip().split('=')[1])
            elif mgfl.startswith('SCANS'):
                idd = mgfl.strip().split('=')[1]
            elif mgfl.startswith('END'):
                if mass in odict.keys():
                    odict[mass] += '|' + idd
                else:
                    odict[mass] = idd
            mgfl = mgf.readline()
    return odict

def match_mass_to_compounds(compd, massd, mgfd, thresh):
    '''
    :param compd: dict, output of read_commpounds
    :param massd: dict of adduct masses for each compound MIM
    :param mgfd: dict, output of read_mgf
    :param thresh: float, what mass error is allowed for two compounds to be the same?
    :returns: 3 col df (mass, our compounds, plant cyc compounds) and a list
    '''
    outdf = pd.DataFrame(columns = ['our_mass', 'our_id',
     'plantcyc_mass', 'plantcyc_cid', 'plantcyc_type', 'plantcyc_name', 'adduct'])
    compdf = pd.DataFrame.from_dict(mgfd, orient = 'index')
    compdf.reset_index(inplace=True)
    compdf.columns = ['our_mass', 'our_id']
    for k, v in compd.items():
        masses = massd[k]
        for adduct, mass in masses.items():
            tmp = compdf[(abs(compdf['our_mass'] - mass) <= thresh)]
            try:
                if tmp.shape[1] > 1:
                    tmp['plantcyc_mass'] = [mass] * tmp.shape[0]
                    tmp['plantcyc_cid'] = [v[0]] * tmp.shape[0]
                    tmp['plantcyc_type'] = [v[1]] * tmp.shape[0]
                    tmp['plantcyc_name'] = [v[2]] * tmp.shape[0]
                    tmp['adduct'] = [adduct] * tmp.shape[0]
                    outdf = outdf.append(tmp, ignore_index = True)
            except: 
                pass
    return outdf

def get_matches_for_our_metabs(mgfd, matchdf):
    '''
    :param mgfd: dict of our mgf
    :param matchdf: df of matches between plantcyc and us
    '''
    mgfdf = pd.DataFrame.from_dict(mgfd, orient = 'index')
    mgfdf.reset_index(inplace = True)
    mgfdf.columns = ['our_mass', 'our_id']

    mgfdf['All_ids'] = mgfdf.our_id.str.split('|')
    unique_ids = set(itertools.chain.from_iterable(mgfdf.All_ids))
    olist = []
    for cid in unique_ids:
        mtch = matchdf.our_id.str.count(cid).sum()
        olist.append(mtch)

    return olist

def plot_num_matches(matches, ofil, thresh):
    '''
    :param matches: list of ints
    :param ofil: string, path to output file
    :param thresh: float, mass error threshold
    '''
    counts = Counter(matches)
    keys = list(counts.keys())
    values = list(counts.values())
    both = list(zip(keys, values))
    toplot = sorted(both, key = lambda x: x[0])
    res = list(zip(*toplot))
    df = pd.DataFrame(res).transpose()
    df.columns = ['x', 'y']
    myplot = sns.barplot(x = 'x', y = 'y', data = df).set_title(f'Number of Matches in '\
        f'our Data to PlantCyc \nMetabolites with Threshold {thresh}')
    fig = myplot.get_figure()
    fig.savefig(ofil)

def main(comp, mgfin, ofil, thresh, ion_mode):
    '''
    All inputs strings
    '''
    compd = read_compounds(comp)
    mgfd = read_mgf(mgfin)
    massd = get_masses_of_compounds(compd, ion_mode)
    matchdf = match_mass_to_compounds(compd, massd, mgfd, thresh)
    numl = get_matches_for_our_metabs(mgfd, matchdf)
    matchdf.to_csv(f'{ofil}_matches{thresh}.tsv', sep = '\t')
    plot_num_matches(numl, f'{ofil}_match_numbers{thresh}.pdf', thresh)

if __name__ == '__main__':
    
    if len(sys.argv) != 6:
        sys.exit('ARGS: 1) compounds.dat 2) our mgf 3) output path and prefix '\
        '4) mass error threshold 5) ion mode, pos OR neg')

    main(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5])
    print('Done!')
