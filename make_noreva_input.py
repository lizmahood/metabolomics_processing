import sys, re, pandas as pd
import numpy as np

'''
read in filtered file with pandas 
make list of new column names ('C' + alignmentID)
get columns that are experimental
make vector of these column names
remove columns that are not these
transpose

read in non-filtered file
get rows of internal standards
get values of these in experimental
columns.
Transpose these, call them IS1 and IS2
Add to df before experimental metabolites
write out
'''
def argparser(argl):
    '''
    argl is a list, returns tuple of args
    '''
    for ar in range(1, len(argl)):
        if argl[ar] == '-ntype':
            ntype = argl[ar+1]
        elif argl[ar] == '-is1':
            is1 = argl[ar+1]
        elif argl[ar] == '-is2':
            is2 = argl[ar+1]
        elif argl[ar] == '-inpfilt':
            inpfilt = argl[ar+1]
        elif argl[ar] == '-inpnonfilt':
            inpnonfilt = argl[ar+1]
    
    if ntype == 'is':
        return(ntype, inpfilt, is1, is2, inpnonfilt)
    else:
        return(ntype, inpfilt)


def change_alignment_names(value):
    '''
    value is a int, ntype is a string
    '''
    return 'C' + str(value)

def get_labels(sampnames):
    '''
    sampnames is a list of strings
    RETURNS: another list of strings
    '''
    return [x.split('_')[0] for x in sampnames]

def get_qc_labels(sampnames):
    '''
    sampnames is a list of strings
    RETURNS: another list of strings
    '''
    olst = []
    for i, nam in enumerate(sampnames):
        nn = nam.replace('.Moghe', '')
        nns = nn.split('_')[0]
        olst.append(f'batch01_{nns}{i+1}')

    return olst


def main(ntype, inpfilt, **kwargs):

    ##getting metabolite values
    filt = pd.read_table(inpfilt)
    filt['Alignment.ID'] = filt['Alignment.ID'].apply(change_alignment_names)
    goodcols = list(filt.columns[32:-1])

    if ntype == 'is':
        metabdf = filt.iloc[:,32:-1].T  
        metabdf.columns = filt['Alignment.ID']
        
        ##getting IS rows
        is1 = kwargs['is1']; is2 = kwargs['is2']
        inpnonfilt = kwargs['inpnonfilt']
        nonfilt = pd.read_table(inpnonfilt)
        newcols = [x.replace(' ', '.').replace('-', '.') for x in nonfilt.columns]
        nonfilt.columns = newcols
        is1row = nonfilt.loc[nonfilt['Alignment.ID'] == int(is1)][goodcols].T
        is2row = nonfilt.loc[nonfilt['Alignment.ID'] == int(is2)][goodcols].T
        is1row.columns = ['IS1']
        is2row.columns = ['IS2']

        ##putting it all together
        label = get_labels(goodcols)
        out = pd.concat([is1row, is2row], axis=1)
        out['SampleName'] = goodcols
        out['Label'] = label
        neworder = ['SampleName', 'Label', 'IS1', 'IS2']
        out = out[neworder]
        finalout = pd.concat([out, metabdf], axis = 1)

    elif ntype == 'qc':
        ##need to change around column order based on run
        goodcols.sort(key = lambda x: int(x.rsplit('_', 1)[1]))
        metabs = filt.iloc[:,32:-1]
        metabs = metabs[goodcols].T
        metabs.columns = filt['Alignment.ID']

        ##making final columns
        sampn = get_qc_labels(goodcols)
        print(sampn)
        metabs.insert(0, 'sample', sampn)
        metabs.insert(1, 'batch', np.ones((len(sampn),), dtype = int))
        Label = get_labels(goodcols)
        metabs.insert(2, 'label', Label)
        metabs.loc[metabs['label'] == 'QC.Moghe', 'label'] = 'NA'
        metabs.insert(3, 'order', np.arange(1, len(sampn)+1))
        finalout = metabs
        print(finalout.iloc[:, 0:5])

    ##writing out
    finalout.to_csv(inpfilt + 'for_noreva.csv', index = False)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.exit('-ntype: is OR qc\n-is1: if ntype = is, what is the is1 Alignment ID?\n'\
            '-is2: what is the is2 Alignment ID?\n-inpfilt: path to filtered input file\n'\
                '-inpnonfilt: if is, what is the path to the original PeakArea file')

    argg = argparser(sys.argv)
    ntype = argg[0]

    if ntype == 'is':
        main(ntype = ntype, inpfilt = argg[1], is1 = argg[2], is2 = argg[3], inpnonfilt = argg[4])
    elif ntype == 'qc':
        main(ntype = ntype, inpfilt = argg[1])
    
    print('Done!')



