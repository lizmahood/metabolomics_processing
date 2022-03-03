import sys, os
import pandas as pd

def get_mgf_id(mgfin):
    '''
    mgfin is a string
    returns: list of ints
    '''
    outl = []
    with open(mgfin, 'r') as mgff: 
        mgfl = mgff.readline()
        while mgfl:
            if 'SCANS' in mgfl:
                idd = int(mgfl.strip().split('=')[1])
                outl.append(idd)
            mgfl = mgff.readline()

    return outl 

def main(inpdfp, inpmgf, ofil):
    '''
    :param inpdfp: string, path to filtered, normalized peak area data frame
    :param inpmgf: string, path tofiltered mgf (NO NOISE FILTER one)
    '''

    idl = get_mgf_id(inpmgf)
    inpdf = pd.read_csv(inpdfp, sep = '\t')
    outdf = inpdf[inpdf['Alignment.ID'].isin(idl)]
    outdf = outdf.replace({0.01: 0})
    outdf.to_csv(ofil, sep = '\t', index = False)

if __name__ == '__main__':

    if len(sys.argv) != 4:
        sys.exit('ARGS: 1) Path to peak Area data frame, normalized and filtered,'\
         '2) Path to filtered mgf (NO NOISE FILTER one) 3) Output file name')

    main(sys.argv[1], sys.argv[2], sys.argv[3])
        
    print('Done!')
