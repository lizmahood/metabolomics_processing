## new script for second part of canopus (first is filtering mgf)
## Takes in mgf and various directories and alignment file
## very similar to old version

import sys, os

def launch_sirius(mgf, odir, mstype):
    '''
    :param mgf: is full path to folder with input ms files
    :param odir: is where sirius output will be
    '''
    command = 'C'

    os.system(command)

def read_alignment(alif):
    '''
    :param alif: path to alignment file
    :returns: dict of metabolites with id as key and ms1 as value
    '''
    odict = {}
    with open(alif) as alifil:
        lin = alifil.readline().strip()
        ct = 0
        while lin:
            linl = lin.split('\t')
            if linl[0][0].isdigit() and linl[31] != '' and '2+' not in linl[4] and '2-' not in linl[4]:
                odict[linl[0]] = [linl[30], linl[31], linl[4], linl[2]]
            ct += 1
            lin = alifil.readline().strip()

    return odict

def make_mgf_for_compound(mgfinpdir, mgfstr, alid, ct):
    '''
    :param mgfinpdir: is a string, path to where mgfs will be output
    :param mgfstr: is a string of lines to output
    :param alid: is a dict, output of read_alignment
    :param ct: int, will be part of directory name 
    output is: a new ms file for the compound in mgf directory
    returns: name of mgf outputted
    '''
    ##making new directory
    mgfl = mgfstr.split('\n')
    scans = mgfl[1].strip().split(': ')[1]

    ##outputting new ms files
    mgfnam = mgfinpdir + '/' + str(ct) + '_scan' + scans + '.ms'
    with open(mgfnam, 'w') as mgfout: 
        mgfout.write('>compound scans' + scans + '\n')
        for i, lin in enumerate(mgfl):
            if 'PRECURSORMZ' in lin:
                mass = lin.split(': ')[1]
                mgfout.write('>parentmass ' + mass + '\n')
            elif 'FORMULA' in lin:
                formula = lin.split(': ')[1]
                mgfout.write('>formula ' + formula + '\n')
            elif 'PRECURSORTYPE' in lin:
                ion = lin.split(': ')[1]
                mgfout.write('>ionization ' + ion + '\n')
            elif 'Num Peaks' in lin:
                mgfout.write('>ms1\n')
                newms1 = alid[scans][0].replace(' ', '\n').replace(':', ' ')
                mgfout.write(newms1 + '\n>collision 40\n')
                numpks = int(lin.split(': ')[1])
                for index in range(i + 1, i + numpks + 1):
                    linel = mgfl[index].split('\t')
                    newline = linel[0] + ' ' + linel[1] + '\n'
                    if 'END' not in newline:
                        mgfout.write(newline)
        mgfout.write('\n')
    return mgfnam

def main(mgfinp, mgfpardir, alifil):

    alid = read_alignment(alifil)

    ct = 1
    with open(mgfinp, 'r') as inp:
        inlin = inp.readline()
        while inlin:
            if inlin.startswith('NAME'):
                ostring = inlin
            elif inlin == '\n':
                ostring += inlin
                scans = ostring.split('\n')[1].strip().split(': ')[1]
                if scans in alid.keys():
                    nam = make_mgf_for_compound(mgfpardir, ostring, alid, ct)
                    ct += 1
                
            else:
                ostring += inlin
            inlin = inp.readline()
        
    #launch_sirius(mgfpardir, sirdir, mstype)
        
if __name__ == '__main__':

    if len(sys.argv) != 4:
        sys.exit('ARGS: 1) full path to input mgf\n '\
                '2) directory for holding individual mgfs (cannot be sirius dir)\n'\
                    '3) Path to alignment file corresponding to input mgf/msp')

    main(sys.argv[1], sys.argv[2], sys.argv[3])   
    print('Done!') 
