## Takes input folder, output folder, and Parameter file as inputs (all full paths)
## launches MsdialConsoleApp.exe on that data
## Takes the alignresult file that was made in the output folder
## Reads it in, finds 34th and 35th column, stores them as list
## Finds correlation between 2 lists, % of first list that is 0 and % of second list that is 0
## Just prints out this info

import sys, os, glob
from scipy.stats.stats import pearsonr

def launch_msdial(input, output, params):
    '''
    '''
    os.system(f'C:/LCMS_Software/MSDIALver4.48/MSDIALver4.48/MsdialConsoleApp.exe '\
    f'lcmsdda -i {input} -o {output} -m {params}')

def parse_msdial_output(output):
    '''
    '''
    ofil = glob.glob(f'{output}/AlignResult*')[0]
    btil = []
    jgil = []
    with open(ofil) as alif:
        lin = alif.readline()
        while lin:
            linl = lin.split('\t')
            if linl[0].startswith('Align'):
                btinam = linl[32]
                jginam = linl[33]
            try:
                mz = float(linl[2])
                if mz < 800 and mz > 100:
                    btil.append(float(linl[32]))
                    jgil.append(float(linl[33]))
            except:
                pass
            lin = alif.readline()

    corr = pearsonr(btil, jgil)
    btimissing = (btil.count(0) / len(btil)) * 100
    jgimissing = (jgil.count(0) / len(jgil)) * 100

    print(f'BTI Sample: {btinam}\nJGI Sample: {jginam}\n'\
        f'Correlation: {corr}\nBTI Missing: {btimissing}%\n'\
            f'JGI Missing: {jgimissing}%')

def main(input, output, params):
    '''
    '''
    launch_msdial(input, output, params)
    parse_msdial_output(output)

if __name__ == '__main__':

    if len(sys.argv) != 4:
        sys.exit('ARGS: 1) Input folder to put into MSDIAL 2) Output folder for MSDIAL '\
            '3) Parameter file for MSDIAL (full path)')

    main(sys.argv[1], sys.argv[2], sys.argv[3])

    print('Done!')
