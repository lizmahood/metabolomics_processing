##Either: takes in complete mgf and launches sirius on that
##Or takes in directory of mgfs and launches sirius on each individual mgf
##Second option is identical to 2make_and_launch_compounds_individually.py in flavo_vs_antho

#1) open input mgf
#2) write out new mgf when the big input file gets to "BEGIN IONS". Write this in its own dir
#3) stop writing out when we get to "END IONS"
#4) make fstring of command for sirius
#5) launch command. 

import sys, os

def make_mgf_for_compound(sirpardir, mgfinpdir, mgfstr, ct):
    '''
    @param: sirpardir is a string, path to parent dir of sirius outputs
    @param: mgfinpdir is a string, path to where mgfs will be output
    @param: mgfstr is a string of lines to output
    @param: ct is an int - counter to include in directory name
    output is: a new mgf for the compound in its own directory
    returns: name of mgf outputted
    '''
    ##making new directory
    print(mgfstr)
    mgfl = mgfstr.split('\n')
    scans = mgfl[1].strip().split('=')[1]
    newdir = f'{sirpardir}/{ct}_scan{scans}/'
    os.mkdir(newdir)

    ##outputting new mgf
    mgfnam = f'{mgfinpdir}/{ct}_scan{scans}.mgf'
    with open(mgfnam, 'w') as mgfout: 
        mgfout.write(mgfstr + '\n')
    
    return mgfnam, newdir

def launch_sirius(mgf, odir):
    '''
    @param: mgf is full path to input mgf
    @param: odir is where sirius output will be
    '''
    command = f'sirius -o {odir} -i {mgf} formula -c 3 -p orbitrap '\
        '--ppm-max-ms2 5.0 --elements-enforced HCNOPS '\
            'structure --database BIO canopus '\

    #   config --IsotopeMs2Settings IGNORE --FormulaSearchDB --FormulaSettings.detectable , --Timeout.secondsPerTree 0 --Timeout.secondsPerInstance 0 , --StructureSearchDB BIO --AdductSettings.fallback [[M + H]+] --FormulaResultThreshold true --RecomputeResults true formula structure canopus
    os.system(command)

def main_individually(mgfinp, sirpardir, mgfpardir):

    ct = 1
    with open(mgfinp, 'r') as inp:
        inlin = inp.readline()
        while inlin:
            if inlin.startswith('BEGIN'):
                ostring = inlin
            elif inlin.startswith('END'):
                ostring += inlin
                inlin = inp.readline()
                ostring += inlin
                nam, sirdir = make_mgf_for_compound(sirpardir, mgfpardir, ostring, ct)
                ct += 1
                launch_sirius(nam, sirdir)
            else:
                ostring += inlin
            inlin = inp.readline()
        
def main_together(mgfinp, sirpardir):

    launch_sirius(mgfinp, sirpardir)

if __name__ == '__main__':

    if len(sys.argv) != 5:
        sys.exit('ARGS: 1) full path to input mgf '\
            '2) sirius dir, in which to output sirius result directories '\
                '3) directory for holding individual mgfs (cannot be sirius dir) '\
                    '4) are mgfs launched together OR separate ?')

    if sys.argv[4] == 'separate':
        main_individually(sys.argv[1], sys.argv[2], sys.argv[3])   
    elif sys.argv[4] == 'together':
        main_together(sys.argv[1], sys.argv[2])
    print('Done!') 
