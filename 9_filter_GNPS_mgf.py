import sys

def get_surviving_ids(fil):
    '''
    fil is a string, path to file with IDs
    RETURNS: list of ids that are in the file 
    and have ms/ms fragments
    '''
    outl = []

    with open(fil, 'r') as msfil:
        msl = msfil.readline()
        while msl:
            mslist = msl.split('\t')
            if mslist[31] != '':
                outl.append(mslist[0])
            msl = msfil.readline()
        
    return outl

def filter_mgf(mgff, msl):
    '''
    mgff is a string, path to mgf
    msl is a list, output of get_surviving_ids
    OUTPUTS: filtered mgf with only 
    instances that survived filtering
    '''
    mgf = open(mgff, 'r')
    filtmgf = open(mgff + '_only_surviving.mgf', 'w')
    mgfln = mgf.readline()
    while mgfln:
        if mgfln.startswith('BEGIN'):
            mgfln = mgf.readline()
            scn = mgfln.split('=')[1]
            if scn.strip() in msl:
                filtmgf.write(f'BEGIN IONS\nSCANS={scn}')
                while not mgfln.strip() == '':
                    mgfln = mgf.readline()
                    filtmgf.write(mgfln)
        mgfln = mgf.readline()

    mgf.close()
    filtmgf.close()

def main():
    if len(sys.argv) != 3:
        sys.exit('ARGS: 1) mgf to filter 2) filtered PeakArea file')

    survl = get_surviving_ids(sys.argv[2])
    filter_mgf(sys.argv[1], survl)

    print('Done!')

if __name__ == '__main__':
    main()


    