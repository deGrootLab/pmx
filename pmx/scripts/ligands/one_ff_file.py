import sys, os
#from pymacs_beta import *
#from pymacs_beta.parser import *
#from pymacs_beta.forcefield import *
#from pymacs_beta.ndx import *
from pmx import *
from pmx.parser import *
from pmx.forcefield import *
from pmx.ndx import *


class FFatom:
    def __init__(self,list):
        self.type = list[0]
        self.sigmaA = list[1]
        self.epsA = list[2]
        self.A = list[3]
        self.sigmaB = list[4]
        self.epsB = list[5]

class FFfile:
    def __init__(self, fname=None):
        if fname is not None:
            foo = fname.split('.')
            bar = re.sub('ff','',foo[0])
            self.name = bar
            self.atoms = []
            self.read_ffitp(fname)

    def read_ffitp(self,file):
        l = open(file).readlines()
        toSkip = "atomtypes"
        for line in l:
            if toSkip not in line:
                self.atoms.append(FFatom(line.split()))

def get_atoms( ffs ):
    atoms = {}
    for ffile in ffs:
        ff = FFfile(ffile)
        for at1 in ff.atoms:
            if at1.type in list(atoms.keys()): # check if atom type already exists
                at2 = atoms[at1.type]
                if (at1.type == at2.type and at1.sigmaA == at2.sigmaA and at1.epsA == at2.epsA and at1.sigmaB == at2.sigmaB and at1.epsB == at2.epsB):
                    continue
                else:
                    sys.stdout.write('Found two atoms of type %s, but they have different parameters, consider renaming atom types\n' % at1.type)
            else:
                atoms[at1.type] = at1

    return (atoms)


def ffwrite(atoms,file):
    fp = open(file,'w')
    print('[ atomtypes ]', file=fp)
    for atype in list(atoms.keys()):
        at = atoms[atype]
        print('   %s      %s      %s      %s      %s      %s' % (at.type,at.sigmaA,at.epsA,at.A,at.sigmaB,at.epsB), file=fp)

def main(argv):

        # define input/output files
    files= [
       FileOption("-ffitp", "r/m",["itp"],"ffMOL.itp",""),
       FileOption("-ffitp_out", "w",["itp"],"ffMOL_out.itp",""),
    ]
    # define options

    options=[]

    help_text = ()

    # pass options, files and the command line to pymacs

    cmdl = Commandline( argv, options = options,
                   fileoptions = files,
                   program_desc = help_text,
                   check_for_existing_files = False )

    ffs = cmdl['-ffitp']

    #get 2 lists of: 1) equivalent atoms in ffMOL1.itp and ffMOL2.itp 2) different atoms
    atoms = get_atoms(ffs)

    #write new ffMOL.itp
    ffwrite(atoms,cmdl['-ffitp_out'])

    sys.stdout.write('All finished correctly\n')

main( sys.argv )
