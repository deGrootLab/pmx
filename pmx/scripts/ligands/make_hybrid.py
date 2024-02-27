import sys, os
from pmx import *
from pmx.parser import *
from pmx.forcefield2 import *
from pmx.ndx import *
import copy as xcopy

def check_if_vsite2_exists( vsites, xs):
    for v in vsites:
        if v[0].id==xs[0].id and v[1].id==xs[1].id and v[2].id==xs[2].id and v[3]==xs[3] and v[4]==xs[4]:
            return(True)
        elif v[0].id==xs[0].id and v[1].id==xs[1].id and v[2].id==xs[2].id:
            print("ERROR: it seems that vsites are to be morphed. Gromacs does not support that")
            sys.exit(0)
    return(False)

def check_if_vsite3_exists( vsites, xs):
    for v in vsites:
        if v[0].id==xs[0].id and v[1].id==xs[1].id and v[2].id==xs[2].id and v[3].id==xs[3].id and v[4]==xs[4] and v[5]==xs[5]:
            return(True)
        elif v[0].id==xs[0].id and v[1].id==xs[1].id and v[2].id==xs[2].id and v[3].id==xs[3].id:
            print("ERROR: it seems that vsites are to be morphed. Gromacs does not support that")
            sys.exit(0)
    return(False)

def check_if_exclusion_valid( exclusions, ex ):
    # return True, if all fine and ex can be added

    # check if the first atom is dummy. If so, all fine
    a0 = ex[0]
    if ('DUM' in a0.atomtype) or ('DUM' in a0.atomtypeB):
        return(True)

    # check if an exclusion already exists
    for excl in exclusions:
        if excl[0]==ex[0]: # some interactions for this atom are already excluded
            if len(excl)!=len(ex): # the exclusion lists are of different length
                print(('ERROR: Something wrong with exclusions: ',excl,ex))
                sys.exit(1)
            else:
                for e in excl[1:]:
                    if e not in ex:
                        print(('ERROR: Something wrong with exclusions: ',excl,ex))
                        sys.exit(1)
    return(True)


def identify_ring_atoms( fname, num ):
    import rdkit
    from rdkit import Chem
    pdbName,atomNameDict = reformatPDB(fname,num)
    mol = Chem.MolFromPDBFile(pdbName,removeHs=False,sanitize=False)
    os.remove(pdbName)
    out = []
    for a in mol.GetAtoms():
        if a.IsInRing()==True:
            out.append(a.GetIdx()+1)
    return(out)

def __atoms_morphe( atoms ):
    for atom in atoms:
        if atom.atomtypeB is not None and (atom.q!=atom.qB or atom.m != atom.mB): return True
    return False

def sum_charge_of_states( itp ):
    qA = 0.0
    qB = 0.0
    for a in itp.atoms:
        qA += a.q
        if __atoms_morphe([a]):
            qB+=a.qB
        else:
            qB+=a.q
    return [qA], [qB]

def findIDinList(ind,decoupAngles):
    for i in decoupAngles:
        if( i==ind ):
            return True
    return False

def adjustCoords(m,mol):
    conf = mol.GetConformer()
    for ai in m.atoms:
        ind = ai.id
#       print ai.x[0]
        posj = conf.GetAtomPosition(ind-1)
        ai.x[0] = posj.x
        ai.x[1] = posj.y
        ai.x[2] = posj.z

def reformatPDB(filename,num,bStrict=False):
#    newname = filename.split(".")[0]+"_tempFormat.pdb"
    newname = "tempFormat"+str(num)+".pdb"
    m = Model().read(filename)

    # adjust atom names and remember the changes
    atomNameDict = {}
    sigmaHoleCounter = 1
    for a in m.atoms:
        newAtomName = a.name
        if 'EP' in a.name:
            newAtomName = 'HSH'+str(sigmaHoleCounter)
            sigmaHoleCounter+=1
        atomNameDict[newAtomName] = a.name
        a.name = newAtomName

    writeFormatPDB(newname,m,bStrict=bStrict)
#    m.write(newname)
    return(newname,atomNameDict)

def restoreAtomNames(mol,atomNameDict):
    for atom in mol.GetAtoms():
        newname = atom.GetMonomerInfo().GetName()
        if newname in list(atomNameDict.keys()):
            oldname = atomNameDict[newname]
            atom.GetMonomerInfo().SetName(oldname)

def writeFormatPDB(fname,m,title="",nr=1,bStrict=False):
    atNameLen = 4
    if bStrict:
        atNameLen = 3
    fp = open(fname,'w')
    for atom in m.atoms:
        foo = xcopy.deepcopy(atom)
        # chlorine
        if( 'CL' in atom.name or 'Cl' in atom.name or 'cl' in atom.name ):
            foo.name = "CL"+"  "
            print(foo, file=fp)
        # bromine
        elif( 'BR' in atom.name or 'Br' in atom.name or 'br' in atom.name ):
            foo.name = "BR"+"  "
            print(foo, file=fp)
        elif( len(atom.name) > atNameLen): # too long atom name
            foo = xcopy.deepcopy(atom)
            foo.name = foo.name[:atNameLen]
            print(foo, file=fp)
        else:
            print(atom, file=fp)
    print('ENDMDL', file=fp)
    fp.close()


def do_log(fp, s):
    l = "make_hybrid__log_> "+s
    print(l, file=sys.stderr)
    print(l, file=fp)


def make_pairs(m1, m2, m3, m4, bFit, bDist, dd, plist = None, grps = None):

    pairs = []
    if plist:
        for n1, n2 in plist:
            a1 = m1.fetch_atoms(n1,how='byid')[0]
            #print a1
            a2 = m2.fetch_atoms(n2,how='byid')[0]
            pairs.append( (a1, a2))
            for atom3 in m3.atoms:
                if a1.id == atom3.id:
                    atom3.x = a2.x
                    if(bFit==True):
                        a4 = m4.fetch_atoms(n2,how='byid')[0]
                        atom3.x = a4.x
        return pairs
    if(grps and bDist):
        lst1 = m1.get_by_id(grps[0])
        lst2 = m2.get_by_id(grps[1])
        for atom in lst1:
            mi = dd # nm
            keep = None
            for at in lst2:
                d = (atom-at)/10.0
                if d < mi:
                    keep = at
                    mi = d
            if keep is not None:
                pairs.append( (atom, keep) )
                for atom3 in m3.atoms:
                    if atom.id == atom3.id:
                        atom3.x = keep.x
        return pairs

    if(bDist):
        for atom in m1.atoms:
            mi = dd # nm
            keep = None
            for at in m2.atoms:
                d = (atom-at)/10.0
                if d < mi:
                    keep = at
                    mi = d
            if keep is not None:
                pairs.append( (atom, keep) )
                for atom3 in m3.atoms:
                    if atom.id == atom3.id:
                        atom3.x = keep.x
    return pairs

def assign_ff(model, itp):
    for i, atom in enumerate(model.atoms):
        at = itp.atoms[i]
        atom.atomtype = at.atomtype
        atom.cgnr = at.cgnr
        atom.q = at.q
        atom.m = at.m
        atom.atomtypeB = at.atomtypeB
        atom.qB = at.qB
        atom.mB = at.mB

# for forcefield2.py
def gen_dih_entry2(a1,a2,a3,a4,dihtype,entry=None,entryB=None,scDumA=1.0,scDumB=1.0):
    dihList = []
    ids = [a1.id,a2.id,a3.id,a4.id]
    if entry!=None:
        entry = entry[0].split()
        entry = [float(foo) for foo in entry]
    if entryB!=None:
        entryB = entryB[0].split()
        entryB = [float(foo) for foo in entryB]
    if (dihtype == 3 ):
        zeroesA = [0,0,0,0,0,0]
        zeroesB = [0,0,0,0,0,0]
    elif (dihtype == 2): # improper has no multiplicity
        angleA = 0
        angleB = 0
        if( entry != None ):
            angleA = entry[0]
        if( entryB != None ):
            angleB = entryB[0]
        zeroesA = [angleA,0]
        zeroesB = [angleB,0]
    else: # all the other types are the same
        multA = 0
        multB = 0
        angleA = 0
        angleB = 0
        if( entry != None ):
            angleA = entry[0]
            multA = entry[2]
        if( entryB != None ):
            angleB = entryB[0]
            multB = entryB[2]
        zeroesA = [angleA,0,multA]
        zeroesB = [angleB,0,multB]

    if( entry != None ):
        if( dihtype == 3 ):
            entry = [float(foo)*scDumA for foo in entry]
        else:
            entry[1] = scDumA*float(entry[1])
        dih = [a1,a2,a3,a4]+[dihtype]+[[dihtype]+entry]+[[dihtype]+zeroesA]
        dihList.append(dih)
    if( entryB != None ):
        if( dihtype == 3 ):
            entryB = [float(foo)*scDumB for foo in entryB]
        else:
            entryB[1] = scDumB*float(entryB[1])
        dih = [a1,a2,a3,a4]+[dihtype]+[[dihtype]+zeroesB]+[[dihtype]+entryB]
        dihList.append(dih)
#    if( entry==None ):
#        dih = ids+entryB+zeroesB
#        dihList.append(dih)
#    if( entryB==None ):
#       dih = ids+zeroesA+entry
#       dihList.append(dih)
    return dihList

# for forcefield.py
def gen_dih_entry(ids,entry=None,entryB=None,scDumA=1.0,scDumB=1.0):
    dihList = []
    if (ids[4] == 3 ):
        zeroesA = [0,0,0,0,0,0]
        zeroesB = [0,0,0,0,0,0]
    elif (ids[4] == 2): # improper has no multiplicity
        angleA = 0
        angleB = 0
        if( entry != None ):
            angleA = entry[0]
        if( entryB != None ):
            angleB = entryB[0]
        zeroesA = [angleA,0]
        zeroesB = [angleB,0]
    else: # all the other types are the same
        multA = 0
        multB = 0
        angleA = 0
        angleB = 0
        if( entry != None ):
            angleA = entry[0]
            multA = entry[2]
        if( entryB != None ):
            angleB = entryB[0]
            multB = entryB[2]
        zeroesA = [angleA,0,multA]
        zeroesB = [angleB,0,multB]

    if( entry != None ):
        if( ids[4] == 3 ):
            entry = [foo*scDumA for foo in entry]
        else:
            entry[1] = scDumA*entry[1]
        dih = ids+entry+zeroesA
        dihList.append(dih)
    if( entryB != None ):
        if( ids[4] == 3 ):
            entryB = [foo*scDumB for foo in entryB]
        else:
            entryB[1] = scDumB*entryB[1]
        dih = ids+zeroesB+entryB
        dihList.append(dih)
#    if( entry==None ):
#        dih = ids+entryB+zeroesB
#        dihList.append(dih)
#    if( entryB==None ):
#       dih = ids+zeroesA+entry
#       dihList.append(dih)
    return dihList

# for the newer forcefield2.py
def get_ff2_entry(ids, itp, gmx45=True, what = 'bond'):
    if what == 'bond':
        for b in itp.bonds:
            if (b[0].id == ids[0] and b[1].id == ids[1]) or \
               (b[1].id == ids[0] and b[0].id == ids[1]):
                out = xcopy.deepcopy(b)
                return out[3:][0]
    elif what == 'angle':
        for b in itp.angles:
            if (b[0].id == ids[0] and b[1].id == ids[1] and b[2].id == ids[2]) or \
               (b[2].id == ids[0] and b[1].id == ids[1] and b[0].id == ids[2]):
                out = xcopy.deepcopy(b)
                return out[4:][0]
    elif what == 'dihedral':
        if (ids[4] == 3):
            for b in itp.dihedrals:
                if (b[0].id == ids[0] and b[1].id == ids[1] and b[2].id == ids[2] and b[3].id == ids[3]) or \
                   (b[3].id == ids[0] and b[2].id == ids[1] and b[1].id == ids[2] and b[0].id == ids[3]):
                    if( ids[4] == b[4] ):
                        itp.dihedrals.remove(b)
                        out = xcopy.deepcopy(b)
                        return out[5:]
        elif (ids[4] == 9):
            for b in itp.dihedrals:
                if (b[0].id == ids[0] and b[1].id == ids[1] and b[2].id == ids[2] and b[3].id == ids[3]) or \
                   (b[3].id == ids[0] and b[2].id == ids[1] and b[1].id == ids[2] and b[0].id == ids[3]):
                    if( ids[4] == b[4] ):
                        itp.dihedrals.remove(b)
                        out = xcopy.deepcopy(b)
                        return out[5:]
        elif (ids[4]==1 and gmx45==True):
            for b in itp.dihedrals:
                if (b[0].id == ids[0] and b[1].id == ids[1] and b[2].id == ids[2] and b[3].id == ids[3]) or \
                   (b[3].id == ids[0] and b[2].id == ids[1] and b[1].id == ids[2] and b[0].id == ids[3]):
                    if( ids[4] == b[4] ):
                        itp.dihedrals.remove(b)
                        out = xcopy.deepcopy(b)
                        return out[5:]
        elif (ids[4]==2 or ids[4]==4 ): # improper
            for b in itp.dihedrals:
                sum = 0
                for b1 in range(0,4):
                    for b2 in range(0,4):
                        if (b[b1].id == ids[b2]):
                            sum += 1
                            break
                if ( (sum ==4) and (b[4]==2 or b[4]==4) ):
                    itp.dihedrals.remove(b)
                    out = xcopy.deepcopy(b)
                    return out[5:]
        elif (ids[4]==1 and gmx45==False ): # improper as proper in version < gmx45
            for b in itp.dihedrals:
                sum = 0
                for b1 in range(0,4):
                    for b2 in range(0,4):
                        if (b[b1].id == ids[b2]):
                            sum += 1
                            break
                if ( (sum == 4) and (b[4]==1) ):
                    itp.dihedrals.remove(b)
                    out = xcopy.deepcopy(b)
                    return out[5:]
    return None

# for the older forcefield.py
def get_ff_entry(ids, itp, gmx45=True, what = 'bond'):
    if what == 'bond':
        for b in itp.bonds:
            if (b[0] == ids[0] and b[1] == ids[1]) or \
               (b[1] == ids[0] and b[0] == ids[1]):
                return b[3:]
    elif what == 'angle':
        for b in itp.angles:
            if (b[0] == ids[0] and b[1] == ids[1] and b[2] == ids[2]) or \
               (b[2] == ids[0] and b[1] == ids[1] and b[0] == ids[2]):
                return b[4:]
    elif what == 'dihedral':
        if (ids[4] == 3):
            for b in itp.dihedrals:
                if (b[0] == ids[0] and b[1] == ids[1] and b[2] == ids[2] and b[3] == ids[3]) or \
                   (b[3] == ids[0] and b[2] == ids[1] and b[1] == ids[2] and b[0] == ids[3]):
                    if( ids[4] == b[4] ):
                        itp.dihedrals.remove(b)
                        return b[5:]
        elif (ids[4] == 9):
            for b in itp.dihedrals:
                if (b[0] == ids[0] and b[1] == ids[1] and b[2] == ids[2] and b[3] == ids[3]) or \
                   (b[3] == ids[0] and b[2] == ids[1] and b[1] == ids[2] and b[0] == ids[3]):
                    if( ids[4] == b[4] ):
                        itp.dihedrals.remove(b)
                        return b[5:]
        elif (ids[4]==1 and gmx45==True):
            for b in itp.dihedrals:
                if (b[0] == ids[0] and b[1] == ids[1] and b[2] == ids[2] and b[3] == ids[3]) or \
                   (b[3] == ids[0] and b[2] == ids[1] and b[1] == ids[2] and b[0] == ids[3]):
                    if( ids[4] == b[4] ):
                        itp.dihedrals.remove(b)
                        return b[5:]
        elif (ids[4]==2 or ids[4]==4 ): # improper
            for b in itp.dihedrals:
                sum = 0
                for b1 in range(0,4):
                    for b2 in range(0,4):
                        if (b[b1] == ids[b2]):
                            sum += 1
                            break
                if ( (sum ==4) and (b[4]==2 or b[4]==4) ):
                    itp.dihedrals.remove(b)
                    return b[5:]
        elif (ids[4]==1 and gmx45==False ): # improper as proper in version < gmx45
            for b in itp.dihedrals:
                sum = 0
                for b1 in range(0,4):
                    for b2 in range(0,4):
                        if (b[b1] == ids[b2]):
                            sum += 1
                            break
                if ( (sum == 4) and (b[4]==1) ):
                    itp.dihedrals.remove(b)
                    return b[5:]
    return None

def read_pairs_file(fn):
    l = open(fn).readlines()
    plst = []
    for line in l:
        plst.append(line.split())
    return plst

################################################################################33


def main(argv):

    version = "1.1"

    # define input/output files
    files= [
       FileOption("-l1", "r",["pdb"],"ligand1.pdb",""),
       FileOption("-l2", "r",["pdb"],"ligand2.pdb",""),
       FileOption("-itp1", "r",["itp"],"lig1.itp",""),
       FileOption("-itp2", "r",["itp"],"lig2.itp",""),
       FileOption("-n1", "r/o",["ndx"],"scaffold1" ,""),
       FileOption("-n2", "r/o",["ndx"],"scaffold2","" ),
       FileOption("-pairs", "r/o",["dat"],"pairs" ,""),
       FileOption("-oa", "w",["pdb"],"mergedA.pdb" ,""),
       FileOption("-ob", "w",["pdb"],"mergedB.pdb" ,""),
       FileOption("-oitp", "w",["itp"],"merged.itp","" ),
       FileOption("-ffitp", "w",["itp"],"ffmerged.itp" ,""),
#          FileOption("-ffitp1", "w/o",["itp"],"ffitp1.itp" ,""),
#          FileOption("-ffitp2", "w/o",["itp"],"ffitp2.itp" ,""),
       FileOption("-log", "w",["log"],"hybrid.log" ,""),
        ]

    # define options
    options=[
#          Option( "-bDist", "bool", "False", "use distance (no alignment) for the morphes"),
       Option( "-d", "float", "0.05", "distance (nm) between atoms to consider them morphable"),
       Option( "-scDUMm", "float", "1.0", "scale dummy masses using the counterpart atoms"),
       Option( "-scDUMa", "float", "1.0", "scale bonded dummy angle parameters"),
       Option( "-scDUMd", "float", "1.0", "scale bonded dummy dihedral parameters"),
       Option( "-deAng", "bool", "false", "decouple angles composed of 1 dummy and 2 non-dummies"),
#          Option( "-GMX45", "bool", "true", "set to noGMX45 for the topologies of earlier gromacs generations"),
       Option( "-fit", "bool", "false", "fit mol2 onto mol1, only works if pairs.dat is provided"),
       Option( "-split", "bool", "false", "split the topology into separate transitions"),
        ]

    help_text = ()

    # pass options, files and the command line to pymacs

    cmdl = Commandline( argv, options = options,
                   fileoptions = files,
                   program_desc = help_text,
                   check_for_existing_files = False )

    # deal with flags
    if cmdl.opt['-pairs'].is_set:
        read_pairs_from_file = True
    else:
        read_pairs_from_file = False

    gmx45 = True
    bFit = False
    deAng = False

    bScaleMass = False #inverse of mDUM
#        if mDUM==True:
#            bScaleMass = False
    if(cmdl['-fit']==True and read_pairs_from_file==True):
        bFit = True
    if(cmdl['-deAng']==True ):
        deAng = True

    if cmdl.opt['-n1'].is_set:
        read_from_idx = True
    else:
        read_from_idx = False
    bSplit = cmdl['-split']

    dist = 0.05
    if( cmdl.opt['-d'].is_set ):
        dist = cmdl['-d']

    scDuma = 1.0
    if( cmdl.opt['-scDUMa'].is_set ):
        scDuma = cmdl['-scDUMa']
    scDumd = 1.0
    if( cmdl.opt['-scDUMd'].is_set ):
        scDumd = cmdl['-scDUMd']
    scDUMm = 1.0
    if( cmdl.opt['-scDUMm'].is_set ):
        scDUMm = cmdl['-scDUMm']

    logfile = open(cmdl['-log'],'w')

    if read_from_idx and read_pairs_from_file:
        do_log(logfile, "Error: Can either read a pair list or scaffold index files!")
        do_log(logfile,"Exiting!")
        sys.exit(1)


    do_log(logfile,'Reading ligand 1 from: "%s"' % cmdl['-l1'])
    do_log(logfile,'Reading ligand 2 from: "%s"' % cmdl['-l2'])

    m1 = Model().read(cmdl['-l1'])
    m2 = Model().read(cmdl['-l2'])

    do_log(logfile,'Reading itp file 1 from: "%s"' % cmdl['-itp1'])
    do_log(logfile,'Reading itp file 2 from: "%s"' % cmdl['-itp2'])

    itp1 = ITPFile(cmdl['-itp1'])
    itp2 = ITPFile(cmdl['-itp2'])

    do_log(logfile,"Assigning forcefield parameters....")
    assign_ff(m1,itp1)
    assign_ff(m2,itp2)
    do_log(logfile,"Making pairs.....")

    if read_pairs_from_file:
        do_log(logfile,'Reading file with atom pairs: "%s"' % cmdl['-pairs'])
        plst = read_pairs_file(cmdl['-pairs'])
    else:
        plst = None

    if read_from_idx:
        do_log(logfile,'Reading scaffold index file: "%s"' % cmdl['-n1'])
#           grp1 = IndexFile(args['-n1']['fns']).dic['scaffold']
        grp1 = IndexFile(cmdl['-n1']).dic['scaffold']
        do_log(logfile,'Reading scaffold index file: "%s"' % cmdl['-n2'])
        grp2 = IndexFile(cmdl['-n2']).dic['scaffold']
        # now we add all atoms with bonds to scaffold atoms
        for b in itp1.bonds:
            if b[0] in grp1.ids and b[1] not in grp1.ids:
                grp1.ids.append(b[1])
                do_log(logfile,'Adding atom %s to scaffold 1' % m1.atoms[b[1]-1].name)
            elif b[1] in grp1.ids and b[0] not in grp1.ids:
                grp1.ids.append(b[0])
                do_log(logfile,'Adding atom %s to scaffold 1' % m1.atoms[b[0]-1].name)
        for b in itp2.bonds:
            if b[0] in grp2.ids and b[1] not in grp2.ids:
                grp2.ids.append(b[1])
                do_log(logfile,'Adding atom %s to scaffold 2' % m2.atoms[b[1]-1].name)
            elif b[1] in grp2.ids and b[0] not in grp2.ids:
                grp2.ids.append(b[0])
                do_log(logfile,'Adding atom %s to scaffold 2' % m2.atoms[b[0]-1].name)
        grps = [grp1.ids, grp2.ids]
    else:
        grps = None

    m3 = m1.copy() #m3 will contain all the atoms from m1, but with the coordinates of the matching atoms from m2
    m4 = m2.copy() #need to copy it when fitting

    # fitting
    if(bFit==True):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        n1 = []
        n2 = []
        for p in plst:
            n1.append(int(p[0])-1)
            n2.append(int(p[1])-1)

        do_log(logfile,'Superimposing mol2 on mol1')
        try:
            pdbName1,atomNameDict1 = reformatPDB(cmdl['-l1'],1)
            pdbName2,atomNameDict2 = reformatPDB(cmdl['-l2'],2)
            mol1 = Chem.MolFromPDBFile(pdbName1,removeHs=False,sanitize=False)
            mol2 = Chem.MolFromPDBFile(pdbName2,removeHs=False,sanitize=False)
            os.remove(pdbName1)
            os.remove(pdbName2)
            Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=list(zip(n2,n1)))
        except:
            pdbName1,atomNameDict1 = reformatPDB(cmdl['-l1'],1,bStrict=True)
            pdbName2,atomNameDict2 = reformatPDB(cmdl['-l2'],2,bStrict=True)
            mol1 = Chem.MolFromPDBFile(pdbName1,removeHs=False,sanitize=False)
            mol2 = Chem.MolFromPDBFile(pdbName2,removeHs=False,sanitize=False)
            os.remove(pdbName1)
            os.remove(pdbName2)
            Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=list(zip(n2,n1)))
        # adjust coordinates of m2
        adjustCoords(m2,mol2)
        restoreAtomNames(mol1,atomNameDict1)
        restoreAtomNames(mol2,atomNameDict2)

        do_log(logfile,'Superimposing mol1 on mol2')
        try:
            pdbName1,atomNameDict1 = reformatPDB(cmdl['-l1'],1)
            pdbName2,atomNameDict2 = reformatPDB(cmdl['-l2'],2)
            mol1 = Chem.MolFromPDBFile(pdbName1,removeHs=False,sanitize=False)
            mol2 = Chem.MolFromPDBFile(pdbName2,removeHs=False,sanitize=False)
            os.remove(pdbName1)
            os.remove(pdbName2)
            Chem.rdMolAlign.AlignMol(mol1,mol2,atomMap=list(zip(n1,n2)))
        except:
            pdbName1,atomNameDict1 = reformatPDB(cmdl['-l1'],1,bStrict=True)
            pdbName2,atomNameDict2 = reformatPDB(cmdl['-l2'],2,bStrict=True)
            mol1 = Chem.MolFromPDBFile(pdbName1,removeHs=False,sanitize=False)
            mol2 = Chem.MolFromPDBFile(pdbName2,removeHs=False,sanitize=False)
            os.remove(pdbName1)
            os.remove(pdbName2)
            Chem.rdMolAlign.AlignMol(mol1,mol2,atomMap=list(zip(n1,n2)))
        # adjust coordinates of m1
        adjustCoords(m3,mol1)
        restoreAtomNames(mol1,atomNameDict1)
        restoreAtomNames(mol2,atomNameDict2)
#       sys.exit(0)

    bDist = True
    pairs = make_pairs(m1, m2, m3, m4, bFit, bDist,dist, plst, grps)

    morphsA = [p[1] for p in pairs]
    morphsB = [p[0] for p in pairs]
    dumsA = []
    dumsA_nofit = []
    if(bFit==False):
        for atom in m2.atoms:
            if atom not in morphsA:
                dumsA.append(atom)
    else:
        for (atom,at) in zip(m2.atoms,m4.atoms):
            if atom not in morphsA:
                dumsA.append(atom)
                dumsA_nofit.append(at)
    dumsB = []
    for atom in m1.atoms:
        if atom not in morphsB:
            dumsB.append(atom)
    do_log(logfile, "Generated %d atom-atom pairs" % len(pairs))
    do_log(logfile,"Dummies in state A: %d" % len(dumsA))
    do_log(logfile,"Dummies in state B: %d" % len(dumsB))



    do_log(logfile,"Making B-states....")
    for a1, a2 in pairs:
        a1.atomtypeB = a2.atomtype
        a2.atomtypeB = a1.atomtype #this is my change to catch the DISAPPEARING dihedrals
        a1.nameB = a2.name
        a1.qB = a2.q
        a1.mB = a2.m
        a1.idB = a2.id
        a2.idB = a1.id #this is my change to catch the DISAPPEARING dihedrals
        a2.mB = a1.m
        do_log(logfile, "Atom....: %4d  %12s | %6.2f | %6.2f -> %12s | %6.2f | %6.2f" %\
               (a1.id, a1.atomtype, a1.q, a1.m, a1.atomtypeB, a1.qB, a1.mB))

    if( bFit==False):
        for atom in dumsA:
            atom.id_old = atom.id
            atom.nameB = atom.name
            atom.name = 'D'+atom.name
            atom.atomtypeB = atom.atomtype
            atom.atomtype = 'DUM_'+atom.atomtype
            atom.qB = atom.q
            atom.q = 0
            atom.mB = atom.m
            atom.m = atom.mB #1.

            atom.m = atom.m*scDUMm
            if( atom.m < 1.0 and atom.mB != 0.0): # exception for virtual particles
                atom.m = 1.0

            m1.residues[0].append(atom)
            m3.residues[0].append(atom)
            do_log(logfile, "Dummy...: %4d  %12s | %6.2f | %6.2f -> %12s | %6.2f | %6.2f" %\
                   (atom.id, atom.atomtype, atom.q, atom.m, atom.atomtypeB, atom.qB, atom.mB))
    else:
        for (atom,at) in zip(dumsA,dumsA_nofit):
            atom.id_old = atom.id
            atom.nameB = atom.name
            atom.name = 'D'+atom.name
            atom.atomtypeB = atom.atomtype
            atom.atomtype = 'DUM_'+atom.atomtype
            atom.qB = atom.q
            atom.q = 0
            atom.mB = atom.m
            atom.m = atom.mB #1.

            atom.m = atom.m*scDUMm
            if( atom.m < 1.0 and atom.mB != 0.0): # exception for virtual particles
                atom.m = 1.0

            m1.residues[0].append(atom)
            at.id_old = at.id
            at.nameB = at.name
            at.name = 'D'+at.name
            at.atomtypeB = at.atomtype
            at.atomtype = 'DUM_'+at.atomtype
            at.qB = at.q
            at.q = 0
            at.mB = at.m
            at.m = at.mB #1.

            at.m = at.m*scDUMm
            if( at.m < 1.0 and at.mB != 0.0): # exception for virtual particles
                at.m = 1.0

            m3.residues[0].append(at)
            do_log(logfile, "Dummy...: %4d  %12s | %6.2f | %6.2f -> %12s | %6.2f | %6.2f" %\
                   (atom.id, atom.atomtype, atom.q, atom.m, atom.atomtypeB, atom.qB, atom.mB))


    for atom in dumsB:
        atom.atomtypeB = 'DUM_'+atom.atomtype
        atom.qB = 0
        atom.mB = atom.m #1.

        atom.mB = atom.mB*scDUMm
        if( atom.mB < 1.0 and atom.m != 0.0): # exception for virtual particles
            atom.mB = 1.0

        do_log(logfile, "Dummy...: %4d  %12s | %6.2f | %6.2f -> %12s | %6.2f | %6.2f" %\
               (atom.id, atom.atomtype, atom.q, atom.m, atom.atomtypeB, atom.qB, atom.mB))

    id_dicAB = {}
    id_dicBA = {}
    for atom in m1.atoms:
        if hasattr(atom,"idB"):
            id_dicAB[atom.id] = atom.idB
            id_dicBA[atom.idB] = atom.id
        if hasattr(atom,"id_old"):
            id_dicAB[atom.id] = atom.id_old
            id_dicBA[atom.id_old] = atom.id

    do_log(logfile, "Generating bonded parameters....")


    # go over bonds
    newbonds = []

    for b in itp1.bonds:
        id1 = b[0].id
        id2 = b[1].id
        a1 = m1.atoms[id1-1]
        a2 = m1.atoms[id2-1]
        bOk = False
        if hasattr(a1,"idB") and hasattr(a2,"idB"):
            idB1 = a1.idB
            idB2 = a2.idB
            entr = get_ff2_entry([idB1, idB2], itp2, gmx45, what= 'bond')
            if entr is not None:
                newbonds.append ([b[0],b[1],b[2],[b[2]]+b[3],[b[2]]+entr])
                bOk = True
            else:
                bOk = False
        elif a1.atomtypeB[:3] == 'DUM' or a2.atomtypeB[:3] == 'DUM':
            entr = get_ff2_entry([a1.id, a2.id], itp1, gmx45, what= 'bond')
            if entr is not None:
                newbonds.append ([b[0],b[1],b[2],[b[2]]+b[3],[b[2]]+entr])
                bOk = True
            else:
                bOk = False
        else:
            newbonds.append(b)
            bOk = True

        if not bOk:
            do_log(logfile, "Error: Something went wrong while assigning bonds!")
            do_log(logfile, "A-> Atom1: %d-%s Atom2: %d-%s" %(a1.id, a1.name, a2.id, a2.name))
            do_log(logfile, "B-> Atom1: %d-%s Atom2: %d-%s" %(a1.idB, a1.nameB, a2.idB, a2.nameB))
            do_log(logfile,"Exiting....")
            sys.exit(1)

    # angles
    newangles = []
    decoupAngles = []
    for b in itp1.angles:
        angtype = b[3]
        ####### explanation for angle parameters ############
        # for angtype 1, entry is [angle,force_const]
        # for angtype 5, entry is [5,angle,force_const,bond,force_const]
        # write_angles() in forcefield2.py expects entry [angtype,angle,force_const,...]
        id1 = b[0].id
        id2 = b[1].id
        id3 = b[2].id
        a1 = m1.atoms[id1-1]
        a2 = m1.atoms[id2-1]
        a3 = m1.atoms[id3-1]
        bOk = False
        if hasattr(a1,"idB") and hasattr(a2,"idB") and hasattr(a3,"idB"):
            idB1 = a1.idB
            idB2 = a2.idB
            idB3 = a3.idB
            entr = get_ff2_entry([idB1, idB2, idB3], itp2, gmx45, what= 'angle')
            if entr is not None:
                if angtype==1:
                    newangles.append ([b[0],b[1],b[2],angtype,[angtype]+b[4],[angtype]+entr])
                else:
                    newangles.append ([b[0],b[1],b[2],angtype,b[4],entr])
                bOk = True
            else:
                bOk = False
        elif a1.atomtypeB[:3] == 'DUM' or \
                 a2.atomtypeB[:3] == 'DUM' or \
                 a3.atomtypeB[:3] == 'DUM':
            entr = get_ff2_entry([a1.id, a2.id, a3.id], itp1, gmx45, what= 'angle')
            if (angtype!=1) and (entr is not None): # remove the angtype from the entr
                entr = entr[1:]
            if entr is not None:
                if( (a1.atomtypeB[:3] != 'DUM' and a2.atomtypeB[:3] != 'DUM') \
                    or (a1.atomtypeB[:3] != 'DUM' and a3.atomtypeB[:3] != 'DUM') \
                    or (a2.atomtypeB[:3] != 'DUM' and a3.atomtypeB[:3] != 'DUM') ):
#                   if( (a1.atomtypeB[:3] != 'DUM') or (a2.atomtypeB[:3] != 'DUM') or (a3.atomtypeB[:3] != 'DUM') ):
                    entr[1] = scDuma*entr[1]
                    if angtype==5:
                        entr[3] = scDuma*entr[3]
                    if(deAng==True):
                        if( a1.atomtypeB[:3] == 'DUM' ):
                            if( findIDinList(id1,decoupAngles)==True ):
                                entr[1] = 0.0
                                if angtype==5:
                                    entr[3] = 0.0
                            else:
                                decoupAngles.append(id1)
                        elif( a2.atomtypeB[:3] == 'DUM' ):
                            if( findIDinList(id2,decoupAngles)==True ):
                                entr[1] = 0.0
                                if angtype==5:
                                    entr[3] = 0.0
                            else:
                                decoupAngles.append(id2)
                        elif( a3.atomtypeB[:3] == 'DUM' ):
                            if( findIDinList(id3,decoupAngles)==True ):
                                entr[1] = 0.0
                                if angtype==5:
                                    entr[3] = 0.0
                            else:
                                decoupAngles.append(id3)
                if angtype==1:
                    newangles.append ([b[0],b[1],b[2],angtype,[angtype]+b[4],[angtype]+entr])
                else:
                    newangles.append ([b[0],b[1],b[2],angtype,b[4],[angtype]+entr])
                bOk = True
            else:
                bOk = False
        else:
            newangles.append(b)
            bOk = True

        if not bOk:
            do_log(logfile, "Error: Something went wrong while assigning angles!")
            do_log(logfile, "A-> Atom1: %d-%s Atom2: %d-%s Atom3: %d-%s" \
                   %(a1.id, a1.name, a2.id, a2.name, a3.id, a3.name))
            do_log(logfile, "B-> Atom1: %d-%s Atom2: %d-%s Atom3: %d-%s" \
                   %(a1.idB, a1.nameB, a2.idB, a2.nameB, a3.idB, a3.nameB))
            do_log(logfile,"Exiting....")
            sys.exit(1)


    #############################
    ############ VG #############
    # COMPLETE DIHEDRAL REWRITE #
    #############################
    newdihedrals = []
    cpItp1 = xcopy.deepcopy(itp1)
    cpItp2 = xcopy.deepcopy(itp2)
    for b in itp1.dihedrals:
        id1 = b[0].id
        id2 = b[1].id
        id3 = b[2].id
        id4 = b[3].id
        dih_type = b[4]
        a1 = m1.atoms[id1-1]
        a2 = m1.atoms[id2-1]
        a3 = m1.atoms[id3-1]
        a4 = m1.atoms[id4-1]
        entrA = get_ff2_entry([id1, id2, id3, id4, dih_type], cpItp1, gmx45, what= 'dihedral')
        bOk = False
        if hasattr(a1,"idB") and hasattr(a2,"idB") and \
               hasattr(a3,"idB") and hasattr(a4,"idB"):
            # switch the A state off
            dih = gen_dih_entry2(a1, a2, a3, a4, dih_type, entrA,None)
            newdihedrals.extend(dih)
            bOk = True
        else:
            # switch the B state on
            if a1.atomtypeB[:3] == 'DUM' or \
                     a2.atomtypeB[:3] == 'DUM' or \
                     a3.atomtypeB[:3] == 'DUM' or \
                     a4.atomtypeB[:3] == 'DUM':
                if entrA is not None:
                    dih = gen_dih_entry2(a1, a2, a3, a4, dih_type,entrA,None)
                    newdihedrals.extend(dih)
#                        if( (a1.atomtypeB[:3] != 'DUM' and a2.atomtypeB[:3] != 'DUM' and a3.atomtypeB[:3] != 'DUM') \
#                            or (a1.atomtypeB[:3] != 'DUM' and a2.atomtypeB[:3] != 'DUM' and a4.atomtypeB[:3] != 'DUM') \
#                            or (a2.atomtypeB[:3] != 'DUM' and a3.atomtypeB[:3] != 'DUM' and a4.atomtypeB[:3] != 'DUM') \
#                            or (a1.atomtypeB[:3] != 'DUM' and a3.atomtypeB[:3] != 'DUM' and a4.atomtypeB[:3] != 'DUM') ):
                    if( a1.atomtypeB[:3] != 'DUM' or a2.atomtypeB[:3] != 'DUM' or a3.atomtypeB[:3] != 'DUM' or a4.atomtypeB[:3] != 'DUM' ):
#                        if( (a1.atomtypeB[:3] != 'DUM' and a2.atomtypeB[:3] != 'DUM') \
#                           or (a1.atomtypeB[:3] != 'DUM' and a3.atomtypeB[:3] != 'DUM') \
#                           or (a1.atomtypeB[:3] != 'DUM' and a4.atomtypeB[:3] != 'DUM') \
#                           or (a2.atomtypeB[:3] != 'DUM' and a3.atomtypeB[:3] != 'DUM') \
#                           or (a2.atomtypeB[:3] != 'DUM' and a4.atomtypeB[:3] != 'DUM') \
#                           or (a3.atomtypeB[:3] != 'DUM' and a4.atomtypeB[:3] != 'DUM') ):
                        if dih_type==2 or dih_type==4:# or dih_type==1: # disable improper for dummy-nondummy
                            dih = gen_dih_entry2(a1, a2, a3, a4, dih_type,None,entrA,1.0,0.0)
                        else:
                            dih = gen_dih_entry2(a1, a2, a3, a4, dih_type,None,entrA,1.0,scDumd)
                    else:
                        dih = gen_dih_entry2(a1, a2, a3, a4, dih_type,None,entrA)
                    newdihedrals.extend(dih)
                    bOk = True
                else:
                    bOk = False
            else:
                newdihedrals.append(b)
                bOk = True

        if not bOk:
            do_log(logfile, "Error: Something went wrong while assigning dihedrals!")
            do_log(logfile, "A-> Atom1: %d-%s Atom2: %d-%s Atom3: %d-%s Atom3: %d-%s" \
                   %(a1.id, a1.name, a2.id, a2.name, a3.id, a3.name, a4.id, a4.name))
            do_log(logfile, "B-> Atom1: %d-%s Atom2: %d-%s Atom3: %d-%s Atom3: %d-%s" \
                   %(a1.idB, a1.nameB, a2.idB, a2.nameB, a3.idB, a3.nameB, a4.idB, a4.nameB))
            do_log(logfile,"Exiting....")
            sys.exit(1)


    # second molecule dihedrals
    for b in itp2.dihedrals:
        id1 = b[0].id
        id2 = b[1].id
        id3 = b[2].id
        id4 = b[3].id
        aB1 = m2.atoms[id1-1]
        aB2 = m2.atoms[id2-1]
        aB3 = m2.atoms[id3-1]
        aB4 = m2.atoms[id4-1]
        newid1 = id_dicBA[b[0].id]
        newid2 = id_dicBA[b[1].id]
        newid3 = id_dicBA[b[2].id]
        newid4 = id_dicBA[b[3].id]
        a1 = m1.atoms[newid1-1]
        a2 = m1.atoms[newid2-1]
        a3 = m1.atoms[newid3-1]
        a4 = m1.atoms[newid4-1]
        dih_type = b[4]
        entrB = get_ff2_entry([b[0].id,b[1].id,b[2].id,b[3].id, dih_type], cpItp2, gmx45, what='dihedral')
        bOk = False
        if hasattr(aB1,"idB") and hasattr(aB2,"idB") and \
               hasattr(aB3,"idB") and hasattr(aB4,"idB"):
            # switch the B state off
            dih = gen_dih_entry2(a1,a2,a3,a4, dih_type,None,entrB)
            newdihedrals.extend(dih)
            bOk = True
        else:
            # switch the A state on
            if a1.atomtype.startswith('DUM') or \
               a2.atomtype.startswith('DUM') or \
               a3.atomtype.startswith('DUM') or \
               a4.atomtype.startswith('DUM'):
                if entrB is not None:
                    dih = gen_dih_entry2(a1,a2,a3,a4,dih_type,None,entrB)
                    newdihedrals.extend(dih)
#                        if( (a1.atomtype.startswith('DUM')==False and a2.atomtype.startswith('DUM')==False and a3.atomtype.startswith('DUM')==False) \
#                           or (a1.atomtype.startswith('DUM')==False and a3.atomtype.startswith('DUM')==False and a4.atomtype.startswith('DUM')==False) \
#                           or (a2.atomtype.startswith('DUM')==False and a3.atomtype.startswith('DUM')==False and a4.atomtype.startswith('DUM')==False) \
#                           or (a1.atomtype.startswith('DUM')==False and a2.atomtype.startswith('DUM')==False and a4.atomtype.startswith('DUM')==False) ):
                    if( a1.atomtype.startswith('DUM')==False or a2.atomtype.startswith('DUM')==False \
                        or a3.atomtype.startswith('DUM')==False or a4.atomtype.startswith('DUM')==False ):
#                        if( (a1.atomtype.startswith('DUM')==False and a2.atomtype.startswith('DUM')==False) \
#                           or (a1.atomtype.startswith('DUM')==False and a3.atomtype.startswith('DUM')==False) \
#                           or (a1.atomtype.startswith('DUM')==False and a4.atomtype.startswith('DUM')==False) \
#                           or (a2.atomtype.startswith('DUM')==False and a3.atomtype.startswith('DUM')==False) \
#                           or (a2.atomtype.startswith('DUM')==False and a4.atomtype.startswith('DUM')==False) \
#                           or (a3.atomtype.startswith('DUM')==False and a4.atomtype.startswith('DUM')==False) ):
                        if dih_type==2 or dih_type==4:# or dih_type==1: # disable improper for dummy-nondummy
                            dih = gen_dih_entry2(a1,a2,a3,a4,dih_type,entrB,None,0.0,1.0)
                        else:
                            dih = gen_dih_entry2(a1,a2,a3,a4,dih_type,entrB,None,scDumd,1.0)
                    else:
                        dih = gen_dih_entry2(a1,a2,a3,a4,dih_type,entrB,None)
                    newdihedrals.extend(dih)
                    bOk = True
                else:
                    bOk = False
            else:
                newdihedrals.append(b)
                bOk = True

        if not bOk:
            do_log(logfile, "Error: Something went wrong while assigning dihedrals!")
            do_log(logfile, "A-> Atom1: %d-%s Atom2: %d-%s Atom3: %d-%s Atom3: %d-%s" \
                   %(a1.id, a1.name, a2.id, a2.name, a3.id, a3.name, a4.id, a4.name))
            do_log(logfile, "B-> Atom1: %d-%s Atom2: %d-%s Atom3: %d-%s Atom3: %d-%s" \
                   %(a1.idB, a1.nameB, a2.idB, a2.nameB, a3.idB, a3.nameB, a4.idB, a4.nameB))
            do_log(logfile,"Exiting....")
            sys.exit(1)
    #############################
    # COMPLETE DIHEDRAL REWRITE #
    #############################


    # vsites2: stateA
    newvsites2 = []
    has_vsites2 = False
#        if(itp1.vsites2):
#            for b in itp1.vsites2:
    if(itp1.virtual_sites2):
        has_vsites2 = True
        for b in itp1.virtual_sites2:
            id1 = b[0].id
            id2 = b[1].id
            id3 = b[2].id
            a1 = m1.atoms[id1-1]
            a2 = m1.atoms[id2-1]
            a3 = m1.atoms[id3-1]
            newvsites2.append(b)

    # vsites3: stateA
    newvsites3 = []
    has_vsites3 = False
    if(itp1.virtual_sites3):
        has_vsites3 = True
        for b in itp1.virtual_sites3:
            id1 = b[0].id
            id2 = b[1].id
            id3 = b[2].id
            id4 = b[3].id
            a1 = m1.atoms[id1-1]
            a2 = m1.atoms[id2-1]
            a3 = m1.atoms[id3-1]
            a4 = m1.atoms[id4-1]
            newvsites3.append(b)

    # exclusions: stateA
    newexclusions = []
    has_exclusions = False
    if(itp1.exclusions):
        has_exclusions = True
        for excl in itp1.exclusions:
            exclToAdd = []
            for ex in excl:
                newid = ex.id
                newa = m1.atoms[newid-1]
                exclToAdd.append(newa)
            newexclusions.append( exclToAdd )

    # now we have all parameter for pairs
    # let's go for the dummies
    for b in itp2.bonds:
        newid1 = id_dicBA[b[0].id]
        newid2 = id_dicBA[b[1].id]
        a1 = m1.atoms[newid1-1]
        a2 = m1.atoms[newid2-1]
        if a1.atomtype.startswith('DUM') or \
           a2.atomtype.startswith('DUM'):
            newbonds.append( [a1, a2, 1, [1]+b[-1], [1]+b[-1]] )

    decoupAngles = []
    for b in itp2.angles:
        angtype = b[3]
        entry = b[-1]
        # for type 1
        paramA = [ 1,entry[0],0.0 ]
        paramAscDum = [ 1,entry[0],scDuma*float(entry[1]) ]
        paramB = [ 1,entry[0],entry[1] ]
        # for type 5
        if angtype==5:
            paramA = [ 5,entry[1],0.0,entry[3],0.0 ]
            paramAscDum = [ 5,entry[1],scDuma*float(entry[2]),entry[3],scDuma*float(entry[4]) ]
            paramB = [ 5,entry[1],entry[2],entry[3],entry[4] ]

        newid1 = id_dicBA[b[0].id]
        newid2 = id_dicBA[b[1].id]
        newid3 = id_dicBA[b[2].id]
        a1 = m1.atoms[newid1-1]
        a2 = m1.atoms[newid2-1]
        a3 = m1.atoms[newid3-1]
        if a1.atomtype.startswith('DUM') or \
           a2.atomtype.startswith('DUM') or \
           a3.atomtype.startswith('DUM'):
            if( (a1.atomtype.startswith('DUM')==False and a2.atomtype.startswith('DUM')==False) \
               or (a1.atomtype.startswith('DUM')==False and a3.atomtype.startswith('DUM')==False) \
               or (a2.atomtype.startswith('DUM')==False and a3.atomtype.startswith('DUM')==False) ):
#               if( a1.atomtype.startswith('DUM')==False or a2.atomtype.startswith('DUM')==False or a3.atomtype.startswith('DUM')==False ):

                if(deAng==True):
                    if( a1.atomtype.startswith('DUM')==True ):
                        if( findIDinList(newid1,decoupAngles)==True ):
                            newangles.append([ a1,a2,a3,angtype, paramA, paramB ])
                        else:
                            newangles.append([ a1,a2,a3,angtype, paramAscDum, paramB ])
                            decoupAngles.append(newid1)
                    elif( a2.atomtype.startswith('DUM')==True ):
                        if( findIDinList(newid2,decoupAngles)==True ):
                            newangles.append([ a1,a2,a3,angtype, paramA, paramB ])
                        else:
                            newangles.append([ a1,a2,a3,angtype, paramAscDum, paramB ])
                            decoupAngles.append(newid2)
                    elif( a3.atomtype.startswith('DUM')==True ):
                        if( findIDinList(newid3,decoupAngles)==True ):
                            newangles.append([ a1,a2,a3,angtype, paramA, paramB ])
                        else:
                            newangles.append([ a1,a2,a3,angtype, paramAscDum, paramB ])
                            decoupAngles.append(newid3)
                else:
                    newangles.append([ a1,a2,a3,angtype, paramAscDum, paramB ])
            else:
                newangles.append([ a1,a2,a3,angtype, paramB, paramB ])

    # dihedrals already accounted for both states
#       cpcpItp2 = xcopy.deepcopy(cpItp2) # dirty hack
#       for b in cpcpItp2.dihedrals:
#           newid1 = id_dicBA[b[0]]
#           newid2 = id_dicBA[b[1]]
#           newid3 = id_dicBA[b[2]]
#           newid4 = id_dicBA[b[3]]
#           a1 = m1.atoms[newid1-1]
#           a2 = m1.atoms[newid2-1]
#           a3 = m1.atoms[newid3-1]
#           a4 = m1.atoms[newid4-1]
#           dih_type = b[4]
#           entrB = get_ff_entry([b[0],b[1],b[2],b[3], dih_type], cpItp2, gmx45, what='dihedral')
#           if a1.atomtype.startswith('DUM') or \
#              a2.atomtype.startswith('DUM') or \
#              a3.atomtype.startswith('DUM') or \
#              a4.atomtype.startswith('DUM'):
#               newdihedrals.append( [newid1, newid2, newid3, newid4, b[4]] + b[5:] + b[5:] )
#                dih = gen_dih_entry([newid1, newid2, newid3, newid4, dih_type],None,entrB)
#               newdihedrals.extend(dih)

    # vsites2: stateB
    if(itp2.virtual_sites2):
        has_vsites2 = True
        for b in itp2.virtual_sites2:
            newid1 = id_dicBA[b[0].id]
            newid2 = id_dicBA[b[1].id]
            newid3 = id_dicBA[b[2].id]
            #VG: this is not tested, use with caution
#               print "UNTESTED PART FOR THE NEW PMX %s %s %s" %(newid1,newid2,newid3)
#                newid1 = id_dicBA[b[0]]
#                newid2 = id_dicBA[b[1]]
#                newid3 = id_dicBA[b[2]]
            a1 = m1.atoms[newid1-1]
            a2 = m1.atoms[newid2-1]
            a3 = m1.atoms[newid3-1]
            vsiteToAdd = [a1,a2,a3,b[3],b[4]]
            if( check_if_vsite2_exists( newvsites2, vsiteToAdd )==False ):
                newvsites2.append( [a1, a2, a3, b[3], b[4]] )
#                newvsites2.append( [newid1, newid2, newid3, b[3], b[4]] )

    # vsites3: stateB
    if(itp2.virtual_sites3):
        has_vsites3 = True
        for b in itp2.virtual_sites3:
            newid1 = id_dicBA[b[0].id]
            newid2 = id_dicBA[b[1].id]
            newid3 = id_dicBA[b[2].id]
            newid4 = id_dicBA[b[3].id]
            a1 = m1.atoms[newid1-1]
            a2 = m1.atoms[newid2-1]
            a3 = m1.atoms[newid3-1]
            a4 = m1.atoms[newid4-1]
            vsiteToAdd = [a1,a2,a3,a4,b[4],b[5]]
            if( check_if_vsite3_exists( newvsites3, vsiteToAdd )==False ):
                newvsites3.append( [a1, a2, a3, a4, b[4], b[5]] )

    # exclusions: stateB
    if(itp2.exclusions):
        has_exclusions = True
        for excl in itp2.exclusions:
            exclToAdd = []
            for ex in excl:
                newid = id_dicBA[ex.id]
                newa = m1.atoms[newid-1]
                exclToAdd.append(newa)
            if( check_if_exclusion_valid( newexclusions, exclToAdd )==True ):
                newexclusions.append( exclToAdd )

    # make pairs
    newpairs = []
    pp = []
    for p in itp1.pairs:
        newpairs.append( p )
        pp.append( (p[0].id,p[1].id) )
    for p in itp2.pairs:
        newid1 = id_dicBA[p[0].id]
        newid2 = id_dicBA[p[1].id]
        a1 = m1.atoms[newid1-1]
        a2 = m1.atoms[newid2-1]
        if (newid1, newid2) not in pp and \
           (newid2, newid1) not in pp:
            newpairs.append([ a1, a2, 1] )

    do_log(logfile, "Generating new itp file")

    newitp = ITPFile(filename=None)
    newitp.name = itp1.name
    newitp.nrexcl = itp1.nrexcl
    newitp.atoms = m1.atoms
    newitp.residues = m1.residues
    for i, atom in enumerate(newitp.atoms):
        atom.cgnr = i +1
    newitp.bonds = newbonds
    newitp.pairs = newpairs
    newitp.angles = newangles
    newitp.dihedrals = newdihedrals
    newitp.virtual_sites2 = newvsites2
    newitp.has_vsites2 = has_vsites2
    newitp.virtual_sites3 = newvsites3
    newitp.has_vsites3 = has_vsites3
    newitp.exclusions = newexclusions
    newitp.has_exclusions = has_exclusions
    # get charges
    qA, qB = sum_charge_of_states( newitp )
    qA_mem = xcopy.deepcopy( qA )
    qB_mem = xcopy.deepcopy( qB )

    do_log(logfile, 'Writing new itp file: "%s"' % cmdl['-oitp'])
    newitp.write(cmdl['-oitp'], target_qB = qB)
    #do_log(logfile, 'Writing new structure file: "%s"' % args['-o']['fns'])
    #m1.write(args['-o'])
    do_log(logfile, 'Writing dummy forcefield file: "%s"' % cmdl['-ffitp'])

    if bSplit==True: # write splitted topology
        root, ext = os.path.splitext(cmdl['-oitp'])
        out_file_qoff = root+'_qoff'+ext
        out_file_vdw = root+'_vdw'+ext
        out_file_qon = root+'_qon'+ext


        print('------------------------------------------------------')
        print('log_> Creating splitted topologies............')
        print('log_> Making "qoff" topology : "%s"' % out_file_qoff)
        contQ = xcopy.deepcopy(qA_mem)
        newitp.write( out_file_qoff, stateQ = 'AB', stateTypes = 'AA', dummy_qB='off',
                      scale_mass = bScaleMass, target_qB = qA, stateBonded = 'AA', full_morphe = False )
        print('log_> Charge of state A: %g' % newitp.qA)
        print('log_> Charge of state B: %g' % newitp.qB)

        print('------------------------------------------------------')
        print('log_> Making "vdw" topology : "%s"' % out_file_vdw)
        contQ = xcopy.deepcopy(qA_mem)
        newitp.write( out_file_vdw, stateQ = 'BB', stateTypes = 'AB', dummy_qA='off', dummy_qB = 'off',
                      scale_mass = bScaleMass, target_qB = contQ, stateBonded = 'AB' , full_morphe = False)
        print('log_> Charge of state A: %g' % newitp.qA)
        print('log_> Charge of state B: %g' % newitp.qB)
        print('------------------------------------------------------')

        print('log_> Making "qon" topology : "%s"' % out_file_qon)
        newitp.write( out_file_qon, stateQ = 'BB', stateTypes = 'BB', dummy_qA='off', dummy_qB = 'on',
                      scale_mass = bScaleMass, target_qB = qB_mem,  stateBonded = 'BB' , full_morphe = False)
        print('log_> Charge of state A: %g' % newitp.qA)
        print('log_> Charge of state B: %g' % newitp.qB)
        print('------------------------------------------------------')


    # write ffitp
    fp = open(cmdl['-ffitp'],'w')
    dd = []
    print('[ atomtypes ]', file=fp)
    for atom in m1.atoms:
        if atom.atomtype.startswith('DUM') and atom.atomtype not in dd:
            print('%8s %12.6f %12.6f %3s %12.6f %12.6f' % \
                  (atom.atomtype, 0, 0, 'A',0,0), file=fp)
            dd.append(atom.atomtype)
        elif atom.atomtypeB.startswith('DUM') and atom.atomtypeB not in dd:
            print('%8s %12.6f %12.6f %3s %12.6f %12.6f' % \
                  (atom.atomtypeB, 0, 0, 'A',0,0), file=fp)
            dd.append(atom.atomtypeB)

    # write merged pdb
    m1.write(cmdl['-oa'])
    m3.write(cmdl['-ob'])

main( sys.argv )
