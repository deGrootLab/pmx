import sys, os
import copy as cp
from pmx import *
from pmx.ndx import *
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit.Chem import MCS
# for new RDKit-2018 use below
#from rdkit.Chem import rdFMCS

from rdkit.Chem import FragmentMatcher, Crippen, rdmolops
from random import randint

def remove_sigmaHoles( n1,n2,sigmaHoleID1,sigmaHoleID2):
    n1new = []
    n2new = []
    for (i,j) in zip(n1,n2):
        if (i+1 not in sigmaHoleID1) and (j+1 not in sigmaHoleID2):
            n1new.append(i)
            n2new.append(j)
    return(n1new,n2new)

def alignOnSubset(mol1,mol2,constrMap):
####### use subsets of atoms for o3a alignment ########
    # only heavy atoms can be constraints
    rem = []
    for c in constrMap:
        a1 = mol1.GetAtomWithIdx(c[0])
        a2 = mol2.GetAtomWithIdx(c[1])
        id1 = a1.GetAtomicNum()
        id2 = a2.GetAtomicNum()
        if((id1==1) or (id2==1)):
            rem.append(c)
    # remove
    for i in rem:
        constrMap.remove(i)

    # align
    mol1_crippen = crippen(mol1)
    mol2_crippen = crippen(mol2)
    pyO3A = AllChem.GetCrippenO3A(mol1,mol2,mol1_crippen,mol2_crippen,constraintMap=constrMap,options=0) # mol1=probe, mol2=ref
    pyO3A.Align()
    return(pyO3A.Score())

def getAttr(n1,n2,iStart,iEnd):
    jStart = None
    jEnd = None
    for foo,bar in zip(n1,n2):
        if(foo==iStart):
            jStart = bar
        if(foo==iEnd):
            jEnd = bar
    return(jStart,jEnd)

def getMapped(n1,n2,ind1):
    ind2 = None
    for id1,id2 in zip(n1,n2):
        if id1==ind1:
            ind2=id2
            return(ind2)
    return(ind2)

def getMappedList(n1,n2,indList1):
    indList2 = []
    for ind1 in indList1:
        for id1,id2 in zip(n1,n2):
            if id1==ind1:
                indList2.append(id2)
                break
    return(indList2)

def checkNeighDist(mol1,mol2,nb1,n1,n2,rem1,rem2):
    c1 = mol1.GetConformer()
    c2 = mol2.GetConformer()
    for neigh1 in nb1:
        ind1 = neigh1.GetIdx()
        if ind1 in n1:
            ind2 = getMapped(n1,n2,ind1)
            if ind2!=None:
                dist = distance_based(mol1,mol2,0.0,id1=[ind1],id2=[ind2],calcOnly=True)
#               print ind1,ind2,dist
                if dist > 0.15: # 0.15 nm distance should be forgiving enough, but also catch the cases of chirality inversions
                    if (ind1 not in rem1) and (ind2 not in rem2):
                        rem1.append(ind1)
                        rem2.append(ind2)

def localEnvironment(mol,i,n):
    a1 = mol.GetAtomWithIdx(i)
    nb1 = a1.GetNeighbors()
    # 1-2 and 1-3
    list12 = []
    list13 = []
    for a2 in nb1:
        ind1 = a2.GetIdx()
        if ind1 in n:
            list12.append(ind1)
        nb2 = a2.GetNeighbors()
        for a3 in nb2:
            ind2 = a3.GetIdx()
            if (ind2 in n) and (ind2 not in list12) and (ind2 not in list13) and (ind2 != i):
                list13.append(ind2)
    # final list
    listFinal = [i] + list12 + list13
    return(listFinal)


def checkChiral(mol1,mol2,n1,n2):
    # create constraint for alignment
    constrMap = list(zip(n1,n2))
    # align on the subset n1,n2
    Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=list(zip(n2,n1)))
    bonds1 = mol1.GetBonds()
    bonds2 = mol2.GetBonds()
    # create two dictionaries for mappings [n1] = n2
    dictn1n2 = mappingDict( n1,n2 )
    dictn2n1 = mappingDict( n2,n1 )

    rem1 = []
    rem2 = []
    for i1,i2 in zip(n1,n2):
        a1 = mol1.GetAtomWithIdx(i1)
        a2 = mol2.GetAtomWithIdx(i2)
        chirality1 = a1.GetChiralTag()
        chirality2 = a2.GetChiralTag()
        # also checking if some bonds are not single
        bonds1 = a1.GetBonds()
        bonds2 = a2.GetBonds()
        bNonsingleBond1 = check_nonsingle_bonds( mol1, mol2, bonds1, n1, n2, dictn1n2 )
        bNonsingleBond2 = check_nonsingle_bonds( mol2, mol1, bonds2, n2, n1, dictn2n1 )
        if (str(chirality1) != "CHI_UNSPECIFIED") or (bNonsingleBond1==True):
#           print "1",i1,chirality1
            # try fitting locally on the 1-2, 1-3 atoms
            localEnv1 = localEnvironment(mol1,i1,n1)
            if len(localEnv1)>2:
                localEnv2 = getMappedList(n1,n2,localEnv1)
                Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=list(zip(localEnv2,localEnv1)))
#               alignOnSubset(mol1,mol2,zip(localEnv1,localEnv2))
            # get the neighbours
            nb1 = a1.GetNeighbors()
            # check if the matched neighbours in the aligned molecules are too far
            checkNeighDist(mol1,mol2,nb1,n1,n2,rem1,rem2)
        if (str(chirality2) != "CHI_UNSPECIFIED") or (bNonsingleBond2==True):
#           print "2",i2,chirality2
            # try fitting locally on the 1-2, 1-3 atoms
            localEnv2 = localEnvironment(mol2,i2,n2)
            if len(localEnv2)>2:
                localEnv1 = getMappedList(n2,n1,localEnv2)
                Chem.rdMolAlign.AlignMol(mol1,mol2,atomMap=list(zip(localEnv1,localEnv2)))
#                alignOnSubset(mol1,mol2,zip(localEnv1,localEnv2))
            # get the neighbours
            nb2 = a2.GetNeighbors()
            # check if the matched neighbours in the aligned molecules are too far
            checkNeighDist(mol2,mol1,nb2,n2,n1,rem2,rem1)

    ####### remove #######
    n1_out = []
    n2_out = []
    for i1,i2 in zip(n1,n2):
        if( (i1 in rem1) or (i2 in rem2) ):
            continue
        n1_out.append(i1)
        n2_out.append(i2)
    return(n1_out,n2_out)

# this module checks for a chirality violation and returns True if smth wrong is found
def bCheckChiralViolation(mol1,mol2,n1,n2):
    bViolation = False
    # if only one atom is in the mapping
    if len(n1)<2 and len(n2)<2:
        return(bViolation)
    # create constraint for alignment
    constrMap = list(zip(n1,n2))
    # align on the subset n1,n2
    Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=list(zip(n2,n1)))
    rem1 = []
    rem2 = []
    bonds1 = mol1.GetBonds()
    bonds2 = mol2.GetBonds()
    # create two dictionaries for mappings [n1] = n2
    dictn1n2 = mappingDict( n1,n2 )
    dictn2n1 = mappingDict( n2,n1 )

    for i1,i2 in zip(n1,n2):
        a1 = mol1.GetAtomWithIdx(i1)
        a2 = mol2.GetAtomWithIdx(i2)
        # checking chirality
        chirality1 = a1.GetChiralTag()
        chirality2 = a2.GetChiralTag()
        # also checking if some bonds are not single
        bonds1 = a1.GetBonds()
        bonds2 = a2.GetBonds()
        bNonsingleBond1 = check_nonsingle_bonds( mol1, mol2, bonds1, n1, n2, dictn1n2 )
        bNonsingleBond2 = check_nonsingle_bonds( mol2, mol1, bonds2, n2, n1, dictn2n1 )
#        print bNonsingleBond1
#        bonds_a1 = find_atom_bonds( i1, a1, bonds1 ) # find all bonds in which a1 takes part
#        sys.exit(0)
        if (str(chirality1) != "CHI_UNSPECIFIED") or (bNonsingleBond1==True):
#           print "1check",i1,chirality1
            localEnv1 = localEnvironment(mol1,i1,n1)
            if len(localEnv1)>2:
                localEnv2 = getMappedList(n1,n2,localEnv1)
                Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=list(zip(localEnv2,localEnv1)))
                #alignOnSubset(mol1,mol2,zip(localEnv1,localEnv2))
            nb1 = a1.GetNeighbors()
            checkNeighDist(mol1,mol2,nb1,n1,n2,rem1,rem2)
        if (str(chirality2) != "CHI_UNSPECIFIED") or (bNonsingleBond2==True):
#           print "2check",i2,chirality2
            localEnv2 = localEnvironment(mol2,i2,n2)
            if len(localEnv2)>2:
                localEnv1 = getMappedList(n2,n1,localEnv2)
                Chem.rdMolAlign.AlignMol(mol1,mol2,atomMap=list(zip(localEnv1,localEnv2)))
                #alignOnSubset(mol1,mol2,zip(localEnv1,localEnv2))
            nb2 = a2.GetNeighbors()
            checkNeighDist(mol2,mol1,nb2,n2,n1,rem2,rem1)
        if len(rem1)+len(rem2)>0:
            return(True)
    return(bViolation)

# bond length table
# http://www.chem.tamu.edu/rgroup/connell/linkfiles/bonds.pdf
# http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
def is_nonsingle( d, anum1, anum2 ):
    if anum1==6 and anum2==6: # C C
        if d < 1.45:
            return(True)
    elif (anum1==6 and anum2==7) or (anum2==6 and anum1==7): # C N
        if d < 1.4:
            return(True)
    elif (anum1==6 and anum2==8) or (anum2==6 and anum1==8): # C O
        if d < 1.35:
            return(True)
    elif (anum1==6 and anum2==16) or (anum2==6 and anum1==16): # C S
        if d < 1.75:
            return(True)
    elif anum1==7 and anum2==7: # N N
        if d < 1.35:
            return(True)
    elif (anum1==7 and anum2==8) or (anum2==7 and anum1==8): # N O
        if d < 1.3:
            return(True)
    elif anum1==8 and anum2==8: # O O
        if d < 1.3:
            return(True)
    elif (anum1==8 and anum2==15) or (anum2==8 and anum1==15): # O P
        if d < 1.575:
            return(True)
    elif anum1==16 and anum2==16: # S S
        if d < 1.75:
            return(True)

    return(False)

def check_nonsingle_bonds( mol1, mol2, bonds1, n1, n2, dictn1n2 ):
    for bond in bonds1:
        foo = bond.GetBeginAtomIdx()
        bar = bond.GetEndAtomIdx()
        if (foo in n1) and (bar in n1): # atoms of the bond considered are both in n1, therefore they are also in n2
            a1 = mol1.GetAtomWithIdx( foo )
            a2 = mol1.GetAtomWithIdx( bar )
            bondLength1 = getBondLength(mol1,foo,bar)
            if is_nonsingle( bondLength1, a1.GetAtomicNum(), a2.GetAtomicNum() )==True:
                return(True)
    return(False)

def mappingDict( n1, n2 ):
    out = {}
    for i1,i2 in zip(n1,n2):
        out[i1] = i2
    return(out)

def checkTopRemove(mol1,mol2,n1,n2,startList,endList):
    rem1 = []
    rem2 = []
    for iStart,iEnd in zip(startList,endList):
        jStart,jEnd = getAttr(n1,n2,iStart,iEnd)
        # count iStart mapped neighbours
        startNeighb = 0
        for b1 in mol1.GetBonds():
            foo = b1.GetBeginAtomIdx()
            bar = b1.GetEndAtomIdx()
            if( iStart==foo ): # atom of interest
                a,b = getAttr(n1,n2,iStart,bar)
                if( (a!=None) and (b!=None) ):
                    startNeighb = startNeighb+1
            elif( iStart==bar ): # atom of interest
                a,b = getAttr(n1,n2,iStart,foo)
                if( (a!=None) and (b!=None) ):
                    startNeighb = startNeighb+1
        # count iEnd mapped neighbour
        endNeighb = 0
        for b1 in mol1.GetBonds():
            foo = b1.GetBeginAtomIdx()
            bar = b1.GetEndAtomIdx()
            if( iEnd==foo ): # atom of interest
                a,b = getAttr(n1,n2,iEnd,bar)
                if( (a!=None) and (b!=None) ):
                    endNeighb = endNeighb+1
            elif( iEnd==bar ): # atom of interest
                a,b = getAttr(n1,n2,iEnd,foo)
                if( (a!=None) and (b!=None) ):
                    endNeighb = endNeighb+1
        # add to remove list
        if( startNeighb < endNeighb ):
            rem1.append(iStart)
            rem2.append(jStart)
        else:
            rem1.append(iEnd)
            rem2.append(jEnd)
    # remove
    for i,j in zip(rem1,rem2):
        n1.remove(i)
        n2.remove(j)
    return(n1,n2)

def getList12(mol,n):
    dict12 = {}
    for a1 in mol.GetAtoms():
        iStart1 = a1.GetIdx()
        if iStart1 not in n:
            continue
        neighbours1 = a1.GetNeighbors()
        for a2 in neighbours1: # 1-2
            iEnd2 = a2.GetIdx()
            if iEnd2 not in n:
                continue
            if iEnd2 == iStart1:
                continue
            if iStart1 in list(dict12.keys()):
                if iEnd2 not in dict12[iStart1]:
                    dict12[iStart1].append(iEnd2)
            else:
                dict12[iStart1] = [iEnd2]
    return(dict12)

def getList13(mol,n):
    dict13 = {}
    for a1 in mol.GetAtoms():
        iStart1 = a1.GetIdx()
        if iStart1 not in n:
            continue
        neighbours1 = a1.GetNeighbors()
        for a2 in neighbours1: # 1-2
            i2 = a2.GetIdx()
#            if i2 not in n:
#                continue
            if i2 == iStart1:
                continue
            neighbours2 = a2.GetNeighbors()
            for a3 in neighbours2: # 1-3
                iEnd3 = a3.GetIdx()
                if iEnd3 not in n:
                    continue
                if (iEnd3==iStart1) or (iEnd3==i2):
                    continue
                if iStart1 in list(dict13.keys()):
                    if iEnd3 not in dict13[iStart1]:
                        dict13[iStart1].append(iEnd3)
                else:
                    dict13[iStart1] = [iEnd3]
    return(dict13)

def getList14(mol,n):
    dict14 = {}
    for a1 in mol.GetAtoms():
        iStart1 = a1.GetIdx()
        if iStart1 not in n:
            continue
        neighbours1 = a1.GetNeighbors()
        for a2 in neighbours1: # 1-2
            i2 = a2.GetIdx()
#            if i2 not in n:
#                continue
            if i2 == iStart1:
                continue
            neighbours2 = a2.GetNeighbors()
            for a3 in neighbours2: # 1-3
                i3 = a3.GetIdx()
#               if i3 not in n:
#                   continue
                if (i3==iStart1) or (i3==i2):
                    continue
                neighbours3 = a3.GetNeighbors()
                for a4 in neighbours3: # 1-4
                    iEnd4 = a4.GetIdx()
                    if iEnd4 not in n:
                        continue
                    if (iEnd4==iStart1) or (iEnd4==i2) or (iEnd4==i3):
                        continue
                    if iStart1 in list(dict14.keys()):
                        if iEnd4 not in dict14[iStart1]:
                            dict14[iStart1].append(iEnd4)
                    else:
                        dict14[iStart1] = [iEnd4]
    return(dict14)

def findProblemsExclusions(n1,n2,dict_mol1,dict_mol2):
    rem_start = []
    rem_end = []
    for iStart in list(dict_mol1.keys()):
        for iEnd in dict_mol1[iStart]:
            jStart,jEnd = getAttr(n1,n2,iStart,iEnd)
            if( (jStart==None) or (jEnd==None) ): # mapped to a dummy, thus no worries
                continue
            if jStart in list(dict_mol2.keys()):
                if jEnd not in dict_mol2[jStart]:
                    # maybe entry already exists
                    if ((jStart in rem_start) or (jStart in rem_end)) and ((jEnd in rem_start) or (jEnd in rem_end)):
                        continue
                    rem_start.append(jStart)
                    rem_end.append(jEnd)
            elif jEnd not in list(dict_mol2.keys()):
                # a weird situation that shouldn't happen
                print("Warning: something wrong in the 1-2, 1-3 or 1-4 lists. Trying to proceed with the warning...")
                rem_start.append(jStart)
                rem_end.append(jEnd)
    return(rem_start,rem_end)

def fixProblemsExclusions(mol1,mol2,n1,n2,startList,endList):
    rem1 = []
    rem2 = []
    for iStart,iEnd in zip(startList,endList):
        jStart,jEnd = getAttr(n1,n2,iStart,iEnd)
        # count iStart mapped neighbours
        startNeighb = 0
        for b1 in mol1.GetBonds():
            foo = b1.GetBeginAtomIdx()
            bar = b1.GetEndAtomIdx()
            if( iStart==foo ): # atom of interest
                a,b = getAttr(n1,n2,iStart,bar)
                if( (a!=None) and (b!=None) ):
                    startNeighb = startNeighb+1
            elif( iStart==bar ): # atom of interest
                a,b = getAttr(n1,n2,iStart,foo)
                if( (a!=None) and (b!=None) ):
                    startNeighb = startNeighb+1
        # count iEnd mapped neighbour
        endNeighb = 0
        for b1 in mol1.GetBonds():
            foo = b1.GetBeginAtomIdx()
            bar = b1.GetEndAtomIdx()
            if( iEnd==foo ): # atom of interest
                a,b = getAttr(n1,n2,iEnd,bar)
                if( (a!=None) and (b!=None) ):
                    endNeighb = endNeighb+1
            elif( iEnd==bar ): # atom of interest
                a,b = getAttr(n1,n2,iEnd,foo)
                if( (a!=None) and (b!=None) ):
                    endNeighb = endNeighb+1
        # add to remove list
        if( startNeighb < endNeighb ):
            rem1.append(iStart)
            rem2.append(jStart)
        else:
            rem1.append(iEnd)
            rem2.append(jEnd)
    # remove
    n1_out = []
    n2_out = []
    for i1,i2 in zip(n1,n2):
        if( (i1 in rem1) or (i2 in rem2) ):
            continue
        n1_out.append(i1)
        n2_out.append(i2)
    return(n1_out,n2_out)

def checkTop(mol1,mol2,n1,n2,bH2H,bH2heavy):
    # 1) generate 1-2, 1-3 and 1-4 lists
    # 2) identify problematic mappings
    # 3) fix the problems: discard the atom with fewer mapped neighbours

    ####### 1-2 #########
    # 1a) 1-2 lists
    dict12_mol1 = getList12(mol1,n1)
    dict12_mol2 = getList12(mol2,n2)
    # 2a) identify problems 1-2; and
    # 3a) fix 1-2
    rem12_mol2_start,rem12_mol2_end = findProblemsExclusions(n1,n2,dict12_mol1,dict12_mol2) # output: indeces of mol2
    n2,n1 = fixProblemsExclusions(mol2,mol1,n2,n1,rem12_mol2_start,rem12_mol2_end)
    rem12_mol1_start,rem12_mol1_end = findProblemsExclusions(n2,n1,dict12_mol2,dict12_mol1) # output: indeces of mol1
    n1,n2 = fixProblemsExclusions(mol1,mol2,n1,n2,rem12_mol1_start,rem12_mol1_end)

    ####### 1-3 #########
    # 1b) 1-3 lists
    dict13_mol1 = getList13(mol1,n1)
    dict13_mol2 = getList13(mol2,n2)
    # 2b) identify problems 1-3 and
    # 3b) fix 1-3
    rem13_mol2_start,rem13_mol2_end = findProblemsExclusions(n1,n2,dict13_mol1,dict13_mol2) # output: indeces of mol2
    n2,n1 = fixProblemsExclusions(mol2,mol1,n2,n1,rem13_mol2_start,rem13_mol2_end)
    rem13_mol1_start,rem13_mol1_end = findProblemsExclusions(n2,n1,dict13_mol2,dict13_mol1) # output: indeces of mol1
    n1,n2 = fixProblemsExclusions(mol1,mol2,n1,n2,rem13_mol1_start,rem13_mol1_end)

    ####### 1-4 #########
    # 1b) 1-4 lists
    dict14_mol1 = getList14(mol1,n1)
    dict14_mol2 = getList14(mol2,n2)
    # 2b) identify problems 1-4 and
    # 3b) fix 1-4
    rem14_mol2_start,rem14_mol2_end = findProblemsExclusions(n1,n2,dict14_mol1,dict14_mol2) # output: indeces of mol2
    n2,n1 = fixProblemsExclusions(mol2,mol1,n2,n1,rem14_mol2_start,rem14_mol2_end)
    rem14_mol1_start,rem14_mol1_end = findProblemsExclusions(n2,n1,dict14_mol2,dict14_mol1) # output: indeces of mol1
    n1,n2 = fixProblemsExclusions(mol1,mol2,n1,n2,rem14_mol1_start,rem14_mol1_end)

    # treat disconnected
 #   n1,n2 = disconnectedMCS(mol1,mol2,n1,n2,bH2H,bH2heavy)
    n1,n2 = disconnectedRecursive(mol1,mol2,n1,n2)
    return(n1,n2)

def checkMCSTop(mol1,mol2,n1,n2,bH2H,bH2heavy):
    # if a bond exists in only one substructure,
    # modify the mapping by discarding one atom participating in the unmatched bond
    # try to discard the atom with fewer mapped neighbours

    # mol1 bonds
    startList = []
    endList = []
    for b1 in mol1.GetBonds():
        iStart = b1.GetBeginAtomIdx()
        iEnd = b1.GetEndAtomIdx()
        jStart,jEnd = getAttr(n1,n2,iStart,iEnd)
        if( (jStart!=None) and (jEnd!=None) ): # not dummies
            bOk = False
            for b2 in mol2.GetBonds():
                foo = b2.GetBeginAtomIdx()
                bar = b2.GetEndAtomIdx()
                if( foo==jStart and bar==jEnd ):
                    bOk = True
                    break
                elif( foo==jEnd and bar==jStart ):
                    bOk = True
                    break
            if(bOk == False):
                startList.append(iStart)
                endList.append(iEnd)
#               return(bOk)
    n1,n2 = checkTopRemove(mol1,mol2,n1,n2,startList,endList,bH2H,bH2heavy)

    # mol2 bonds
    startList = []
    endList = []
    for b1 in mol2.GetBonds():
        iStart = b1.GetBeginAtomIdx()
        iEnd = b1.GetEndAtomIdx()
        jStart,jEnd = getAttr(n2,n1,iStart,iEnd)
        if( (jStart!=None) and (jEnd!=None) ): # not dummies
            bOk = False
            for b2 in mol1.GetBonds():
                foo = b2.GetBeginAtomIdx()
                bar = b2.GetEndAtomIdx()
                if( foo==jStart and bar==jEnd ):
                    bOk = True
                    break
                elif( foo==jEnd and bar==jStart ):
                    bOk = True
                    break
            if(bOk == False):
                startList.append(iStart)
                endList.append(iEnd)
#               return(bOk)
    n1,n2 = checkTopRemove(mol2,mol1,n2,n1,startList,endList)

    return(True)

def writeFormatPDB(fname,m,title="",nr=1):
    fp = open(fname,'w')
    for atom in m.atoms:
        foo = cp.deepcopy(atom)
        # chlorine
        if( 'CL' in atom.name or 'Cl' in atom.name or 'cl' in atom.name ):
            foo.name = "Cl"#+"  "
            print(foo, file=fp)
        # bromine
        elif( 'BR' in atom.name or 'Br' in atom.name or 'br' in atom.name ):
            foo.name = "Br"#+"  "
            print(foo, file=fp)
        elif( len(atom.name) >= 4): # too long atom name
            foo = cp.deepcopy(atom)
            foo.name = foo.name[:3]
            print(foo, file=fp)
        else:
            print(atom, file=fp)
    print('ENDMDL', file=fp)
    fp.close()
#    sys.exit(0)

def reformatPDB(filename,num,randint=42):
    newname = "tempFormat_"+str(randint)+'_'+str(num)+".pdb"
    m = Model().read(filename)

    # adjust atom names and remember the changes
    atomNameID = {}
    sigmaHoleCounter = 1
    sigmaHoleID = []
    for a in m.atoms:
        newAtomName = a.name
        if ('EP' in a.name) or ('LPH' in a.name):
            newAtomName = 'HSH'+str(sigmaHoleCounter)
            sigmaHoleCounter+=1
            sigmaHoleID.append(a.id)
        atomNameID[a.id] = a.name
        a.name = newAtomName

    writeFormatPDB(newname,m)
    return(newname,atomNameID,sigmaHoleID)

def restoreAtomNames(mol,atomNameID):

    for atom in mol.GetAtoms():
        newname = atom.GetMonomerInfo().GetName()
        ind = atom.GetIdx()+1
        oldname = atomNameID[ind]
#        print oldname,newname,len(oldname),len(newname)
#        if (newname.strip() != oldname.strip()) and ('EP' in oldname):
        nametoset = "{:<4}".format(oldname[:4])
        atom.GetMonomerInfo().SetName(nametoset)

def write_pairs(n1,n2,pairsFilename):
    fp = open(pairsFilename,"w")
    for i1,i2 in zip(n1,n2):
        foo = i1 + 1
        bar = i2 + 1
        fp.write("%s    %s\n" % (foo,bar) )
    fp.close()

def calcScore(mol1,mol2,n1,n2,bH2H,bH2heavy):
    res = 0.0
    nn1 = len(n1)
    nn2 = len(n2)
    res = (nn1+nn2)/2.0
    if( bH2H==True or bH2heavy==True): # consider hydrogens
        na1 = mol1.GetNumAtoms()
        na2 = mol2.GetNumAtoms()
    else: # no hydrogens
        na1 = mol1.GetNumHeavyAtoms()
        na2 = mol2.GetNumHeavyAtoms()
    res = 1.0 - res/(na1+na2-res)
    return(res)

def distance_based(mol1, mol2, d, id1=None, id2=None, calcOnly=False):
    pairs1 = []
    pairs2 = []

    # to choose one MCS out of many
    if(calcOnly==True):
        dist = 0.0
        c1 = mol1.GetConformer()
        c2 = mol2.GetConformer()
        for ind1,ind2 in zip(id1,id2):
            pos1 = c1.GetAtomPosition(ind1)
            pos2 = c2.GetAtomPosition(ind2)
            dist = dist + 0.1*pos1.Distance(pos2) # Angstroms in pdb files
        return(dist)

    # o3a
    if(id1==None or id2==None):
        c1 = mol1.GetConformer()
        c2 = mol2.GetConformer()
        for a1 in mol1.GetAtoms():
            pos1 = c1.GetAtomPosition(a1.GetIdx())
            dd = d*10.0 # Angstroms in pdb files
            keep1 = None
            keep2 = None
            for a2 in mol2.GetAtoms():
                pos2 = c2.GetAtomPosition(a2.GetIdx())
                dist = pos1.Distance(pos2)
                if(dist < dd):
                    dd = dist
                    keep1 = a1.GetIdx()
                    keep2 = a2.GetIdx()
            if( (keep1 is not None) and (keep2 is not None) ):
                pairs1.append(keep1)
                pairs2.append(keep2)
        return(pairs1,pairs2)

    # mcs
    for ind1 in id1:
        c1 = mol1.GetConformer()
        pos1 = c1.GetAtomPosition(ind1)
        dd = d*10.0 # Angstroms in pdb files
        keep1 = None
        keep2 = None
        for ind2 in id2:
            c2 = mol2.GetConformer()
            pos2 = c2.GetAtomPosition(ind2)
            dist = pos1.Distance(pos2)
            if(dist < dd):
                dd = dist
                keep1 = ind1
                keep2 = ind2
        if( (keep1 is not None) and (keep2 is not None) ):
            pairs1.append(keep1)
            pairs2.append(keep2)
    return(pairs1,pairs2)

def chargesTypesMMFF(mol):
    mmff = AllChem.MMFFGetMoleculeProperties(mol,mmffVariant='MMFF94')
    return mmff

def crippen(mol):
    crippen = Crippen._GetAtomContribs(mol)
    return crippen

def o3a_alignment(mol1, mol2, bH2H, bH2Hpolar, bH2heavy, bRingsOnly, sigmaHoleID1,sigmaHoleID2, bCrippen=True, d=0.05):
    n1 = []
    n2 = []
####################################
# prepare molecules and parameters #
####################################
    if( bRingsOnly==True ):
        submol1 = subMolRing(mol1)
        submol2 = subMolRing(mol2)
###################
#### now align ####
###################
    if( bCrippen==True ): # always Crippen
        mol1_crippen = crippen(mol1)
        mol2_crippen = crippen(mol2)
        if( bRingsOnly==True ):
            submol1_crippen = crippen(submol1)
            submol2_crippen = crippen(submol2)
            pyO3A = AllChem.GetCrippenO3A(submol1,submol2,submol1_crippen,submol2_crippen,options=0)
            rmsd = pyO3A.Align()
            # now align the full molecule
            pyO3A = AllChem.GetCrippenO3A(mol1,submol1,mol1_crippen,submol1_crippen,options=0)
            pyO3A.Align()
        else:
            pyO3A = AllChem.GetCrippenO3A(mol1,mol2,mol1_crippen,mol2_crippen,options=0) # mol1=probe, mol2=ref
            pyO3A.Align()
        # distances
        n1,n2 = distance_based(mol1,mol2,d)
# hydrogen rule
    n1,n2 = removeH(mol1,mol2,n1,n2,bH2H,bH2Hpolar,bH2heavy)
# remove sigma hole virtual particles
    n1,n2 = remove_sigmaHoles( n1,n2,sigmaHoleID1,sigmaHoleID2)
# remove polar hydrogen mappings
#    n1,n2 = removePolarHmappings(mol1,mol2,n1,n2,bH2Hpolar)
# triple bond rule: for simulation stability do not allow morphing an atom that is involved in a triple bond
# into an atom that is involved in a non-triple bond (and vice versa)
#    n1,n2 = tripleBond(mol1,mol2,n1,n2)
# rings
    n1,n2 = matchRings(mol1,mol2,n1,n2)
# checking possible issues with the 1-2, 1-3 and 1-4 interactions
# this is done before the ringCheck and repeated after
    n1,n2 = checkTop(mol1,mol2,n1,n2,bH2H,bH2heavy)
# do not break rings
    bBreakRings = False
    if( bBreakRings==False ):
        print("Avoiding breaking rings.\n")
        n1,n2 = matchFullRings(mol1,mol2,n1,n2)
#    sys.exit(0)
# treat disconnected
#    n1,n2 = disconnectedMCS(mol1,mol2,n1,n2,bH2H,bH2heavy)
    n1,n2 = disconnectedRecursive(mol1,mol2,n1,n2)
# n1 and n2 are not sorted by pairs at this point
    n1,n2 = sortInd(mol1,mol2,n1,n2)
# checking possible issues with the 1-2, 1-3 and 1-4 interactions
    n1,n2 = checkTop(mol1,mol2,n1,n2,bH2H,bH2heavy)

    return(n1,n2,pyO3A)

# checking that a non-ring atom would not be morphed into a ring atom
def matchRings(mol1,mol2,nfoo,nbar):
    newn1 = []
    newn2 = []
    for n1,n2 in zip(nfoo,nbar):
        a1 = mol1.GetAtomWithIdx(n1)
        a2 = mol2.GetAtomWithIdx(n2)
#        arom1 = a1.GetIsAromatic()
#        arom2 = a2.GetIsAromatic()
        ring1 = a1.IsInRing()
        ring2 = a2.IsInRing()
        if(ring1==True and ring2==False):
            continue
        if(ring1==False and ring2==True):
            continue
        newn1.append(n1)
        newn2.append(n2)

    # only one atom morphed in a ring will not work, check for that
    n1,n2 = oneAtomInRing(mol1,mol2,newn1,newn2)
    return(n1,n2)

def oneAtomInRing(mol1,mol2,n1,n2):
    newn1 = []
    newn2 = []
    for i,j in zip(n1,n2):
        a1 = mol1.GetAtomWithIdx(i)
        a2 = mol2.GetAtomWithIdx(j)
        ring1 = a1.IsInRing()
        ring2 = a2.IsInRing()
        if( ring1==True and ring2==True ):
            bonds1 = a1.GetBonds()
            bonds2 = a2.GetBonds()
            found = 0
            for b1 in bonds1:
                id1 = b1.GetEndAtomIdx()
                at1 = b1.GetEndAtom()
                if( b1.GetEndAtomIdx()==i ):
                    id1 = b1.GetBeginAtomIdx()
                    at1 = b1.GetBeginAtom()
                for b2 in bonds2:
                    id2 = b2.GetEndAtomIdx()
                    at2 = b2.GetEndAtom()
                    if( b2.GetEndAtomIdx()==j ):
                        id2 = b2.GetBeginAtomIdx()
                        at2 = b2.GetBeginAtom()
                    if(at1.IsInRing()==True and at2.IsInRing()==True):
                        if( (id1 in n1) and (id2 in n2) ):
                            found = 1
                            break
                if(found==1):
                    break
            if(found==1):
                newn1.append(i)
                newn2.append(j)
        else:
            newn1.append(i)
            newn2.append(j)
    return(newn1,newn2)

def carbonize(mol, bH2heavy=False, bRingsOnly=None):
    for atom in mol.GetAtoms():
        if(atom.GetAtomicNum() != 1):
            atom.SetAtomicNum(6)
        elif(bH2heavy == True):
            atom.SetAtomicNum(6)
        if( (bRingsOnly!=None) and (atom.IsInRing()==False) ):
            atom.SetAtomicNum(bRingsOnly)

def carbonizeCrippen(mol, bH2heavy=False, bRingsOnly=None):
    crippen = []
    for atom in mol.GetAtoms():
        if(atom.GetAtomicNum() != 1):
            foo = (0.1441,2.503)
            crippen.append(foo)
        elif(bH2heavy == True):
            foo = (0.1441,2.503)
            crippen.append(foo)
        if( (bRingsOnly!=None) and (atom.IsInRing()==False) ):
            foo = (bRingsOnly*(-1),bRingsOnly)
            crippen.append(foo)
    return(crippen)

def getBondLength(mol,id1,id2):
    conf = mol.GetConformer()
    pos1 = conf.GetAtomPosition(id1)
    pos2 = conf.GetAtomPosition(id2)
    return(pos1.Distance(pos2))

def isTriple(mol,a):
    bTriple = False
    neighb = a.GetNeighbors()
    # analyze C,N (S not considered, because somewhat exotic) atoms for triple bonds
    # C
    if( a.GetAtomicNum()==6 ):
        if( len(a.GetNeighbors())==2 ):
            for neighb in a.GetNeighbors():
                if( neighb.GetAtomicNum()==7 ):
                    if( len(neighb.GetNeighbors())==1 ):
                        bTriple=True
                if( neighb.GetAtomicNum()==6 ):
                    if( len(neighb.GetNeighbors())==2 ):
                        if( getBondLength(mol,a.GetIdx(),neighb.GetIdx())<1.25 ): # need to check bond length (in Angstroms)
                            bTriple=True
    # N
    elif( a.GetAtomicNum()==7 ):
        if( len(a.GetNeighbors())==1 ):
            bTriple=True

   #    Chem.MolToMolFile(mol1,"foomol.mol")
#    foo = Chem.MolFromMolFile("foomol.mol")
#    Chem.MolToMolFile(mol2,"barmol.mol")
#    bar = Chem.MolFromMolFile("barmol.mol")
#    os.remove("foomol.mol")
#    os.remove("barmol.mol")

    return(bTriple)

def isPolarH( a ):
    neighb = a.GetNeighbors()
    for a2 in neighb:
        anum2 = a2.GetAtomicNum()
        if anum2!=6:
            return(True)
    return(False)

def removePolarHmappings(mol1,mol2,nfoo,nbar,bH2Hpolar):
    if bH2Hpolar==True:
        return(nfoo,nbar)

    newn1 = []
    newn2 = []
    for n1,n2 in zip(nfoo,nbar):
        a1 = mol1.GetAtomWithIdx(n1)
        a2 = mol2.GetAtomWithIdx(n2)
        anum1 = a1.GetAtomicNum()
        anum2 = a2.GetAtomicNum()

        # remove polar H mappings
        bPolar1 = False
        bPolar2 = False
        if anum1==1:
            bPolar1 = isPolarH( a1 )
        if anum2==1:
            bPolar2 = isPolarH( a2 )

        if(bPolar1==True or bPolar2==True):
            continue

        newn1.append(n1)
        newn2.append(n2)

    return(newn1,newn2)


def tripleBond(mol1,mol2,nfoo,nbar):
    newn1 = []
    newn2 = []
    for n1,n2 in zip(nfoo,nbar):
        a1 = mol1.GetAtomWithIdx(n1)
        a2 = mol2.GetAtomWithIdx(n2)
        bTriple1 = False
        bTriple2 = False
        # identify if bTriple is True/False
        bTriple1 = isTriple(mol1,a1)
        bTriple2 = isTriple(mol2,a2)
        if(bTriple1==True and bTriple2==False):
            continue
        elif(bTriple2==True and bTriple1==False):
            continue
        newn1.append(n1)
        newn2.append(n2)
    return(newn1,newn2)

def removeH(mol1,mol2,nfoo,nbar,bH2H,bH2Hpolar,bH2heavy):
    newn1 = []
    newn2 = []
    for n1,n2 in zip(nfoo,nbar):
        a1 = mol1.GetAtomWithIdx(n1)
        a2 = mol2.GetAtomWithIdx(n2)
        id1 = a1.GetAtomicNum()
        id2 = a2.GetAtomicNum()
        bPolar1 = False
        bPolar2 = False
        if( id1==1 and id2==1):
            bPolar1 = isPolarH( a1 )
            bPolar2 = isPolarH( a2 )
            if(bPolar1==True and bPolar2==True):
                if( bH2Hpolar==False ):
                    continue
            elif( bH2H==False ):
                continue
        elif(bH2heavy==False and ( (id1==1) ^ (id2==1) ) ): # ^ := xor
            continue
        newn1.append(n1)
        newn2.append(n2)
    return(newn1,newn2)

def subMolByIndex(mol,ind):
    copyMol = copy.deepcopy(mol)
    editMol = Chem.EditableMol(copyMol)
    indRm = []
#   create an inverted list of ind, i.e. indRm
    for a in mol.GetAtoms():
        found = 0
        for i in ind:
            if(i == a.GetIdx()):
                found = 1
                break
        if(found == 0):
            indRm.append(a.GetIdx())
#   remove the indRm atoms
    indRm.sort(reverse=True)
    for i in indRm:
        editMol.RemoveAtom(i)
    copyMol = editMol.GetMol()
    return(copyMol)

def subMolRing(mol):
    copyMol = copy.deepcopy(mol)
    editMol = Chem.EditableMol(copyMol)
    indRm = []
#   create an inverted list of ind, i.e. indRm
    for a in mol.GetAtoms():
        if(a.IsInRing()==False):
            indRm.append(a.GetIdx())
#   remove the indRm atoms
    indRm.sort(reverse=True)
    for i in indRm:
        editMol.RemoveAtom(i)
    copyMol = editMol.GetMol()
    return(copyMol)

def calcRMSD(mol1,mol2,ind1,ind2):
    rmsd = 0.0
    c1 = mol1.GetConformer()
    c2 = mol2.GetConformer()
    for id1 in ind1:
        pos1 = c1.GetAtomPosition(id1)
        toAdd = 999.999
        for id2 in ind2:
            pos2 = c2.GetAtomPosition(id2)
            if(pos1.Distance(pos2)<toAdd):
                toAdd = pos1.Distance(pos2)
        rmsd = rmsd + pow(toAdd,2)
    rmsd = sqrt(rmsd)
    return(rmsd)

def matchIDbyRMSD(subMol,ind,mol):
    res_ind = []
    subc = subMol.GetConformer()
    c = mol.GetConformer()
    for i in ind:
        pos1 = subc.GetAtomPosition(i)
        rmsd = 999.999
        keep = -1
        for atom in mol.GetAtoms():
            pos2 = c.GetAtomPosition(atom.GetIdx())
            if(pos1.Distance(pos2)<rmsd):
                rmsd = pos1.Distance(pos2)
                keep = atom.GetIdx()
        res_ind.append(keep)
    return(res_ind)


def disconnectedMCS(mol1,mol2,ind1,ind2,bH2H=True,bH2Hpolar=True,bH2heavy=True):
    subMol1 = subMolByIndex(mol1,ind1)
    subMol2 = subMolByIndex(mol2,ind2)
    res = MCS.FindMCS([subMol1,subMol2],ringMatchesRingOnly=True, completeRingsOnly=True, atomCompare='any', bondCompare='any')
    p = Chem.FragmentMatcher
    pp = p.FragmentMatcher()
    n1_orig = []
    n2_orig = []
    try:
        pp.Init(res.smarts)
    except:
        if len(ind1)>1:
            print("WARNING: the mapping may (but not necessarily) contain disconnected fragments. Proceed with caution.")
        return(ind1,ind2)
#       return(n1_orig,n2_orig)
    n1_list = pp.GetMatches(subMol1)
    n2_list = pp.GetMatches(subMol2)
# remove all matched hydrogens accordingly
    n1_list,n2_list = mcsHremove(subMol1,subMol2,n1_list,n2_list,bH2H,bH2Hpolar,bH2heavy)
# find out which of the generated list pairs has the smallest rmsd
    minRMSD = 999999.99
    for nl1 in n1_list:
        for nl2 in n2_list:
            rmsd = calcRMSD(subMol1,subMol2,nl1,nl2)
            if(rmsd < minRMSD):
                minRMSD = rmsd
                n1 = nl1
                n2 = nl2
# match indices n1,n2 to the original molecule ind1,ind2
    n1_orig = matchIDbyRMSD(subMol1,n1,mol1)
    n2_orig = matchIDbyRMSD(subMol2,n2,mol2)
    return(n1_orig,n2_orig)

def recursiveFragmentSearch(mol,atom,fragment,n_list):
    ind = atom.GetIdx()
    if (ind not in fragment) and (ind in n_list):
        fragment.append(ind)
        a = mol.GetAtomWithIdx(ind)
        nb = a.GetNeighbors()
        for neigh in nb:
            recursiveFragmentSearch(mol,neigh,fragment,n_list)

def disconnectedRecursive(mol1,mol2,ind1,ind2):
###################################
### find disconnected fragments ###
    # molecule 1
    n1_newlist = []
    n1_fragments = []
    for ind in ind1:
        if ind in n1_newlist:
            continue
        fragment = []
        atom = mol1.GetAtomWithIdx(ind)
        recursiveFragmentSearch(mol1,atom,fragment,ind1)
        n1_newlist.extend(fragment)
        n1_fragments.append(fragment)
    # molecule 2
    n2_newlist = []
    n2_fragments = []
    for ind in ind2:
        if ind in n2_newlist:
            continue
        fragment = []
        atom = mol2.GetAtomWithIdx(ind)
        recursiveFragmentSearch(mol2,atom,fragment,ind2)
        n2_newlist.extend(fragment)
        n2_fragments.append(fragment)
###########################################
#### identify all the fragment matches ####
# for testing purposes
#    n2_fragments = [[0, 25, 24, 23, 22, 21, 20, 19, 14], [13, 12, 3, 2, 1, 27, 28, 5, 6, 7, 9, 11, 10, 8, 33, 35, 36, 34, 32, 30, 31, 29, 4, 18, 17, 16, 38, 39, 40, 41, 42, 43, 37, 15, 44, 45, 46, 26, 50, 47, 48, 49], [51]]
    matchDict1 = {}
    matchDict2 = {}
    for id1,id2 in zip(ind1,ind2):
        # find id1
        key1 = ''
        for i in range(0,len(n1_fragments)):
            if id1 in n1_fragments[i]:
                key1 =  str(i)
                break
        # find id2
        key2 = ''
        for i in range(0,len(n2_fragments)):
            if id2 in n2_fragments[i]:
                key2 =  str(i)
                break
        key = key1+'_'+key2
        if key in list(matchDict1.keys()):
            matchDict1[key].append(id1)
            matchDict2[key].append(id2)
        else:
            matchDict1[key] = [id1]
            matchDict2[key] = [id2]
##################################
#### find the largest matches ####
    maxMatchSize = -1
    minMatchRMSD = 99999.999
    maxMatchKey = ''
    for key in matchDict1:
        if len(matchDict1[key]) > maxMatchSize:
            maxMatchSize = len(matchDict1[key])
            maxMatchKey = key
            minMatchRMSD = Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=list(zip(matchDict2[key],matchDict1[key])))
        elif len(matchDict1[key]) == maxMatchSize:
            rmsd = Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=list(zip(matchDict2[key],matchDict1[key])))
            if rmsd < minMatchRMSD:
                minMatchRMSD = rmsd
                maxMatchKey = key
#########################
######## output #########
    if maxMatchKey == '':
        return([],[])
    return(matchDict1[maxMatchKey],matchDict2[maxMatchKey])

def sortInd(mol1,mol2,ind1,ind2):
    n1out = []
    n2out = []
    c1 = mol1.GetConformer()
    c2 = mol2.GetConformer()
    for id1 in ind1:
        pos1 = c1.GetAtomPosition(id1)
        minDist = 9999.999
        keep = -1
        for id2 in ind2:
            pos2 = c2.GetAtomPosition(id2)
            if(pos1.Distance(pos2)<minDist):
                minDist = pos1.Distance(pos2)
                keep = id2
        n1out.append(id1)
        n2out.append(keep)
    return(n1out,n2out)

def disconnected(mol,ind):
#   create a molecule of the ind atoms only
    copyMol = copy.deepcopy(mol)
    editMol = Chem.EditableMol(copyMol)
    indRm = []
#   create an inverted list of ind, i.e. indRm
    for a in mol.GetAtoms():
        found = 0
        for i in ind:
            if(i == a.GetIdx()):
                found = 1
                break
        if(found == 0):
            indRm.append(a.GetIdx())
#   remove the indRm atoms
#   print len(indRm),len(ind),copyMol.GetNumAtoms()
    indRm.sort(reverse=True)
    for i in indRm:
        editMol.RemoveAtom(i)
    copyMol = editMol.GetMol()
#   get the disconnected fragments
    indLists = Chem.GetMolFrags(copyMol,asMols=False,sanitizeFrags=False)
#   only keep the largest lists
#    maxSize = getLargestList(indLists)
#    idLists = []
#    for l in indLists:
#       if(len(l) == maxSize):
#           idLists.append(l)
#   match the IDs of the idLists to the original ind
    resID = []
    for l in indLists:
        ll = [ind[i] for i in l]
        resID.append(ll)
    return(resID)

def getLargestList(lists):
    res = 0
    for l in lists:
        if(len(l) > res):
            res = len(l)
    return(res)

def genFilename(filename,counter):
    if(counter == 0):
        name = filename
    else:
        name = os.path.splitext(filename)[0]+"_"+str(counter)+os.path.splitext(filename)[1]
    return(name)

def incrementByOne(foo):
    bar = []
    for l in foo:
        l = l+1
        bar.append(l)
    return(bar)

def mcsHremove(mol1,mol2,n1_list,n2_list,bH2H,bH2Hpolar,bH2heavy):
    n1 = []
    n2 = []
    for nfoo in n1_list:
        for nbar in n2_list:
            foo,bar = removeH(mol1,mol2,nfoo,nbar,bH2H,bH2Hpolar,bH2heavy)
            n1.append(foo)
            n2.append(bar)
    return(n1,n2)

def mcsDist(mol1,mol2,n1_list,n2_list,d,bH2H,bH2heavy):
    n1 = []
    n2 = []
    # distances
    maxMCS = 0 # size of the largest MCS fulfilling the distances
    for nfoo,nbar in zip(n1_list,n2_list):
        alignID = list(zip(nfoo,nbar))
##########################################
###### o3a alignment may work better #####
        rmsd = Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=list(zip(nbar,nfoo)))
#       rmsd = alignOnSubset(mol1,mol2,alignID) # but it has some dependence on the molecule sequence, not sure if I trust it
#        print "RMSD after alignment: %f Angstroms" %rmsd
        x,y = distance_based(mol1,mol2,d,nfoo,nbar)
        n1.append(x)
        n2.append(y)
        if(len(x)>maxMCS):
            maxMCS = len(x)

    nn1 = []
    nn2 = []
    # match rings and remove disconnected
    maxMCS = 0
    for nfoo,nbar in zip(n1,n2):
        # rings
        x,y = matchRings(mol1,mol2,nfoo,nbar)
        # disconnected
        #x,y = disconnectedMCS(mol1,mol2,x,y,bH2H,bH2heavy)
        x,y = disconnectedRecursive(mol1,mol2,x,y)
        nn1.append(x)
        nn2.append(y)
        if(len(x)>maxMCS):
            maxMCS = len(x)

    print("maxMCS after distance treatment: %d" % maxMCS)

    n1 = []
    n2 = []
    # only keep the largest MCSs
    maxSize = getLargestList(nn1+nn2)
    for nfoo,nbar in zip(nn1,nn2):
        if(len(nfoo)==maxSize and len(nbar)==maxSize):
            n1.append(nfoo)
            n2.append(nbar)
#           print "foo",nfoo,nbar
    return(n1,n2)

def selectOneMCS(n1_list,n2_list,mol1,mol2):
    n1 = n1_list[0]
    n2 = n2_list[0]

    # select the largest MCSs
    largestSize = -42
    for nfoo,nbar in zip(n1_list,n2_list):
        if len(nfoo)>largestSize:
            largestSize = len(nfoo)
    n1_largest = []
    n2_largest = []
    for nfoo,nbar in zip(n1_list,n2_list):
        if len(nfoo)==largestSize:
            n1_largest.append(nfoo)
            n2_largest.append(nbar)

    # of the largest ones select the minRMSD
    rmsdMin = 9999.999
    for nfoo,nbar in zip(n1_largest,n2_largest):
        # align
        alignID = list(zip(nbar,nfoo))
        try:
            rmsd = Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=alignID)
        except:
            rmsd = rmsdMin*10.0
#            x,y = distance_based(mol1,mol2,d,nfoo,nbar,True)
            # compare
        if( rmsd < rmsdMin ):
            rmsdMin = rmsd
            n1 = nfoo
            n2 = nbar
    return(n1,n2)

def checkSingleAtomInRingMap( rem1,n1,n2,mol1,mol2,r1,r2,dontRem1,dontRem2 ):

    for i,j in zip(n1,n2):
        if i in rem1:
            a1 = mol1.GetAtomWithIdx(i)
            a2 = mol2.GetAtomWithIdx(j)
            ring1 = a1.IsInRing()
            ring2 = a2.IsInRing()
            identifiedRings = []
            if ring1==False and ring2==True:
                for ar2 in r2:#.AtomRings(): # go over the rings
                    for at2 in ar2: # go over atoms
#                        print "vg ",i,ar2,at2,j
                        if at2==j: # found a ring of interest
#                            print "vg FOUND"
                            identifiedRings.append(ar2)
                            break

#            print "VG ir",identifiedRings
            # now go over all the rings, where atom a2 is found
            # check if there is any other mapping between the atoms of those rings
            bRem = False
            for ring2 in identifiedRings:
#                print "vg ",r2
                bFoo,mappedInd = countMapped( ring2,n2 )
#               print "vg ",bFoo,mappedInd,len(mappedInd)
                if len(mappedInd)>1:
 #                   print "HERE"
                    bRem = True
                    break

            if bRem==False: # do not need to remove these atoms
                dontRem1.append(i)
                dontRem2.append(j)


def matchFullRings(mol1,mol2,n1,n2):
    r1 = mol1.GetRingInfo()
    r2 = mol2.GetRingInfo()
    if isinstance(r1,tuple)==False:
        r1 = r1.AtomRings()
    if isinstance(r2,tuple)==False:
        r2 = r2.AtomRings()
    rem1 = []
    rem2 = []
    # investigate the rings: soft-check
#    print "Soft-check for ring matching\n"
    for i,j in zip(n1,n2):
        a1 = mol1.GetAtomWithIdx(i)
        a2 = mol2.GetAtomWithIdx(j)
        ring1 = a1.IsInRing()
        ring2 = a2.IsInRing()
        if( (ring1==True) and (ring2==False) ):
            rem1.append(i)
            rem2.append(j)
        elif( (ring1==False) and (ring2==True) ):
            rem1.append(i)
            rem2.append(j)
        elif( (ring1==True) and (ring2==True) ):
            mapped1 = False
            mapped2 = False
            # here only checking if, given one morphable atom in a ring,
            # is there at least one ring which would harbor this atom
            # and all the other atoms in that ring would also be morphable
            for ar1 in r1:#.AtomRings():
                if( i in ar1 ):
                    mapped1 = isMapped(ar1,n1)
                if( mapped1 == True):
                    break
            for ar2 in r2:#.AtomRings():
                if( j in ar2 ):
                    mapped2 = isMapped(ar2,n2)
                if( mapped2 == True):
                    break
            if( (mapped1==False) or (mapped2==False) ):
                rem1.append(i)
                rem2.append(j)

    # before removing, check for a special case:
    # a single non-ring atom could be allowed to map to a single ring atom
    # if the rest of the ring atoms are not mapped
    dontRem1 = []
    dontRem2 = []
    checkSingleAtomInRingMap( rem1,n1,n2,mol1,mol2,r1,r2,dontRem1,dontRem2 )
    checkSingleAtomInRingMap( rem2,n2,n1,mol2,mol1,r2,r1,dontRem2,dontRem1 )

#    print "VG ",dontRem1
#    print "VG ",dontRem2

    # remove
    n1_out = []
    n2_out = []
    for i1,i2 in zip(n1,n2):
        if( (i1 in rem1) or (i2 in rem2) ):
            if( (i1 not in dontRem1) and (i2 not in dontRem2) ):
                continue
        n1_out.append(i1)
        n2_out.append(i2)

    # investigate the rings: strict check
    # after the soft-check only those morphable atoms have survived
    # which belong to rings in both structures
    # check if more than two atoms in a ring are morphed, while the rest are not
    bFound = True
    while( bFound==True ):
        n1 = cp.deepcopy(n1_out)
        n2 = cp.deepcopy(n2_out)
        bRem1 = False
        bRem2 = False
        minRemSize = 999
        minRem = []

        # go over the rings in the first molecule
        for ar1 in r1:#.AtomRings():
            bRem1,mappedInd = countMapped(ar1,n1)
#            print ar1,mappedInd,bRem1
            # if an inappropriate mapping is found
            if( bRem1 == True ):
                # identify the other ring that participates in this mapping
                # more precisely, find the smallest number of atoms to remove
                for ar in r1:#.AtomRings():
                    if ar==ar1:
                        continue
                    bAddThisRing = False
                    for a in ar:
                        if a in mappedInd:
                            bAddThisRing = True
                            break
                    foo = []
                    if bAddThisRing==True:
                        for a in ar:
                            foo.append(a)
                    if len(foo)>0 and len(foo)<minRemSize:
                        minRemSize = len(foo)
                        minRem = cp.deepcopy(foo)
                # identify the counterpart indices
                minRemB = []
                for i1,i2 in zip(n1,n2):
                    if i1 in minRem:
                        minRemB.append(i2)
                # remove the atoms
                n1_out,n2_out = removeInd( n1,n2,minRem,minRemB )
#                print n1,minRem
                break
        if minRemSize>0 and minRemSize<999:
            continue
#        sys.exit(0)

        # go over the rings in the second molecule
        for ar2 in r2:#.AtomRings():
            bRem2,mappedInd = countMapped(ar2,n2)
            # if an inappropriate mapping is found
            if( bRem2 == True ):
                # identify the other ring that participates in this mapping
                # more precisely, find the smallest number of atoms to remove
                for ar in r2:#.AtomRings():
                    if ar==ar2:
                        continue
                    bAddThisRing = False
                    for a in ar:
                        if a in mappedInd:
                            bAddThisRing = True
                            break
                    foo = []
                    if bAddThisRing==True:
                        for a in ar:
                            foo.append(a)
                    if len(foo)>0 and len(foo)<minRemSize:
                        minRemSize = len(foo)
                        minRem = cp.deepcopy(foo)
                # identify the counterpart indices
                minRemB = []
                for i1,i2 in zip(n1,n2):
                    if i2 in minRem:
                        minRemB.append(i1)
                # remove the atoms
                n1_out,n2_out = removeInd( n1,n2,minRemB,minRem )
                break
        if minRemSize>0 and minRemSize<999:
            continue

        bFound = False

    return(n1_out,n2_out)

def removeInd( n1,n2,rem1,rem2 ):
    n1_out = []
    n2_out = []
    for i1,i2 in zip(n1,n2):
        if( (i1 in rem1) or (i2 in rem2) ):
            continue
        n1_out.append(i1)
        n2_out.append(i2)
    return(n1_out,n2_out)

def countMapped(ring,ind):
    countMapped = 0
    countUnmapped = 0
    mapped = []
    for a in ring:
        if( a in ind ):
            countMapped += 1
            mapped.append(a)
        else:
            countUnmapped += 1
#    print ring,countMapped,countUnmapped,ind
    if (countMapped > 2) and (countUnmapped > 1):
        return(True,mapped)
    return(False,mapped)

def isMapped(ring,ind):
    for a in ring:
        if( a in ind):
            continue
        else:
            return(False)
    return(True)

def mcs(mol1, mol2, molForMcs1,molForMcs2, bH2H, bH2Hpolar, bH2heavy, bdMCS, bRingsOnly, d, bChiral, sigmaHoleID1, sigmaHoleID2, t=None):
    # make all atoms into carbon
    foo = copy.deepcopy(molForMcs1)
    bar = copy.deepcopy(molForMcs2)
    if( bRingsOnly==True ):
        carbonize(foo,bH2heavy,42)
        carbonize(bar,bH2heavy,43)
    else:
        carbonize(foo,bH2heavy)
        carbonize(bar,bH2heavy)
    mols = [foo,bar]
    print("Searching...")
    res = MCS.FindMCS(mols,ringMatchesRingOnly=True, completeRingsOnly=True, atomCompare='elements', bondCompare='any', timeout=int(t), maximize='bonds')
    # for new RDKit-2018 use below
#    res = rdFMCS.FindMCS(mols,ringMatchesRingOnly=True, completeRingsOnly=True, timeout=int(t), maximizeBonds=True, bondCompare=rdFMCS.BondCompare.CompareAny, atomCompare=rdFMCS.AtomCompare.CompareElements)
    p = Chem.FragmentMatcher
    pp = p.FragmentMatcher()
    n1_list = []
    n2_list = []
    try:
        pp.Init(res.smarts)
    except:
        return(n1_list,n2_list)
    n1_list = pp.GetMatches(foo)
    n2_list = pp.GetMatches(bar)
    print('Found %d MCSs in total (mol1: %d, mol2: %d), each with %d atoms and %d bonds' % (len(n1_list)*len(n2_list),len(n1_list),len(n2_list),res.numAtoms,res.numBonds))
    # if hydrogens to be removed
    n1_list,n2_list = mcsHremove(mol1,mol2,n1_list,n2_list,bH2H,bH2Hpolar,bH2heavy)
    # from this point n1_list and n2_list elements must match 1to1, i.e. the number of elements in the lists is the same

# triple bond rule: for simulation stability do not allow morphing an atom that is involved in a triple bond
# into an atom that is involved in a non-triple bond (and vice versa)
# also checking possible issues with the 1-2, 1-3 and 1-4 interactions
    n1_foo = []
    n2_foo = []
    for n1,n2 in zip(n1_list,n2_list):
#        n1,n2 = removePolarHmappings(mol1,mol2,n1,n2,bH2Hpolar)
#        n1,n2 = tripleBond(mol1,mol2,n1,n2)
        n1,n2 = checkTop(mol1,mol2,n1,n2,bH2H,bH2heavy)
        #n1,n2 = disconnectedMCS(mol1,mol2,n1,n2,bH2H,bH2heavy)
        n1,n2 = disconnectedRecursive(mol1,mol2,n1,n2)
        n1_foo.append(n1)
        n2_foo.append(n2)
    n1_list = cp.copy(n1_foo)
    n2_list = cp.copy(n2_foo)


# test
#    for n1,n2 in zip(n1_list,n2_list):
#        print "foo"
#       for i1,i2 in zip(n1,n2):
#           print i1+1,i2+1
#       print "\n";
#    sys.exit(0)
#############################################
######### chirality check ###################
#    bChiral = True#False#True
    if( bChiral==True ):
        print("Chirality check.")
        n1_foo = []
        n2_foo = []
        for n1,n2 in zip(n1_list,n2_list):
            while(bCheckChiralViolation(mol1,mol2,n1,n2)==True):
                n1,n2 = checkChiral(mol1,mol2,n1,n2)
                n1,n2 = disconnectedRecursive(mol1,mol2,n1,n2)
            n1_foo.append(n1)
            n2_foo.append(n2)
        n1_list = cp.copy(n1_foo)
        n2_list = cp.copy(n2_foo)

# checking possible issues with the 1-2, 1-3 and 1-4 interactions
# this is done before the ringCheck and repeated after
    n1_foo = []
    n2_foo = []
    for n1,n2 in zip(n1_list,n2_list):
        n1,n2 = checkTop(mol1,mol2,n1,n2,bH2H,bH2heavy)
        n1_foo.append(n1)
        n2_foo.append(n2)
    n1_list = cp.copy(n1_foo)
    n2_list = cp.copy(n2_foo)


#############################################
######### do not break rings ################
    bBreakRings = False
    if( bBreakRings==False ):
        print("Avoiding breaking rings.")
        n1_foo = []
        n2_foo = []
        for n1,n2 in zip(n1_list,n2_list):
            n1,n2 = matchFullRings(mol1,mol2,n1,n2)
            n1,n2 = disconnectedRecursive(mol1,mol2,n1,n2)
            n1_foo.append(n1)
            n2_foo.append(n2)
        n1_list = cp.copy(n1_foo)
        n2_list = cp.copy(n2_foo)
# test
#    for n1,n2 in zip(n1_list,n2_list):
#        print "foo"
#       for i1,i2 in zip(n1,n2):
#           print i1+1,i2+1
#       print "\n";
#    sys.exit(0)


    # if distances to be compared
    if(bdMCS==True):
        n1_list,n2_list = mcsDist(mol1,mol2,n1_list,n2_list,d,bH2H,bH2heavy)
        # due to meeting distance criterium
        # the rings may be broken
        # and disconnected fragments may appear
        bBreakRings = False
        if( bBreakRings==False ):
            print("Avoiding breaking rings after meeting distance criterium.\n")
            n1_foo = []
            n2_foo = []
            for n1,n2 in zip(n1_list,n2_list):
                n1,n2 = matchFullRings(mol1,mol2,n1,n2)
                n1,n2 = disconnectedRecursive(mol1,mol2,n1,n2)
                n1_foo.append(n1)
                n2_foo.append(n2)
            n1_list = cp.copy(n1_foo)
            n2_list = cp.copy(n2_foo)

    # if there are several MCSs, select the one yielding the smallest RMSD
    n1,n2 = selectOneMCS(n1_list,n2_list,mol1,mol2)

    # one more final check for possible issues with the 1-2, 1-3 and 1-4 interactions
    n1,n2 = checkTop(mol1,mol2,n1,n2,bH2H,bH2heavy)

    # remove sigma hole virtual particles
    n1,n2 = remove_sigmaHoles( n1,n2,sigmaHoleID1,sigmaHoleID2)

    print('Final MCS that survived after pruning: %d atoms' % (len(n1)))

    return n1,n2

def checkRingsOnlyFlag(mol1,mol2):
    flag = 0
    for atom in mol1.GetAtoms():
        if(atom.IsInRing()==True):
            flag = flag + 1
            break
    for atom in mol2.GetAtoms():
        if(atom.IsInRing()==True):
            flag = flag + 1
            break
    if(flag==2):
        return(True)
    else:
        return(False)

def main(argv):

    desc=("Provided two structures find atoms to be morphed.")

# define input/output files

    files= [
        FileOption("-i1", "r",["pdb"],"input1.pdb", "input"),
        FileOption("-i2", "r",["pdb"],"input2.pdb", "input"),
        FileOption("-o", "w",["dat"],"pairs.dat", "output"),
        FileOption("-opdb1", "w/o",["pdb"],"out_pdb1.pdb", "optional output of superimposed structure 1"),
        FileOption("-opdb2", "w/o",["pdb"],"out_pdb2.pdb", "optional output of superimposed structure 2"),
        FileOption("-opdbm1", "w/o",["pdb"],"out_pdb_morphe1.pdb", "optional output of the morphable atoms str1"),
        FileOption("-opdbm2", "w/o",["pdb"],"out_pdb_morphe2.pdb", "optional output of the morphable atoms str2"),
        FileOption("-score", "w/o",["dat"],"out_score.dat", "optional output: score of the morph"),
        FileOption("-n1", "r/o",["ndx"],"scaffold1" ,"optionally read index of atoms to consider: mol1"),
        FileOption("-n2", "r/o",["ndx"],"scaffold2","optionally read index of atoms to consider: mol2" ),
        ]

# define options

    options=[
        Option( "-alignment", "bool", "True", "method 1: 3D alignment"),
        Option( "-mcs", "bool", "True", "method 2: maximum common substructure"),
        Option( "-H2H", "bool", "False", "should non-polar hydrogen be morphed into any other hydrogen?"),
        Option( "-H2Hpolar", "bool", "False", "should polar hydrogens be morphed into polar hydrogens?"),
        Option( "-H2heavy", "bool", "False", "should hydrogen be morphed into a heavy atom (also sets -H2H to true)"),
        Option( "-RingsOnly", "bool", "False", "should rings only be used in the MCS search or alignment"),
        Option( "-dMCS", "bool", "False", "find MCS, superimpose the structures based on the MCS, apply distance in\
                Cartesian coordinate space to define the morphes"),
        Option( "-chirality", "bool", "True", "perform chirality check for MCS mapping"),
        Option( "-d", "float", "0.05", "distance (nm) between atoms to consider them morphable"),
        Option( "-timeout", "float", "None", "maximum time (s) for an MCS search"),
        ]

    help_text = ()

# pass options, files and the command line to pymacs

    cmdl = Commandline( argv, options = options, fileoptions = files, program_desc = help_text, check_for_existing_files = False, version = "0.0" )

    # deal with the flags
    bH2H = False
    bH2Hpolar = False
    bChiral = True
    d = 0.05
    timeout = None
    if(cmdl['-H2H']==True):
        bH2H = True
    if(cmdl['-H2Hpolar']==True):
        bH2Hpolar = True
    bH2heavy = False
    if(cmdl['-H2heavy']==True):
        bH2heavy = True
        bH2H = True
    if(cmdl.opt['-d'].is_set):
        d = cmdl['-d']
    if(cmdl.opt['-timeout'].is_set):
        timeout = cmdl['-timeout']
    if(cmdl['-chirality']==False):
        bChiral = False

#####################################
#### identify the methods to use ####
    bAlignment = True
    bMCS = True
    if cmdl.opt['-alignment'].is_set:
        bAlignment = cmdl['-alignment']
        if cmdl['-alignment']==True:
            if cmdl.opt['-mcs'].is_set and cmdl['-mcs']==False:
                bMCS = False
            if cmdl.opt['-mcs'].is_set==False:
                bMCS = False
    if cmdl.opt['-mcs'].is_set:
        bMCS = cmdl['-mcs']
        if cmdl['-mcs']==True:
            if cmdl.opt['-alignment'].is_set and cmdl['-alignment']==False:
                bAlignment = False
            if cmdl.opt['-alignment'].is_set==False:
                bAlignment = False
    if bMCS==False and bAlignment==False:
        print("No method (alignment, mcs) was selected.")
        sys.exit(0)
    print("Morphable atoms will be identified using the following methods:")
    print("Alignment: ",bAlignment)
    print("MCS: ",bMCS)
    print("\n")
######################################

    # read index
    read_from_idx1 = False
    read_from_idx2 = False
    if cmdl.opt['-n1'].is_set:
        read_from_idx1 = True
        idx1 = IndexFile(cmdl['-n1']).dic['scaffold']
    if cmdl.opt['-n2'].is_set:
        read_from_idx2 = True
        idx2 = IndexFile(cmdl['-n2']).dic['scaffold']

#args=parse_common_args(sys.argv,files,options, desc)

    # reformat PDB to read two letter atoms properly (why does it have to be so painful?)
    #randNum = randint(0, 9999)
    pid = os.getpid()
    pdbName1,atomNameID1,sigmaHoleID1 = reformatPDB(cmdl['-i1'],1,pid)
    pdbName2,atomNameID2,sigmaHoleID2 = reformatPDB(cmdl['-i2'],2,pid)

    mol1 = Chem.MolFromPDBFile(pdbName1,removeHs=False,sanitize=False)
    mol2 = Chem.MolFromPDBFile(pdbName2,removeHs=False,sanitize=False)

#    Chem.SanitizeMol(mol1)
#    Chem.SanitizeMol(mol2)

    os.remove(pdbName1)
    os.remove(pdbName2)

    molForMcs1 = cp.deepcopy(mol1)
    molForMcs2 = cp.deepcopy(mol2)
    try:
        rdmolops.AssignAtomChiralTagsFromStructure(mol1)
        rdmolops.AssignAtomChiralTagsFromStructure(mol2)
        molForMcs1 = cp.deepcopy(mol1)
        molForMcs2 = cp.deepcopy(mol2)
        rdmolops.AssignStereochemistry(mol1)
        rdmolops.AssignStereochemistry(mol2)
    except:
        print("Chirality not assigned")

#    mol1 = Chem.SDMolSupplier(cmdl['-i1'],removeHs=False,sanitize=True)
#    mol2 = Chem.SDMolSupplier(cmdl['-i2'],removeHs=False,sanitize=True)
#    mol1 = Chem.SDMolSupplier(cmdl['-i1'])
#    mol2 = Chem.SDMolSupplier(cmdl['-i2'])
#    mol1 = mol1[0]
#    mol2 = mol2[0]

    # deal with the -RingsOnly flag
    bRingsOnly = False
    bYesRings = checkRingsOnlyFlag(mol1,mol2)
    if(cmdl['-RingsOnly']==True):
        if( bYesRings==True ):
            bRingsOnly=True
        else:
            print("-RingsOnly flag is unset, because one (or both) molecule has no rings\n")

    n1 = []
    n2 = []

#########################
## make molecule copies #
    molcp1 = cp.deepcopy(mol1)
    molcp2 = cp.deepcopy(mol2)


#*******************************************************************************#
#*******************************************************************************#
#*******************************************************************************#
###########################################
########### mcs ###########################
    bMCSfailed = False
    n1mcs = []
    n2mcs = []
    if(bMCS==True):
        print("The topology matching approach will be used (MCS)")
        print("fmcs module: Copyright (c) 2012 Andrew Dalke Scientific AB\n")
        if(bRingsOnly==True):
            n1mcs,n2mcs = mcs(mol1,mol2,molForMcs1,molForMcs2,bH2H,bH2Hpolar,bH2heavy,cmdl['-dMCS'],bRingsOnly,d,bChiral,sigmaHoleID1,sigmaHoleID2,timeout)
        else:
            print("Trying to run an MCS using all atoms...")
            n1A,n2A = mcs(mol1,mol2,molForMcs1,molForMcs2,bH2H,bH2Hpolar,bH2heavy,cmdl['-dMCS'],False,d,bChiral,sigmaHoleID1,sigmaHoleID2,timeout)
            n1mcs = n1A
            n2mcs = n2A
            print("Size of mapping: ",len(n1mcs))
            if( bYesRings==True ):
                print("Trying to run an MCS using rings only...")
                n1B,n2B = mcs(mol1,mol2,molForMcs1,molForMcs2,bH2H,bH2Hpolar,bH2heavy,cmdl['-dMCS'],True,d,bChiral,sigmaHoleID1,sigmaHoleID2,timeout)
                if(len(n1A)<=len(n1B)):
                    print("Using ring only MCS result.")
                    n1mcs = n1B
                    n2mcs = n2B
                else:
                    print("Using all atom MCS result.")
        if len(n1mcs)==0:
            bMCSfailed = True
########################################
####### alignment ######################
    n1align = []
    n2align = []
    if( (bAlignment==True) or (bMCSfailed==True) ):
        if bMCSfailed==True:
            print("The MCS approach did not find any mapping")
        print("\nThe alignment approach will be used")
        print("Tosco, P., Balle, T. & Shiri, F. Open3DALIGN: an open-source software aimed at unsupervised ligand alignment. J Comput Aided Mol Des 25:777-83 (2011)")
        print("Alignment is based on atom logP contributions: S. A. Wildman and G. M. Crippen JCICS _39_ 868-873 (1999)\n")
        # only use rings
        if(bRingsOnly==True):
            n1align,n2align,pyO3A = o3a_alignment(mol1,mol2,bH2H,bH2Hpolar,bH2heavy,bRingsOnly,sigmaHoleID1,sigmaHoleID2,True,d)
        # else try both options and choose better
        else:
            print("Trying to align all atoms...")
            n1align,n2align,pyO3A = o3a_alignment(mol1,mol2,bH2H,bH2Hpolar,bH2heavy,False,sigmaHoleID1,sigmaHoleID2,True,d)
            print("Size of mapping: ",len(n1align))
            if( bYesRings==True ):
                print("Trying to align rings only...")
                mol1 = cp.deepcopy(molcp1)
                mol2 = cp.deepcopy(molcp2)
                n1B,n2B,pyO3A = o3a_alignment(mol1,mol2,bH2H,bH2Hpolar,bH2heavy,True,sigmaHoleID1,sigmaHoleID2,True,d)
                print("Size of mapping: ",len(n1B))
                if(len(n1align)<=len(n1B)):
                    print("Using ring only alignment result.")
                    n1align = n1B
                    n2align = n2B
                else:
                    print("Using all atom alignment result.")
    # select the better result: mcs or align
    if len(n1align)>=len(n1mcs):
        print("The final result is based on the O3A alignment.")
        n1 = n1align
        n2 = n2align
    else:
        print("The final result is based on the MCS.")
        n1 = n1mcs
        n2 = n2mcs
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
# also try the same procedure by inverting the ordering of mol1 and mol2 #
###########################################
########### mcs ###########################
    print("\n********************************************************************************************************************************")
    print("To ensure that the mapping is symmetric, i.e. mol1->mol2=mol2->mol1, we repeat the same procedure by swapping the molecule order.\n")
    print("********************************************************************************************************************************")
    bMCSfailed = False
    n1mcs = []
    n2mcs = []
    if(bMCS==True):
        print("Running MCS:")
        if(bRingsOnly==True):
            n2mcs,n1mcs = mcs(mol2,mol1,molForMcs2,molForMcs1,bH2H,bH2Hpolar,bH2heavy,cmdl['-dMCS'],bRingsOnly,d,bChiral,sigmaHoleID2,sigmaHoleID1,timeout)
        else:
            n2A,n1A = mcs(mol2,mol1,molForMcs2,molForMcs1,bH2H,bH2Hpolar,bH2heavy,cmdl['-dMCS'],False,d,bChiral,sigmaHoleID2,sigmaHoleID1,timeout)
            n1mcs = n1A
            n2mcs = n2A
            if( bYesRings==True ):
                n2B,n1B = mcs(mol2,mol1,molForMcs2,molForMcs1,bH2H,bH2Hpolar,bH2heavy,cmdl['-dMCS'],True,d,bChiral,sigmaHoleID2,sigmaHoleID1,timeout)
                if(len(n1A)<=len(n1B)):
                    n1mcs = n1B
                    n2mcs = n2B
        if len(n1mcs)==0:
            bMCSfailed = True
########################################
####### alignment ######################
    n1align = []
    n2align = []
    if( (bAlignment==True) or (bMCSfailed==True) ):
        print("\nRunning alignment:")
        if(bRingsOnly==True):
            n2align,n1align,pyO3A = o3a_alignment(mol2,mol1,bH2H,bH2Hpolar,bH2heavy,bRingsOnly,sigmaHoleID2,sigmaHoleID1,True,d)
        else:
            n2align,n1align,pyO3A = o3a_alignment(mol2,mol1,bH2H,bH2Hpolar,bH2heavy,False,sigmaHoleID2,sigmaHoleID1,True,d)
            print("Size of mapping: ",len(n1align))
            if( bYesRings==True ):
                mol1 = cp.deepcopy(molcp1)
                mol2 = cp.deepcopy(molcp2)
                n2B,n1B,pyO3A = o3a_alignment(mol2,mol1,bH2H,bH2Hpolar,bH2heavy,True,sigmaHoleID2,sigmaHoleID1,True,d)
                print("Size of mapping: ",len(n1B))
                if(len(n1align)<=len(n1B)):
                    n1align = n1B
                    n2align = n2B
    # select the better result (mcs or align) and compare with the non-inverted mapping
    if (len(n1align)<=len(n1)) and (len(n1mcs)<=len(n1)):
        print("Swapping of the molecules did not yield a better atom mapping.")
    if (len(n1align)>=len(n1mcs)) and (len(n1align)>len(n1)):
        print("The final result is based on the O3A alignment after swapping the molecule order (mol2-mol1).")
        n1 = n1align
        n2 = n2align
    elif (len(n1align)<len(n1mcs)) and (len(n1mcs)>len(n1)):
        print("The final result is based on the MCS after swapping the molecule order (mol2-mol1).")
        n1 = n1mcs
        n2 = n2mcs
    print("\n")
#*******************************************************************************#
#*******************************************************************************#
#*******************************************************************************#

    # a check
    if( len(n1) != len(n2) ):
        print("Warning: something went wrong.")
        print("Number of the morphable atoms in the ligands does not match.\n")

    # calculate score
    score = calcScore(mol1,mol2,n1,n2,bH2H,bH2heavy)

    # print some output
    print("FINAL RESULTS")
    if( bH2H==True or bH2heavy==True ):
        print("Atoms considered in mol1: ",mol1.GetNumAtoms())
        print("Atoms considered in mol2: ",mol2.GetNumAtoms())
    else:
        print("Atoms considered in mol1: ",mol1.GetNumHeavyAtoms())
        print("Atoms considered in mol2: ",mol2.GetNumHeavyAtoms())
    print("Morphable atoms in both molecules: ",len(n1),len(n2))
    print("Dissimilarity (distance) score: %.4f\n" % score)
    if(cmdl['-score']):
        fp = open(cmdl['-score'],'w')
        fp.write("Score: %.4f\n" % score)
        fp.close()

    restoreAtomNames(mol1,atomNameID1)
    restoreAtomNames(mol2,atomNameID2)

    # output
    if cmdl.opt['-opdb1'].is_set:
#        Chem.rdMolAlign.AlignMol(mol1,molcp1,atomMap=zip(n1,n1))
        Chem.MolToPDBFile(molcp1,cmdl['-opdb1'])
    if cmdl.opt['-opdb2'].is_set:
#       if( cmdl['-mcs']==True ):
        try:
#            Chem.rdMolAlign.AlignMol(mol2,mol1,atomMap=zip(n2,n1))
            Chem.rdMolAlign.AlignMol(mol2,molcp1,atomMap=list(zip(n2,n1)))
        except:
            print("Cannot superimpose -opdb2 structure. Maybe no morphable atoms have been found\n")
        Chem.MolToPDBFile(mol2,cmdl['-opdb2'])
    if cmdl.opt['-opdbm1'].is_set:
        mol = subMolByIndex(molcp1,n1)
        Chem.MolToPDBFile(mol,cmdl['-opdbm1'])
    if cmdl.opt['-opdbm2'].is_set:
        mol = subMolByIndex(mol2,n2)
        Chem.MolToPDBFile(mol,cmdl['-opdbm2'])

    # write out pairs
    pairsFile = cmdl['-o']
    write_pairs(n1,n2,pairsFile)

main( sys.argv )
