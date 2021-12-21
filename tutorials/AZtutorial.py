import pmx
from pmx.utils import create_folder
from pmx import gmx, ligand_alchemy, jobscript
import sys
import os,shutil
import re
import subprocess
import glob
import random
import pandas as pd
import numpy as np

class AZtutorial:
    """Class contains parameters for setting up free energy calculations

    Parameters
    ----------
    ...

    Attributes
    ----------
    ....

    """

    def __init__(self, **kwargs):
        
        # set gmxlib path
        gmx.set_gmxlib()
        
        # the results are summarized in a pandas framework
        self.resultsAll = pd.DataFrame()
        self.resultsSummary = pd.DataFrame()
        
        # paths
        self.workPath = './'
        self.mdpPath = '{0}/mdp'.format(self.workPath)
        self.proteinPath = None
        self.ligandPath = None
        
        # information about inputs
        self.protein = {} # protein[path]=path,protein[str]=pdb,protein[itp]=[itps],protein[posre]=[posres],protein[mols]=[molnames]
        self.ligands = {} # ligands[ligname]=path
        self.edges = {} # edges[edge_lig1_lig2] = [lig1,lig2]
        
        # parameters for the general setup
        self.replicas = 3        
        self.simTypes = ['em','eq','transitions']
        self.states = ['stateA','stateB']
        self.thermCycleBranches = ['water','protein']
                
        # simulation setup
        self.ff = 'amber99sb-star-ildn-mut.ff'
        self.boxshape = 'dodecahedron'
        self.boxd = 1.5
        self.water = 'tip3p'
        self.conc = 0.15
        self.pname = 'NaJ'
        self.nname = 'ClJ'
        
        # job submission params
        self.JOBqueue = 'SGE' # could be SLURM
        self.JOBsimtime = 24 # hours
        self.JOBsimcpu = 8 # CPU default
        self.JOBbGPU = True
        self.JOBmodules = []
        self.JOBsource = []
        self.JOBexport = []
        self.JOBgmx = 'gmx mdrun'

        for key, val in kwargs.items():
            setattr(self,key,val)
            
    def prepareFreeEnergyDir( self ):
        
        # protein={}
        # protein[path] = [path], protein[str] = pdb, protein[itp] = [itp], protein[posre] = [posre]
        self.proteinPath = self._read_path( self.proteinPath )
        self._protein = self._read_protein()
        
        # read ligands
        self.ligandPath = self._read_path( self.ligandPath )
        self._read_ligands()
        
        # read edges (directly or from a file)
        self._read_edges()
        
        # read mdpPath
        self.mdpPath = self._read_path( self.mdpPath )
        
        # workpath
        self.workPath = self._read_path( self.workPath )
        create_folder( self.workPath )
        
        # create folder structure
        self._create_folder_structure( )
        
        # print summary
        self._print_summary( )
                        
        # print folder structure
        self._print_folder_structure( )    
        
        print('DONE')
        
        
    # _functions to quickly get a path at different levels, e.g wppath, edgepath... like in _create_folder_structure
    def _get_specific_path( self, edge=None, bHybridStrTop=False, wp=None, state=None, r=None, sim=None ):
        if edge==None:
            return(self.workPath)       
        edgepath = '{0}/{1}'.format(self.workPath,edge)
        
        if bHybridStrTop==True:
            hybridStrPath = '{0}/hybridStrTop'.format(edgepath)
            return(hybridStrPath)

        if wp==None:
            return(edgepath)
        wppath = '{0}/{1}'.format(edgepath,wp)
        
        if state==None:
            return(wppath)
        statepath = '{0}/{1}'.format(wppath,state)
        
        if r==None:
            return(statepath)
        runpath = '{0}/run{1}'.format(statepath,r)
        
        if sim==None:
            return(runpath)
        simpath = '{0}/{1}'.format(runpath,sim)
        return(simpath)
                
    def _read_path( self, path ):
        return(os.path.abspath(path))
        
    def _read_ligands( self ):
        # read ligand folders
        ligs = glob.glob('{0}/*'.format(self.ligandPath))
        # get ligand names
        for l in ligs:
            lname = l.split('/')[-1]
            lnameTrunc = lname
            if lname.startswith('lig_'):
                lnameTrunc = lname[4:]
            elif lname.startswith('lig'):
                lnameTrunc = lname[3:]
            lpath = '{0}/{1}'.format(self.ligandPath,lname)
            self.ligands[lnameTrunc] = os.path.abspath(lpath)
 
    def _read_protein( self ):
        # read protein folder
        self.protein['path'] = os.path.abspath(self.proteinPath)
        # get folder contents
        self.protein['posre'] = []
        self.protein['itp'] = []
        self.protein['mols'] = [] # mols to add to .top
        self.protein['str'] = ''
        for l in glob.glob('{0}/*'.format(self.proteinPath)):
            fname = l.split('/')[-1]
            if '.itp' in fname: # posre or top
                if 'posre' in fname:
                    self.protein['posre'].append(os.path.abspath(l))
                else:
                    self.protein['itp'].append(os.path.abspath(l))
                    if fname.startswith('topol_'):
                        self.protein['mols'].append(fname[6:-4])
                    else:
                        self.protein['mols'].append(fname[:-4])                        
            if '.pdb' in fname:
                self.protein['str'] = fname
        self.protein['mols'].sort()
                
    def _read_edges( self ):
        # read from file
        try:
            if os.path.isfile( self.edges ):
                self._read_edges_from_file( self )
        # edge provided as an array
        except: 
            foo = {}
            for e in self.edges:
                key = 'edge_{0}_{1}'.format(e[0],e[1])
                foo[key] = e
            self.edges = foo
            
    def _read_edges_from_file( self ):
        self.edges = 'Edges read from file'
        
        
    def _create_folder_structure( self, edges=None ):
        # edge
        if edges==None:
            edges = self.edges        
        for edge in edges:
            print(edge)            
            edgepath = '{0}/{1}'.format(self.workPath,edge)
            create_folder(edgepath)
            
            # folder for hybrid ligand structures
            hybridTopFolder = '{0}/hybridStrTop'.format(edgepath)
            create_folder(hybridTopFolder)
            
            # water/protein
            for wp in self.thermCycleBranches:
                wppath = '{0}/{1}'.format(edgepath,wp)
                create_folder(wppath)
                
                # stateA/stateB
                for state in self.states:
                    statepath = '{0}/{1}'.format(wppath,state)
                    create_folder(statepath)
                    
                    # run1/run2/run3
                    for r in range(1,self.replicas+1):
                        runpath = '{0}/run{1}'.format(statepath,r)
                        create_folder(runpath)
                        
                        # em/eq_posre/eq/transitions
                        for sim in self.simTypes:
                            simpath = '{0}/{1}'.format(runpath,sim)
                            create_folder(simpath)
                            
    def _print_summary( self ):
        print('\n---------------------\nSummary of the setup:\n---------------------\n')
        print('   workpath: {0}'.format(self.workPath))
        print('   mdp path: {0}'.format(self.mdpPath))
        print('   protein files: {0}'.format(self.proteinPath))
        print('   ligand files: {0}'.format(self.ligandPath))
        print('   number of replicase: {0}'.format(self.replicas))        
        print('   edges:')
        for e in self.edges.keys():
            print('        {0}'.format(e))    
            
    def _print_folder_structure( self ):
        print('\n---------------------\nDirectory structure:\n---------------------\n')
        print('{0}/'.format(self.workPath))
        print('|')
        print('|--edge_X_Y')
        print('|--|--water')
        print('|--|--|--stateA')
        print('|--|--|--|--run1/2/3')
        print('|--|--|--|--|--em/eq_posre/eq/transitions')        
        print('|--|--|--stateB')
        print('|--|--|--|--run1/2/3')
        print('|--|--|--|--|--em/eq_posre/eq/transitions')       
        print('|--|--protein')
        print('|--|--|--stateA')
        print('|--|--|--|--run1/2/3')
        print('|--|--|--|--|--em/eq_posre/eq/transitions')        
        print('|--|--|--stateB')
        print('|--|--|--|--run1/2/3')
        print('|--|--|--|--|--em/eq_posre/eq/transitions')    
        print('|--|--hybridStrTop')        
        print('|--edge_..')
        
    def _be_verbose( self, process, bVerbose=False ):
        out = process.communicate()            
        if bVerbose==True:
            printout = out[0].splitlines()
            for o in printout:
                print(o)
        # error is printed every time                  
        printerr = out[1].splitlines()                
        for e in printerr:
            print(e)              
        
    def atom_mapping( self, edges=None, bVerbose=False ):
        print('-----------------------')
        print('Performing atom mapping')
        print('-----------------------')
        
        if edges==None:
            edges = self.edges        
        for edge in edges:
            print(edge)
            lig1 = self.edges[edge][0]
            lig2 = self.edges[edge][1]
            lig1path = '{0}/lig_{1}'.format(self.ligandPath,lig1)
            lig2path = '{0}/lig_{1}'.format(self.ligandPath,lig2)
            outpath = self._get_specific_path(edge=edge,bHybridStrTop=True)
            
            # params
            i1 = '{0}/mol_gmx.pdb'.format(lig1path)
            i2 = '{0}/mol_gmx.pdb'.format(lig2path)
            o1 = '{0}/pairs1.dat'.format(outpath)
            o2 = '{0}/pairs2.dat'.format(outpath)            
            opdb1 = '{0}/out_pdb1.pdb'.format(outpath)
            opdb2 = '{0}/out_pdb2.pdb'.format(outpath)
            opdbm1 = '{0}/out_pdbm1.pdb'.format(outpath)
            opdbm2 = '{0}/out_pdbm2.pdb'.format(outpath)
            score = '{0}/score.dat'.format(outpath)
            log = '{0}/mapping.log'.format(outpath)
            
            process = subprocess.Popen(['pmx','atomMapping',
                                '-i1',i1,
                                '-i2',i2,
                                '-o1',o1,
                                '-o2',o2,
                                '-opdb1',opdb1,
                                '-opdb2',opdb2,                                        
                                '-opdbm1',opdbm1,
                                '-opdbm2',opdbm2,
                                '-score',score,
                                '-log',log],
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)

            self._be_verbose( process, bVerbose=bVerbose )              
                
            process.wait()      
        print('DONE')            
            
            
    def hybrid_structure_topology( self, edges=None, bVerbose=False ):
        print('----------------------------------')
        print('Creating hybrid structure/topology')
        print('----------------------------------')

        if edges==None:
            edges = self.edges        
        for edge in edges:
            print(edge)
            lig1 = self.edges[edge][0]
            lig2 = self.edges[edge][1]
            lig1path = '{0}/lig_{1}'.format(self.ligandPath,lig1)
            lig2path = '{0}/lig_{1}'.format(self.ligandPath,lig2)
            outpath = self._get_specific_path(edge=edge,bHybridStrTop=True)
            
            # params
            i1 = '{0}/mol_gmx.pdb'.format(lig1path)
            i2 = '{0}/mol_gmx.pdb'.format(lig2path)
            itp1 = '{0}/MOL.itp'.format(lig1path)
            itp2 = '{0}/MOL.itp'.format(lig2path)            
            pairs = '{0}/pairs1.dat'.format(outpath)            
            oA = '{0}/mergedA.pdb'.format(outpath)
            oB = '{0}/mergedB.pdb'.format(outpath)
            oitp = '{0}/merged.itp'.format(outpath)
            offitp = '{0}/ffmerged.itp'.format(outpath)
            log = '{0}/hybrid.log'.format(outpath)
            
            process = subprocess.Popen(['pmx','ligandHybrid',
                                '-i1',i1,
                                '-i2',i2,
                                '-itp1',itp1,
                                '-itp2',itp2,
                                '-pairs',pairs,
                                '-oA',oA,  
                                '-oB',oB,
                                '-oitp',oitp,
                                '-offitp',offitp,
                                '-log',log],
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)

            self._be_verbose( process, bVerbose=bVerbose )                    
                
            process.wait()    
        print('DONE')
            
            
    def _make_clean_pdb(self, fnameIn,fnameOut,bAppend=False):
        # read 
        fp = open(fnameIn,'r')
        lines = fp.readlines()
        out = []
        for l in lines:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                out.append(l)
        fp.close()
        
        # write
        if bAppend==True:
            fp = open(fnameOut,'a')
        else:
            fp = open(fnameOut,'w')
        for l in out:
            fp.write(l)
        fp.close()
            
    def assemble_systems( self, edges=None ):
        print('----------------------')
        print('Assembling the systems')
        print('----------------------')

        if edges==None:
            edges = self.edges        
        for edge in edges:
            print(edge)            
            lig1 = self.edges[edge][0]
            lig2 = self.edges[edge][1]
            lig1path = '{0}/lig_{1}'.format(self.ligandPath,lig1)
            lig2path = '{0}/lig_{1}'.format(self.ligandPath,lig2)
            hybridStrTopPath = self._get_specific_path(edge=edge,bHybridStrTop=True)                    
            outLigPath = self._get_specific_path(edge=edge,wp='water')
            outProtPath = self._get_specific_path(edge=edge,wp='protein')
                        
            # Ligand structure
            self._make_clean_pdb('{0}/mergedA.pdb'.format(hybridStrTopPath),'{0}/init.pdb'.format(outLigPath))
            
            # Ligand+Protein structure
            self._make_clean_pdb('{0}/{1}'.format(self.proteinPath,self.protein['str']),'{0}/init.pdb'.format(outProtPath))
            self._make_clean_pdb('{0}/mergedA.pdb'.format(hybridStrTopPath),'{0}/init.pdb'.format(outProtPath),bAppend=True)
    
            # Ligand topology
            # ffitp
            ffitpOut = '{0}/ffmerged.itp'.format(hybridStrTopPath)
            ffitpIn1 = '{0}/ffMOL.itp'.format(lig1path)
            ffitpIn2 = '{0}/ffMOL.itp'.format(lig2path)
            ffitpIn3 = '{0}/ffmerged.itp'.format(hybridStrTopPath)
            pmx.ligand_alchemy._merge_FF_files( ffitpOut, ffsIn=[ffitpIn1,ffitpIn2,ffitpIn3] )        
            # top
            ligTopFname = '{0}/topol.top'.format(outLigPath)
            ligFFitp = '{0}/ffmerged.itp'.format(hybridStrTopPath)
            ligItp ='{0}/merged.itp'.format(hybridStrTopPath)
            itps = [ligFFitp,ligItp]
            systemName = 'ligand in water'
            self._create_top(fname=ligTopFname,itp=itps,systemName=systemName)
            
            # Ligand+Protein topology
            # top
            protTopFname = '{0}/topol.top'.format(outProtPath)
            mols = []
            for m in self.protein['mols']:
                mols.append([m,1])
            mols.append(['MOL',1])
            itps.extend(self.protein['itp'])
            systemName = 'protein and ligand in water'
            self._create_top(fname=protTopFname,itp=itps,mols=mols,systemName=systemName)            
        print('DONE')            
        
            
    def _create_top( self, fname='topol.top',  
                   itp=['merged.itp'], mols=[['MOL',1]],
                   systemName='simulation system',
                   destination='',toppaths=[] ):

        fp = open(fname,'w')
        # ff itp
        fp.write('#include "%s/forcefield.itp"\n' % self.ff)
        # additional itp
        for i in itp:
            fp.write('#include "%s"\n' % i) 
        # water itp
        fp.write('#include "%s/%s.itp"\n' % (self.ff,self.water)) 
        # ions
        fp.write('#include "%s/ions.itp"\n\n' % self.ff)
        # system
        fp.write('[ system ]\n')
        fp.write('{0}\n\n'.format(systemName))
        # molecules
        fp.write('[ molecules ]\n')
        for mol in mols:
            fp.write('%s %s\n' %(mol[0],mol[1]))
        fp.close()

        
    def _clean_backup_files( self, path ):
        toclean = glob.glob('{0}/*#'.format(path)) 
        for clean in toclean:
            os.remove(clean)        
    
    def boxWaterIons( self, edges=None, bBoxLig=True, bBoxProt=True, 
                                        bWatLig=True, bWatProt=True,
                                        bIonLig=True, bIonProt=True):
        print('----------------')
        print('Box, water, ions')
        print('----------------')
        
        if edges==None:
            edges = self.edges
        for edge in edges:
            print(edge)            
            outLigPath = self._get_specific_path(edge=edge,wp='water')
            outProtPath = self._get_specific_path(edge=edge,wp='protein')
            
            # box ligand
            if bBoxLig==True:
                inStr = '{0}/init.pdb'.format(outLigPath)
                outStr = '{0}/box.pdb'.format(outLigPath)
                gmx.editconf(inStr, o=outStr, bt=self.boxshape, d=self.boxd, other_flags='')                
            # box protein
            if bBoxProt==True:            
                inStr = '{0}/init.pdb'.format(outProtPath)
                outStr = '{0}/box.pdb'.format(outProtPath)
                gmx.editconf(inStr, o=outStr, bt=self.boxshape, d=self.boxd, other_flags='')
                
            # water ligand
            if bWatLig==True:            
                inStr = '{0}/box.pdb'.format(outLigPath)
                outStr = '{0}/water.pdb'.format(outLigPath)
                top = '{0}/topol.top'.format(outLigPath)
                gmx.solvate(inStr, cs='spc216.gro', p=top, o=outStr)
            # water protein
            if bWatProt==True:            
                inStr = '{0}/box.pdb'.format(outProtPath)
                outStr = '{0}/water.pdb'.format(outProtPath)
                top = '{0}/topol.top'.format(outProtPath)
                gmx.solvate(inStr, cs='spc216.gro', p=top, o=outStr)                
            
            # ions ligand
            if bIonLig:
                inStr = '{0}/water.pdb'.format(outLigPath)
                outStr = '{0}/ions.pdb'.format(outLigPath)
                mdp = '{0}/em_l0.mdp'.format(self.mdpPath)
                tpr = '{0}/tpr.tpr'.format(outLigPath)
                top = '{0}/topol.top'.format(outLigPath)
                mdout = '{0}/mdout.mdp'.format(outLigPath)
                gmx.grompp(f=mdp, c=inStr, p=top, o=tpr, maxwarn=4, other_flags=' -po {0}'.format(mdout))        
                gmx.genion(s=tpr, p=top, o=outStr, conc=self.conc, neutral=True, 
                      other_flags=' -pname {0} -nname {1}'.format(self.pname, self.nname))
            # ions protein
            if bIonProt:
                inStr = '{0}/water.pdb'.format(outProtPath)
                outStr = '{0}/ions.pdb'.format(outProtPath)
                mdp = '{0}/em_l0.mdp'.format(self.mdpPath)
                tpr = '{0}/tpr.tpr'.format(outProtPath)
                top = '{0}/topol.top'.format(outProtPath)
                mdout = '{0}/mdout.mdp'.format(outProtPath)
                gmx.grompp(f=mdp, c=inStr, p=top, o=tpr, maxwarn=4, other_flags=' -po {0}'.format(mdout))        
                gmx.genion(s=tpr, p=top, o=outStr, conc=self.conc, neutral=True, 
                      other_flags=' -pname {0} -nname {1}'.format(self.pname, self.nname))                
                     
            # clean backed files
            self._clean_backup_files( outLigPath )
            self._clean_backup_files( outProtPath )
        print('DONE')
            
    def _prepare_single_tpr( self, simpath, toppath, state, simType, empath=None, frameNum=0 ):
        
        mdpPrefix = ''
        if simType=='em':
            mdpPrefix = 'em'
        elif simType=='eq':
            mdpPrefix = 'eq'
        elif simType=='transitions':
            mdpPrefix = 'ti'        
            
        top = '{0}/topol.top'.format(toppath)
        tpr = '{0}/tpr.tpr'.format(simpath)
        mdout = '{0}/mdout.mdp'.format(simpath)
        # mdp
        if state=='stateA':
            mdp = '{0}/{1}_l0.mdp'.format(self.mdpPath,mdpPrefix)
        else:
            mdp = '{0}/{1}_l1.mdp'.format(self.mdpPath,mdpPrefix)
        # str
        if simType=='em':
            inStr = '{0}/ions.pdb'.format(toppath)
        elif simType=='eq':
            inStr = '{0}/confout.gro'.format(empath)
        elif simType=='transitions':
            inStr = '{0}/frame{1}.gro'.format(simpath,frameNum)
            tpr = '{0}/ti{1}.tpr'.format(simpath,frameNum)
            
        gmx.grompp(f=mdp, c=inStr, p=top, o=tpr, maxwarn=4, other_flags=' -po {0}'.format(mdout))
        self._clean_backup_files( simpath )
                    
         
    def prepare_simulation( self, edges=None, simType='em', bLig=True, bProt=True ):
        print('-----------------------------------------')
        print('Preparing simulation: {0}'.format(simType))
        print('-----------------------------------------')
        
        mdpPrefix = ''
        if simType=='em':
            mdpPrefix = 'em'
        elif simType=='eq':
            mdpPrefix = 'eq'
        elif simType=='transitions':
            mdpPrefix = 'ti'
        
        if edges==None:
            edges = self.edges
        for edge in edges:
            print(edge)
            ligTopPath = self._get_specific_path(edge=edge,wp='water')
            protTopPath = self._get_specific_path(edge=edge,wp='protein')            
            
            for state in self.states:
                for r in range(1,self.replicas+1):
                    
                    # ligand
                    if bLig==True:
                        wp = 'water'
                        simpath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim=simType)
                        empath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim='em')
                        toppath = ligTopPath
                        self._prepare_single_tpr( simpath, toppath, state, simType, empath )
                    
                    # protein
                    if bProt==True:
                        wp = 'protein'
                        simpath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim=simType)
                        empath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim='em')
                        toppath = protTopPath
                        self._prepare_single_tpr( simpath, toppath, state, simType, empath )    
        print('DONE')
        

    def _run_mdrun( self, tpr=None, ener=None, confout=None, mdlog=None, 
                    cpo=None, trr=None, xtc=None, dhdl=None, bVerbose=False):
        # EM
        if xtc==None:
            process = subprocess.Popen(['gmx','mdrun',
                                '-s',tpr,
                                '-e',ener,
                                '-c',confout,
                                '-o',trr,                                        
                                '-g',mdlog],
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
            self._be_verbose( process, bVerbose=bVerbose )                    
            process.wait()           
        # other FE runs
        else:
            process = subprocess.Popen(['gmx','mdrun',
                                '-s',tpr,
                                '-e',ener,
                                '-c',confout,
                                '-dhdl',dhdl,
                                '-x',xtc,
                                '-o',trr,
                                '-cpo',cpo,                                        
                                '-g',mdlog],
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
            self._be_verbose( process, bVerbose=bVerbose )                    
            process.wait()           
            

    def run_simulation_locally( self, edges=None, simType='em', bLig=True, bProt=True, bVerbose=False ):
        print('-------------------------------------------')
        print('Run simulation locally: {0}'.format(simType))
        print('-------------------------------------------')
        
        if edges==None:
            edges = self.edges
        for edge in edges:
            
            for state in self.states:
                for r in range(1,self.replicas+1):            
                    
                    # ligand
                    if bLig==True:
                        wp = 'water'
                        print('Running: LIG {0} {1} run{2}'.format(edge,state,r))
                        simpath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim=simType)
                        tpr = '{0}/tpr.tpr'.format(simpath)
                        ener = '{0}/ener.edr'.format(simpath)
                        confout = '{0}/confout.gro'.format(simpath)
                        mdlog = '{0}/md.log'.format(simpath)
                        trr = '{0}/traj.trr'.format(simpath)                        
                        self._run_mdrun(tpr=tpr,trr=trr,ener=ener,confout=confout,mdlog=mdlog,bVerbose=bVerbose)
                        self._clean_backup_files( simpath )
                    
                    # protein
                    if bProt==True:
                        wp = 'protein'
                        print('Running: PROT {0} {1} run{2}'.format(edge,state,r))
                        simpath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim=simType)
                        tpr = '{0}/tpr.tpr'.format(simpath)
                        ener = '{0}/ener.edr'.format(simpath)
                        confout = '{0}/confout.gro'.format(simpath)
                        mdlog = '{0}/md.log'.format(simpath)
                        trr = '{0}/traj.trr'.format(simpath)                                                
                        self._run_mdrun(tpr=tpr,trr=trr,ener=ener,confout=confout,mdlog=mdlog,bVerbose=bVerbose)
                        self._clean_backup_files( simpath )
        print('DONE')
 
    def prepare_jobscripts( self, edges=None, simType='em', bLig=True, bProt=True ):
        print('---------------------------------------------')
        print('Preparing jobscripts for: {0}'.format(simType))
        print('---------------------------------------------')
        
        jobfolder = '{0}/{1}_jobscripts'.format(self.workPath,simType)
        os.system('mkdir {0}'.format(jobfolder))
        
        if edges==None:
            edges = self.edges
            
        counter = 0
        for edge in edges:
            
            for state in self.states:
                for r in range(1,self.replicas+1):            
                    
                    # ligand
                    if bLig==True:
                        wp = 'water'
                        simpath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim=simType)
                        jobfile = '{0}/jobscript{1}'.format(jobfolder,counter)
                        jobname = 'lig_{0}_{1}_{2}_{3}'.format(edge,state,r,simType)
                        job = pmx.jobscript.Jobscript(fname=jobfile,
                                        queue=self.JOBqueue,simcpu=self.JOBsimcpu,
                                        jobname=jobname,modules=self.JOBmodules,source=self.JOBsource,
                                        gmx=self.JOBgmx)
                        cmd1 = 'cd {0}'.format(simpath)
                        cmd2 = '$GMXRUN -s tpr.tpr'
                        job.cmds = [cmd1,cmd2]                        
                        if simType=='transitions':
                            self._commands_for_transitions( simpath, job )
                        job.create_jobscript()
                        counter+=1
                    
                    # protein
                    if bProt==True:
                        wp = 'protein'
                        simpath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim=simType)
                        jobfile = '{0}/jobscript{1}'.format(jobfolder,counter)
                        jobname = 'prot_{0}_{1}_{2}_{3}'.format(edge,state,r,simType)
                        job = pmx.jobscript.Jobscript(fname=jobfile,
                                        queue=self.JOBqueue,simcpu=self.JOBsimcpu,
                                        jobname=jobname,modules=self.JOBmodules,source=self.JOBsource,
                                        gmx=self.JOBgmx)
                        cmd1 = 'cd {0}'.format(simpath)
                        cmd2 = '$GMXRUN -s tpr.tpr'
                        job.cmds = [cmd1,cmd2]
                        if simType=='transitions':
                            self._commands_for_transitions( simpath, job )                        
                        job.create_jobscript()
                        counter+=1
                        
        #######
        self._submission_script( jobfolder, counter, simType )
        print('DONE')
        
    def _commands_for_transitions( self, simpath, job ):
        for i in range(1,81):
            if self.JOBqueue=='SGE':
                cmd1 = 'cd $TMPDIR'
                cmd2 = 'cp {0}/ti$SGE_TASK_ID.tpr tpr.tpr'.format(simpath)
                cmd3 = '$GMXRUN -s tpr.tpr -dhdl dhdl$SGE_TASK_ID.xvg'.format(simpath)
                cmd4 = 'cp dhdl$SGE_TASK_ID.xvg {0}/.'.format(simpath)
                job.cmds = [cmd1,cmd2,cmd3,cmd4]
            elif self.JOBqueue=='SLURM':
                cmd1 = 'cd {0}'.format(simpath)
                cmd2 = 'for i in {1..81};do'
                cmd3 = '$GMXRUN -s ti$i.tpr -dhdl dhdl$i'
                cmd4 = 'done'
                job.cmds = [cmd1,cmd2,cmd3,cmd4]
        
        
    def _submission_script( self, jobfolder, counter, simType='eq' ):
        fname = '{0}/submit.py'.format(jobfolder)
        fp = open(fname,'w')
        fp.write('import os\n')
        fp.write('for i in range(0,{0}):\n'.format(counter))
        if self.JOBqueue=='SGE':
            cmd = '\'qsub jobscript{0}\'.format(i)'
            if simType=='transitions':
                cmd = '\'qsub -t 1-80:1 jobscript{0}\'.format(i)'
        elif self.JOBqueue=='SLURM':
            cmd = '\'sbatch jobscript{0}\'.format(i)'
        fp.write('    os.system({0})\n'.format(cmd))
        fp.close()

    def _extract_snapshots( self, eqpath, tipath ):
        tpr = '{0}/tpr.tpr'.format(eqpath)
        trr = '{0}/traj.trr'.format(eqpath)
        frame = '{0}/frame.gro'.format(tipath)
        
        gmx.trjconv(s=tpr,f=trr,o=frame, sep=True, ur='compact', pbc='mol', other_flags=' -b 2250')
        # move frame0.gro to frame80.gro
        cmd = 'mv {0}/frame0.gro {0}/frame80.gro'.format(tipath)
        os.system(cmd)
        
        self._clean_backup_files( tipath )
        
        
    def prepare_transitions( self, edges=None, bLig=True, bProt=True, bGenTpr=True ):
        print('---------------------')
        print('Preparing transitions')
        print('---------------------')
        
        if edges==None:
            edges = self.edges
        for edge in edges:
            ligTopPath = self._get_specific_path(edge=edge,wp='water')
            protTopPath = self._get_specific_path(edge=edge,wp='protein')            
            
            for state in self.states:
                for r in range(1,self.replicas+1):
                    
                    # ligand
                    if bLig==True:
                        print('Preparing: LIG {0} {1} run{2}'.format(edge,state,r))
                        wp = 'water'
                        eqpath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim='eq')
                        tipath = simpath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim='transitions')
                        toppath = ligTopPath
                        self._extract_snapshots( eqpath, tipath )
                        if bGenTpr==True:
                            for i in range(1,81):
                                self._prepare_single_tpr( tipath, toppath, state, simType='transitions',frameNum=i )
                    
                    # protein
                    if bProt==True:
                        print('Preparing: PROT {0} {1} run{2}'.format(edge,state,r))
                        wp = 'protein'
                        eqpath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim='eq')
                        tipath = simpath = self._get_specific_path(edge=edge,wp=wp,state=state,r=r,sim='transitions')                        
                        toppath = protTopPath
                        self._extract_snapshots( eqpath, tipath )
                        if bGenTpr==True:
                            for i in range(1,81):
                                self._prepare_single_tpr( tipath, toppath, state, simType='transitions',frameNum=i )
        print('DONE')  
        
        
    def _run_analysis_script( self, analysispath, stateApath, stateBpath, bVerbose=False ):
        fA = ' '.join( glob.glob('{0}/*xvg'.format(stateApath)) )
        fB = ' '.join( glob.glob('{0}/*xvg'.format(stateBpath)) )
        oA = '{0}/integ0.dat'.format(analysispath)
        oB = '{0}/integ1.dat'.format(analysispath)
        wplot = '{0}/wplot.png'.format(analysispath)
        o = '{0}/results.txt'.format(analysispath)

        cmd = 'pmx analyse -fA {0} -fB {1} -o {2} -oA {3} -oB {4} -w {5} -t {6} -b {7}'.format(\
                                                                            fA,fB,o,oA,oB,wplot,298,100) 
        os.system(cmd)
        
        if bVerbose==True:
            fp = open(o,'r')
            lines = fp.readlines()
            fp.close()
            bPrint = False
            for l in lines:
                if 'ANALYSIS' in l:
                    bPrint=True
                if bPrint==True:
                    print(l,end='')
        
    def run_analysis( self, edges=None, bLig=True, bProt=True, bParseOnly=False, bVerbose=False ):
        print('----------------')
        print('Running analysis')
        print('----------------')
        
        if edges==None:
            edges = self.edges
        for edge in edges:
            print(edge)
            
            for r in range(1,self.replicas+1):
                
                # ligand
                if bLig==True:
                    wp = 'water'
                    analysispath = '{0}/analyse{1}'.format(self._get_specific_path(edge=edge,wp=wp),r)
                    create_folder(analysispath)
                    stateApath = self._get_specific_path(edge=edge,wp=wp,state='stateA',r=r,sim='transitions')
                    stateBpath = self._get_specific_path(edge=edge,wp=wp,state='stateB',r=r,sim='transitions')
                    self._run_analysis_script( analysispath, stateApath, stateBpath, bVerbose=bVerbose )
                    
                # protein
                if bProt==True:
                    wp = 'protein'
                    analysispath = '{0}/analyse{1}'.format(self._get_specific_path(edge=edge,wp=wp),r)
                    create_folder(analysispath)
                    stateApath = self._get_specific_path(edge=edge,wp=wp,state='stateA',r=r,sim='transitions')
                    stateBpath = self._get_specific_path(edge=edge,wp=wp,state='stateB',r=r,sim='transitions')
                    self._run_analysis_script( analysispath, stateApath, stateBpath, bVerbose=bVerbose )
        print('DONE')
        
        
    def _read_neq_results( self, fname ):
        fp = open(fname,'r')
        lines = fp.readlines()
        fp.close()
        out = []
        for l in lines:
            l = l.rstrip()
            foo = l.split()
            if 'BAR: dG' in l:
                out.append(float(foo[-2]))
            elif 'BAR: Std Err (bootstrap)' in l:
                out.append(float(foo[-2]))
            elif 'BAR: Std Err (analytical)' in l:
                out.append(float(foo[-2]))      
            elif '0->1' in l:
                out.append(int(foo[-1]))      
            elif '1->0' in l:
                out.append(int(foo[-1]))
        return(out)         
    
    def _fill_resultsAll( self, res, edge, wp, r ):
        rowName = '{0}_{1}_{2}'.format(edge,wp,r)
        self.resultsAll.loc[rowName,'val'] = res[2]
        self.resultsAll.loc[rowName,'err_analyt'] = res[3]
        self.resultsAll.loc[rowName,'err_boot'] = res[4]
        self.resultsAll.loc[rowName,'framesA'] = res[0]
        self.resultsAll.loc[rowName,'framesB'] = res[1]
        
    def _summarize_results( self, edges ):
        bootnum = 1000
        for edge in edges:
            for wp in ['water','protein']:
                dg = []
                erra = []
                errb = []
                distra = []
                distrb = []
                for r in range(1,self.replicas+1):
                    rowName = '{0}_{1}_{2}'.format(edge,wp,r)
                    dg.append( self.resultsAll.loc[rowName,'val'] )
                    erra.append( self.resultsAll.loc[rowName,'err_analyt'] )
                    errb.append( self.resultsAll.loc[rowName,'err_boot'] )
                    distra.append(np.random.normal(self.resultsAll.loc[rowName,'val'],self.resultsAll.loc[rowName,'err_analyt'] ,size=bootnum))
                    distrb.append(np.random.normal(self.resultsAll.loc[rowName,'val'],self.resultsAll.loc[rowName,'err_boot'] ,size=bootnum))
                  
                rowName = '{0}_{1}'.format(edge,wp)
                distra = np.array(distra).flatten()
                distrb = np.array(distrb).flatten()

                if self.replicas==1:
                    self.resultsAll.loc[rowName,'val'] = dg[0]                              
                    self.resultsAll.loc[rowName,'err_analyt'] = erra[0]
                    self.resultsAll.loc[rowName,'err_boot'] = errb[0]
                else:
                    self.resultsAll.loc[rowName,'val'] = np.mean(dg)
                    self.resultsAll.loc[rowName,'err_analyt'] = np.sqrt(np.var(distra)/float(self.replicas))
                    self.resultsAll.loc[rowName,'err_boot'] = np.sqrt(np.var(distrb)/float(self.replicas))
                    
            #### also collect resultsSummary
            rowNameWater = '{0}_{1}'.format(edge,'water')
            rowNameProtein = '{0}_{1}'.format(edge,'protein')            
            dg = self.resultsAll.loc[rowNameProtein,'val'] - self.resultsAll.loc[rowNameWater,'val']
            erra = np.sqrt( np.power(self.resultsAll.loc[rowNameProtein,'err_analyt'],2.0) \
                            - np.power(self.resultsAll.loc[rowNameWater,'err_analyt'],2.0) )
            errb = np.sqrt( np.power(self.resultsAll.loc[rowNameProtein,'err_boot'],2.0) \
                            - np.power(self.resultsAll.loc[rowNameWater,'err_boot'],2.0) )
            rowName = edge
            self.resultsSummary.loc[rowName,'val'] = dg
            self.resultsSummary.loc[rowName,'err_analyt'] = erra
            self.resultsSummary.loc[rowName,'err_boot'] = errb
            
                    
    def analysis_summary( self, edges=None ):
        if edges==None:
            edges = self.edges
            
        for edge in edges:
            for r in range(1,self.replicas+1):
                for wp in ['water','protein']:
                    analysispath = '{0}/analyse{1}'.format(self._get_specific_path(edge=edge,wp=wp),r)
                    resultsfile = '{0}/results.txt'.format(analysispath)
                    res = self._read_neq_results( resultsfile )
                    self._fill_resultsAll( res, edge, wp, r )
        
        # the values have been collected now
        # let's calculate ddGs
        self._summarize_results( edges )

