
## Summary ##
This is a collection of the input files for 13 protein-ligand systems (482 ligand perturbations):
PDE2, CMET, Bace (Hunt et al), Bace p2, Galectin, JNK1, Tyk2, Bace, MCL1, CDK2, Thrombin, PTP1b, P38

## Citation: ##
  Large Scale Relative Protein Ligand Binding Affinities Using Non-Equilibrium Alchemy.  
  Gapsys, Pérez-Benito, Aldeghi, Seeliger, Van Vlijmen, Tresadern and de Groot.  
  2020. Chemical Science. 10.1039/C9SC03754C

## Folder structure: ##
   ### For every protein-ligand dataset ###
	- ligands_gaff2: for every ligand two topology files are present. 
	ffMOL.itp contains the atomtypes, the rest of topology parameters are in MOL.itp. 
	Structure is in the mol_gmx.pdb file. Force field: Gaff v2.1
	
	- ligands_cgenff: topology and structure for cgenff. 
	CGenFF v3.0.1 with the atom typing based on MATCH was used for most of the ligands, 
	except for Bace inhibitors, where paramchem.org web-server was used in combination with CGenFF v4.1.
	
	- protein_amber: structure and topology for protein and, 
	if available, co-crystallized waters and ions. Force field: amber99sb*ILDN
	
	- protein_charmm: structure and topology for Charmm36m force field
	
	- transformations_gaff2: edge information and hybrid structures/topologies for gaff.
	
	- transformations_cgenff: edge information for cgenff. 
	
  ddg_data: the folder contains calculated ddG values for all the protein-ligand datasets

  dg_data_allRepeats: gaff and cgenff force field dG data for transitions in water and protein separately for each repeat

  mdp: simulation parameter files
  ```
	- em_l0.mdp and em_l1.mdp: energy minimization for the states A and B
	- eq_nvt_l0.mdp and eq_nvt_l1.mdp: 10 ps NVT equilibration for the states A and B
	- eq_l0.mdp and eq_l1.mdp: 6 ns equilibrium run for the states A and B
	- ti_l0_gmx46.mdp and ti_l1_gmx46.mdp: 50 ps transition A->B and B->A. 
These files are meant to be used with a custom gmx4.6 version with an alternative soft-core formulation: 
pmx.mpibpc.mpg.de/gromacs462_newsc.tar
	- ti_l0_gmx2018.mdp and ti_l1_gmx2018.mdp: 50 ps transition A->B and B->A. 
Files to be used with gmx2018 and higher versions
```
