Note: the formatting of this document may look odd if you do not have "word wrap" switched on in your text editor. 

This repo is split into the following subfolders:


########################################### KEMP VARIANTS ###############################################

######### Classical_MD_Simulations: #########
The starting structures used to perform all MD simulations in our paper are collected at Starting_Structures folder: 
        * GNCA4-WT.pdb
        * GNCA4-12.pdb
        * V4.pdb
        * V4-4.pdb

Alongside the input files used for the equilibration and production of the MD simulations are at Inputs_MD folder. Simulations were performed with Amber and the final step used to prepare each structure was with tleap (part of amber tools).

The parameters for the substrate (5-nitrobenzioxazole) are at Parameter_SUB folder. 

Please see the paper for details on how each structure was prepared. 


######### EVB_Simulations ######### 
EVB simulations were performed using Q6. 
We have simulated four systems: 

	* GNCA4-WT
	* GNCA4-12
	* V4	
	* V4-4

The order of operations for the equil files is as follows:
"dyn_min_wat_noshake", "dyn_min_wat", "dyn_min_wat_eq", "dyn_cool", "dyn_min", "dyn_warm30k", "dyn_warm150k", "dyn_warm300k", "dyn_eq300k_{1..6}".

With 30 replicas performed per complex. 

The production runs were performed with the input called "qfep.inp"


########################################### PTE System ###############################################

######### Classical_MD_Simulations: #########
The starting structures used to perform all MD simulations in our paper are collected at Starting_Structures folder:
        * 4pcp_R0_clear_propka_chimhis_2NH_complex.pdb
        * 4pcp_R0_clear_propka_chimhis_complex.pdb

Alongside the input files used for the equilibration and production of the MD simulations are at Inputs_MD folder. Simulations were performed with Amber and the final step used to prepare each structure was with tleap (part of amber tools).

The parameters for the substrate (2-naphthyl-hexanoate, 2NH) and carboxylysine group (KCX) are at Additional_Parameter folder.

Please see the paper for details on how each structure and parameters were prepared. 

########################################### CHI System ###############################################

######### Classical_MD_Simulations: #########
The starting structures used to perform all MD simulations in our paper are collected at Starting_Structures folder:
        * 5WKR_APO.pdb
        * 5WKR_SUB.pdb

Alongside the input files used for the equilibration and production of the MD simulations are at Inputs_MD folder. Simulations were performed with Amber and the final step used to prepare each structure was with tleap (part of amber tools).

The parameters for the substrate (chalconaringenin) Parameter_SUB folder.

Please see the paper for details on how each structures were prepared. 


