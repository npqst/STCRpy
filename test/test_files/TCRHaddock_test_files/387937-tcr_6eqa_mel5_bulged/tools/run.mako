! run.cns
!    The file containing all parameters for HADDOCK
!
! ***********************************************************************
! * Copyright 2003-2018 Alexandre Bonvin, Utrecht University.           *
! * Originally adapted from Aria 1.2 from Nilges and Linge, EMBL.       *
! * All rights reserved.                                                *
! * This code is part of the HADDOCK software and governed by its       *
! * license. Please see the LICENSE file that should have been included *
! * as part of this package.                                            *
! ***********************************************************************
!
module(
iteration;
filenames;
data;
iterations;
saprotocol;
refine;
toppar;
analysis;
)

{+ File: run.cns +}
{+ Description: this file contains all necessary information to run HADDOCK. +}

{+ Authors: Ezgi Karaca, Joao Rodrigues, Mikael Trellet, Alexandre Bonvin<br>
Adapted from HADOOCK version 2.4 <br><br>
Initially adapted from ARIA of Nilges and Linge +}

{+ Please cite the following references when using this protocol: +}
{+ reference: Cyril Dominguez, Rolf Boelens and Alexandre M.J.J. Bonvin (2003).  HADDOCK: a protein-protein docking approach
based on biochemical and/or biophysical information. <i>J. Am. Chem. Soc.</i> <b>125</b>, 1731-1737.
<p>
<b>When using <i>residual dipolar couplings</i> in HADDOCK cite in addition:</b><p>
<LI>A.D.J. van Dijk, D. Fushman and A.M.J.J. Bonvin (2005). Various strategies of using residual dipolar
couplings in NMR-driven protein docking: Application to Lys48-linked di-ubiquitin and validation against
15N-relaxation data. <EM>Proteins: Struc. Funct. & Bioinformatics</EM>, <STRONG>60</STRONG>, 367-381.</li>
<p>
<b>When using <i>diffusion anisotropy data</i> in HADDOCK cite in addition:</b><p>
<li>A.D.J. van Dijk, R. Kaptein, R. Boelens and A.M.J.J. Bonvin (2006). Combining NMR relaxation with
chemical shift perturbation data to drive protein-protein docking. <EM>J. Biomol. NMR</EM>,
<STRONG>34</STRONG>, 237-244.</li>
<p>
<b>When using <i>solvated docking</i> in HADDOCK cite in addition:</b><p>
<li>A.D.J. van Dijk and A.M.J.J. Bonvin (2006). Solvated docking: introducing water into the modelling
of biomolecular complexes. <EM>Bioinformatics</EM>,  <STRONG>22</STRONG> 2340-2347.
<p>
<b>When performing <i>flexible protein-DNA docking</i> using HADDOCK cite in addition:</b><p>
<li>M. van Dijk, A.D.J. van Dijk, V. Hsu, R. Boelens and  A.M.J.J. Bonvin (2006).
Information-driven Protein-DNA Docking using HADDOCK: it is a matter of flexibility.
<EM>Nucl. Acids Res.</EM>, <STRONG>34</STRONG> 3317-3325.</li>
<p>
<b>When performing the Nmolecule integrative modelling protocol please cite:</b><p>
<li>Ezgi Karaca, Joao P.G.L.M. Rodrigues, Andrea Graziadei, Alexandre M.J.J. Bonvin, Teresa Carlomagno (2017).
An Integrative Framework for Structure Determination of Molecular Machines.
<EM>Nature Methods</EM>, Advanced Online Publication.</li>
+}

{- Guidelines for using this file:
   - all strings must be quoted by double-quotes
   - logical variables (true/false) are not quoted
   - do not remove any evaluate statements from the file
   - pathnames should not exceed 80 characters -}
{- begin block parameter definition -} define(


{======== number of molecules for docking ==================}
{* number of components *}
<% ncomp = len(partners) %>
{===>} ncomponents=${ncomp};

{======================= filenames =========================}
{*  the name of your current project *}
{*  this will be used as name for the generated structures *}
{===>} fileroot="${fileroot}";

{* RUN directory *}
{*  the absolute path of your current run, e.g. /home/haddock/run1*}
{===>} run_dir="${run_dir}";
## Consider having a protein dictionary with the different values ({'A': {'pdb_file': 'pdb_name', 'psf_file': 'psf_name',
## 'root': 'protein1', dna: False })
<% nb_mol_max = 20 %>
% for n in range(1,ncomp+1):
{* PDB file of molecule ${n} *}
{===>} prot_coor_mol${n}="${partners[n]['pdb_file']}";
{* PSF file of molecule ${n} *}
{===>} prot_psf_mol${n}="${partners[n]['psf_file']}";
{* segid of molecule ${n} *}
{===>} prot_segid_mol${n}="${partners[n]['segid']}";
{* fileroot of molecule ${n} *}
{===>} prot_root_mol${n}="${partners[n]['root']}";
{* Fix Molecule at Origin during it0 *}
{+ choice: true false +}
{===>} fix_origin_mol${n}=${partners[n]['fix_origin']};
{* Is molecule ${n} nucleic acid? *}
{+ choice: true false +}
{===>} dna_mol${n}=${partners[n]['dna']};
{* Is molecule ${n} a cyclic peptide? *}
{+ choice: true false +}
{===>} cyclicpept_mol${n}=${partners[n]['cyclic']};
{* Is molecule ${n} a shape? *}
{+ choice: true false +}
{===>} shape_mol${n}=${partners[n]['shape']};
{* Coarse grained molecule? *}
{+ choice: true false +}
{===>} cg_mol${n}=${partners[n]['cg']};
% endfor

{* Remove non-polar hydrogens? *}
{+ choice: true false +}
{===>} delenph=${delenph};

{* HADDOCK directory *}
{*  the absolute path of the HADDOCK program files *}
{===>} haddock_dir=${haddock_dir};

{* Logfile directory *}
{* specify a directory for the large CNS log files *}
{===>} temptrash_dir=${temptrash_dir};

{==================== histidine patches =====================}
## partners[p]['his_patch'] data structure example -> {12: 'HISE', 23: 'HISD', 175:'HIS+'}
## Output: {===>} 1_hisd_resid_1=12;
<%
    numhisd_all = []
    numhise_all = []
    for n in range(0, ncomp):
        numhisd_all.append([])
        numhise_all.append([])
        if partners[n+1]['his_patch']:
            for k, v in partners[n+1]['his_patch'].items():
                if v == "HISD":
                    numhisd_all[n].append(k)
                elif v == "HISE":
                    numhise_all[n].append(k)
%>

{==================== histidine patches =====================}
{* Automatically define histidine protonation state based on energetics *}
{===>} autohis=false;

{* Patch to change doubly protonated HIS to singly protonated histidine (HD1) *}
{* just give the residue number of the histidines for the HISD patch, set them to zero if you don't want them *}
<% c = 1 %>
% for n in range(0, ncomp):
{* Number of HISD for molecule ${c} *}
numhisd_${n+1}=${len(numhisd_all[n])};
    % for hisd in numhisd_all[n]:
{===>} hisd_${n+1}_${c}=${hisd};
    <% c += 1 %>
    % endfor
% endfor

{* Patch to change doubly protonated HIS to singly protonated histidine (HE2) *}
{* just give the residue number of the histidines for the HISE patch, set them to zero if you don't want them *}
<% c = 1 %>
% for n in range(0, ncomp):
{* Number of HISE for molecule ${c} *}
numhise_${n+1}=${len(numhise_all[n])};
    % for hise in numhise_all[n]:
{===>} hise_${n+1}_${c}=${hise};
    <% c += 1 %>
    % endfor
% endfor

{========= Definition of semi-flexible interface ============}
{* Define the interface of each molecule.*}
{* Side-chains and backbone of these residues will be allowed to move during semi-flexible refinement*}
{* Distance cutoff in A for the automatic definition of flexible segments based on intermolecular residues contacts *}
{===>} flcut_nb=${flcut_nb};
## semi_flex data structure example -> {'A': [[12,20], [45,60]], 'B': None}
## Output: {===>} A_start_seg_1=12; \n{===>} A_end_seg_1=20; ...
% for n in range(0,ncomp):
{* number of semi-flexible segments for molecule  ${n+1} (-1 for automated mode) *}
{* Note that current max is 10 (edit the run.cns to add more segments *}
    % if partners[n+1]['semi_flex']:
{===>} nseg_${n+1}=${len(partners[n+1]['semi_flex'])};

{* Residues of molecule  ${n+1} at interface *}
{+ table: rows=10 "segment 1" "segment 2" "segment 3" "segment 4" "segment 5" "segment 6" "segment 7" "segment 8" "segment 9" "segment 10" cols=2 "Start residue" "End residue" +}\
        <% c = 1 %>
        % for seq in partners[n+1]['semi_flex']:
{===>} start_seg_${n+1}_${c}=${seq[0]};
{===>} end_seg_${n+1}_${c}=${seq[1]};
            <% c+= 1 %>
        % endfor
    %else:
{===>} nseg_${n+1}=-1;
    % endif
% endfor

{=========== Definition of fully flexible segments ==========}
{* Define the fully flexible segment of each molecule.*}
{* These segments will be allowed to move at all stages of it1 *}

## fully_flex data structure example -> {1: [12,20], [45,60]}, 2: None}
## Output: {===>} 1_start_fle_1=12; \n{===>} 1_end_fle_1=20; ...

% for n in range(0,ncomp):
{* Number of fully flexible segments for molecule  ${n+1}            *}
{* Note that current max is 5 (edit the run.cns to add more segments     *}
    % if partners[n+1]['fully_flex']:
{===>} nfle_${n+1}=${len(partners[n+1]['fully_flex'])};

{* Fully flexible segments of molecule  ${n+1} *}
{+ table: rows=5 "segment 1" "segment 2" "segment 3" "segment 4" "segment 5" cols=2 "Start residue" "End residue" +}\
        <% c = 1 %>
        % for seq in partners[n+1]['fully_flex']:
{===>} start_fle_${n+1}_${c}=${seq[0]};
{===>} end_fle_${n+1}_${c}=${seq[1]};\
            <% c += 1 %>
        % endfor
    % else:
{===>} nfle_${n+1}=0;
    %endif
% endfor

{==================== membrane positioning restraints  ==================}
{* Do you want to use membrane positioning restraints ? *}
{+ choice: true false +}
{===>} zres_on=${str(zres_on).lower()};

{* Number of membrane positioning restrained segments *}
{===>} numzres=${numzres};
% if zres_on:
{* Force constant for membrane positioning restraints ? *}
{===>} kzres=${kzres};

{* Maximum z value for membrane positioning restraints ? *}
{===>} zresmax=${zresmax};

{* Minimum z value for membrane positioning restraints ? *}
{===>} zresmin=${zresmin};

{* Define the segment for membrane positioning restraints *}
{+ table: rows=10 "seg 1" "seg 2" "seg 3" "seg 4" "seg 5" "seg 6" "seg 7" "seg 8" "seg 9" "seg 10" cols=4 "Start res seg1" "End res seg1" "Segid seg1" "inside/outside" +}
    % for n in range(1, numzres+1):
{===>} zres_sta_${n}=${zres[n]['sta']};
{===>} zres_end_${n}=${zres[n]['end']};
{===>} zres_seg_${n}=${zres[n]['seg']};
{+ choice: "inside" "outside"+}
{===>} zres_type_${n}=${zres[n]['zres_type']};
    % endfor
% endif

{====================== NCS restraints  =====================}
{* Do you want to use NCS restraints? *}
{+ choice: true false +}
{===>} ncs_on=${str(ncs_on).lower()};

{* Number of NCS pairs *}
{===>} numncs=${numncs};
% if ncs_on:
{* Force constant for NCS restraints *}
{===>} kncs=${kncs};

{* Define the segments pairs for NCS restraints *}
{+ table: rows=5 "pair 1" "pair 2" "pair 3" "pair 4" "pair 5" cols=6 "Start res seg1" "End res seg1" "Segid seg1" "Start res seg2" "End res seg2" "Segid seg2" +}
    % for n in range(numncs):
{===>} ncs_sta1_${n+1}="${ncs[n]['sta1']}";
{===>} ncs_end1_${n+1}="${ncs[n]['end1']}";
{===>} ncs_seg1_${n+1}="${ncs[n]['seg1']}";
{===>} ncs_sta2_${n+1}="${ncs[n]['sta2']}";
{===>} ncs_end2_${n+1}="${ncs[n]['end2']}";
{===>} ncs_seg2_${n+1}="${ncs[n]['seg2']}";
    % endfor
% endif

{==================== Symmetry restraints  ==================}
{* Do you want to use symmetry restraints ? *}
{+ choice: true false +}
{===>} sym_on=${str(sym_on).lower()};

{* Force constant for symmetry restraints ? *}
{===>} ksym=${ksym};
{* Number of C2 symmetry pairs *}
{===>} numc2sym=${numc2sym};
% if sym_on and numc2sym > 0:
{* Define the segment pairs C2 symmetry restraints *}
{+ table: rows=10 "pair 1" "pair 2" "pair 3" "pair 4" "pair 5" "pair 6" "pair 7" "pair 8" "pair 9" "pair 10" cols=6 "Start res seg1" "End res seg1" "Segid seg1" "Start res seg2" "End res seg2" "Segid seg2" +}
    % for n in range(numc2sym):
{===>} c2sym_sta1_${n+1}="${c2sym[n]['sta1']}";
{===>} c2sym_end1_${n+1}="${c2sym[n]['end1']}";
{===>} c2sym_seg1_${n+1}="${c2sym[n]['seg1']}";
{===>} c2sym_sta2_${n+1}="${c2sym[n]['sta2']}";
{===>} c2sym_end2_${n+1}="${c2sym[n]['end2']}";
{===>} c2sym_seg2_${n+1}="${c2sym[n]['seg2']}";
    % endfor
% endif
{* Number of C3 symmetry triples*}
{===>} numc3sym=${numc3sym};
% if sym_on and numc3sym > 0:
{* Define the segment triples for C3 symmetry restraints *}
{+ table: rows=2 "triple 1" "triple 2" cols=9 "Start res seg1" "End res seg1" "Segid seg1" "Start res seg2" "End res seg2" "Segid seg2" "Start res seg3" "End res seg3" "Segid seg3" +}
    % for n in range(numc3sym):
{===>} c3sym_sta1_${n+1}="${c3sym[n]['sta1']}";
{===>} c3sym_end1_${n+1}="${c3sym[n]['end1']}";
{===>} c3sym_seg1_${n+1}="${c3sym[n]['seg1']}";
{===>} c3sym_sta2_${n+1}="${c3sym[n]['sta2']}";
{===>} c3sym_end2_${n+1}="${c3sym[n]['end2']}";
{===>} c3sym_seg2_${n+1}="${c3sym[n]['seg2']}";
{===>} c3sym_sta3_${n+1}="${c3sym[n]['sta3']}";
{===>} c3sym_end3_${n+1}="${c3sym[n]['end3']}";
{===>} c3sym_seg3_${n+1}="${c3sym[n]['seg3']}";
    % endfor
% endif
{* Number of S3 symmetry triples*}
{===>} nums3sym=${nums3sym};
% if sym_on and nums3sym > 0:
{* Define the segment triples for S3 symmetry restraints *}
{+ table: rows=4 "triple 1" "triple 2" "triple 3" "triple 4" cols=9 "Start res seg1" "End res seg1" "Segid seg1" "Start res seg2" "End res seg2" "Segid seg2" "Start res seg3" "End res seg3" "Segid seg3" +}
    % for n in range(nums3sym):
{===>} s3sym_sta1_${n+1}="${s3sym[n]['sta1']}";
{===>} s3sym_end1_${n+1}="${s3sym[n]['end1']}";
{===>} s3sym_seg1_${n+1}="${s3sym[n]['seg1']}";
{===>} s3sym_sta2_${n+1}="${s3sym[n]['sta2']}";
{===>} s3sym_end2_${n+1}="${s3sym[n]['end2']}";
{===>} s3sym_seg2_${n+1}="${s3sym[n]['seg2']}";
{===>} s3sym_sta3_${n+1}="${s3sym[n]['sta3']}";
{===>} s3sym_end3_${n+1}="${s3sym[n]['end3']}";
{===>} s3sym_seg3_${n+1}="${s3sym[n]['seg3']}";
    % endfor
% endif
{* Number of C4 symmetry quadruples *}
{===>} numc4sym=${numc4sym};
% if sym_on and numc4sym > 0:
{* Define the segment quadruples for C4 symmetry restraints *}
{+ table: rows=2 "quadruples 1" "quadruples 2" cols=12 "Start res seg1" "End res seg1" "Segid seg1" "Start res seg2" "End res seg2" "Segid seg2" "Start res seg3" "End res seg3" "Segid seg3" "Start res seg4" "End res seg4" "Segid seg4" +}
    % for n in range(numc4sym):
{===>} c4sym_sta1_${n+1}="${c4sym[n]['sta1']}";
{===>} c4sym_end1_${n+1}="${c4sym[n]['end1']}";
{===>} c4sym_seg1_${n+1}="${c4sym[n]['seg1']}";
{===>} c4sym_sta2_${n+1}="${c4sym[n]['sta2']}";
{===>} c4sym_end2_${n+1}="${c4sym[n]['end2']}";
{===>} c4sym_seg2_${n+1}="${c4sym[n]['seg2']}";
{===>} c4sym_sta3_${n+1}="${c4sym[n]['sta3']}";
{===>} c4sym_end3_${n+1}="${c4sym[n]['end3']}";
{===>} c4sym_seg3_${n+1}="${c4sym[n]['seg3']}";
{===>} c4sym_sta4_${n+1}="${c4sym[n]['sta4']}";
{===>} c4sym_end4_${n+1}="${c4sym[n]['end4']}";
{===>} c4sym_seg4_${n+1}="${c4sym[n]['seg4']}";
    % endfor
% endif
{* Number of C5 symmetry *}
{===>} numc5sym=${numc5sym};
% if sym_on and numc5sym > 0:
{* Define the segments for C5 symmetry restraints *}
{+ table: rows=5 "Segment1" "Segment2" "Segment3" "Segment4" "Segment5" cols=3 "Start residue" "End residue" "Segid" +}
    % for n in range(numc5sym):
{===>} c5sym_sta1_${n+1}="${c5sym[n]['sta1']}";
{===>} c5sym_end1_${n+1}="${c5sym[n]['end1']}";
{===>} c5sym_seg1_${n+1}="${c5sym[n]['seg1']}";
{===>} c5sym_sta2_${n+1}="${c5sym[n]['sta2']}";
{===>} c5sym_end2_${n+1}="${c5sym[n]['end2']}";
{===>} c5sym_seg2_${n+1}="${c5sym[n]['seg2']}";
{===>} c5sym_sta3_${n+1}="${c5sym[n]['sta3']}";
{===>} c5sym_end3_${n+1}="${c5sym[n]['end3']}";
{===>} c5sym_seg3_${n+1}="${c5sym[n]['seg3']}";
{===>} c5sym_sta4_${n+1}="${c5sym[n]['sta4']}";
{===>} c5sym_end4_${n+1}="${c5sym[n]['end4']}";
{===>} c5sym_seg4_${n+1}="${c5sym[n]['seg4']}";
{===>} c5sym_sta5_${n+1}="${c5sym[n]['sta5']}";
{===>} c5sym_end5_${n+1}="${c5sym[n]['end5']}";
{===>} c5sym_seg5_${n+1}="${c5sym[n]['seg5']}";
    % endfor
% endif

{* Number of c6 symmetry *}
{===>} numc6sym=${numc6sym};
% if sym_on and numc6sym > 0:
{* Define the segments for c6 symmetry restraints *}
{+ table: rows=6 "Segment1" "Segment2" "Segment3" "Segment4" "Segment5" "Segment6" cols=3 "Start residue" "End residue" "Segid" +}
    % for n in range(numc6sym):
{===>} c6sym_sta1_${n+1}="${c6sym[n]['sta1']}";
{===>} c6sym_end1_${n+1}="${c6sym[n]['end1']}";
{===>} c6sym_seg1_${n+1}="${c6sym[n]['seg1']}";
{===>} c6sym_sta2_${n+1}="${c6sym[n]['sta2']}";
{===>} c6sym_end2_${n+1}="${c6sym[n]['end2']}";
{===>} c6sym_seg2_${n+1}="${c6sym[n]['seg2']}";
{===>} c6sym_sta3_${n+1}="${c6sym[n]['sta3']}";
{===>} c6sym_end3_${n+1}="${c6sym[n]['end3']}";
{===>} c6sym_seg3_${n+1}="${c6sym[n]['seg3']}";
{===>} c6sym_sta4_${n+1}="${c6sym[n]['sta4']}";
{===>} c6sym_end4_${n+1}="${c6sym[n]['end4']}";
{===>} c6sym_seg4_${n+1}="${c6sym[n]['seg4']}";
{===>} c6sym_sta5_${n+1}="${c6sym[n]['sta5']}";
{===>} c6sym_end5_${n+1}="${c6sym[n]['end5']}";
{===>} c6sym_seg5_${n+1}="${c6sym[n]['seg5']}";
{===>} c6sym_sta6_${n+1}="${c6sym[n]['sta6']}";
{===>} c6sym_end6_${n+1}="${c6sym[n]['end6']}";
{===>} c6sym_seg6_${n+1}="${c6sym[n]['seg6']}";
    % endfor
% endif

{=========================== Distance restraints  ========================}
{* Turn on/off and energy constants for distance restraints *}
{+ table: rows=3 "distances" "AIR (ambig)" "hbonds" cols=6 "firstIteration" "lastIteration" "hot" "cool1" "cool2" "cool3"+}

{===>} unamb_firstit=${unamb_firstit};
{===>} unamb_lastit=${unamb_lastit};
{===>} unamb_hot=${unamb_hot};
{===>} unamb_cool1=${unamb_cool1};
{===>} unamb_cool2=${unamb_cool2};
{===>} unamb_cool3=${unamb_cool3};
{===>} amb_firstit=${amb_firstit};
{===>} amb_lastit=${amb_lastit};
{===>} amb_hot=${amb_hot};
{===>} amb_cool1=${amb_cool1};
{===>} amb_cool2=${amb_cool2};
{===>} amb_cool3=${amb_cool3};
{===>} hbond_firstit=${hbond_firstit};
{===>} hbond_lastit=${hbond_lastit};
{===>} hbond_hot=${hbond_hot};
{===>} hbond_cool1=${hbond_cool1};
{===>} hbond_cool2=${hbond_cool2};
{===>} hbond_cool3=${hbond_cool3};

{* Do you want to randomly exclude a fraction of the ambiguous restraints (AIRs)? *}
{+ choice: true false +}
{===>} noecv=${noecv};

{* Number of partitions for random exclusion (%excluded=100/number of partitions)? *}
{===>} ncvpart=${ncvpart};

{* Do you want to use hydrogen bond restraints? *}
{+ choice: true false +}
{===>} hbonds_on=${hbonds_on};

{* Do you want to define randomly ambiguous interaction restraints from accessible residues? *}
{* Only residues in the defined flexible segments will be considered *}
{* Note that this option is exclusive with any other distance restraints and only for it0    *}
{+ choice: true false +}
{===>} ranair=${ranair};

{* Do you want to define center of mass restraints to enforce contact between the molecules? *}
{* Note that these are only active during it0 and it1 *}
{+ choice: true false +}
{===>} cmrest=${cmrest};

{* Define tight CM restraints? *}
{+ choice: true false +}
{===>} cmtight=${cmtight};

{* Force constant for center of mass restraints *}
{===>} kcont=${kcont};

{* Do you want to define surface contact restraints to enforce contact between the molecules? *}
{* Note that these are only active during it0 and it1 *}
{+ choice: true false +}
{===>} surfrest=${surfrest};

{* Force constant for surface contact restraints *}
{===>} ksurf=${ksurf};

{ Use automated distance restraints weighting }
{ choice: true false }
air_scaling=false;
{ Define the number of distance restraints for automated weighting }
tot_unamb=25;
{ Define the number of AIR restraints for automated weighting }
tot_amb=0;
{ potential shape }
mrswi_hot=0.5;
mrswi_cool1=0.5;
mrswi_cool2=0.5;
mrswi_cool3=0.5;
rswi_hot=0.5;
rswi_cool1=0.5;
rswi_cool2=0.5;
rswi_cool3=0.5;
masy_hot=-1.0;
masy_cool1=-1.0;
masy_cool2=-0.1;
masy_cool3=-0.1;
asy_hot=1.0;
asy_cool1=1.0;
asy_cool2=0.1;
asy_cool3=0.1;

{=========================== radius of gyration restraint  ============}
{* Turn on/off and energy constants for Rg restraints *}
{* Do you want to define a radius of gyration restraint (e.g. from SAXS)? *}
{+ choice: true false +}
{===>} rgrest=${str(rgrest).lower()};

% if rgrest:
{* Radius of gyration *}
{===>} rgtarg=${rgtarg};

{* Force constant for radius of gyration restraint *}
{===>} krg_hot=${krg_hot};
{===>} krg_cool1=${krg_cool1};
{===>} krg_cool2=${krg_cool2};
{===>} krg_cool3=${krg_cool3};

{* Atom selections for the radius of gyration restraint *}
{===>} rgsele="${rgsele}";
% endif

{======================DNA-RNA restraints ============================}
{* Use DNA/RNA restraints (dna-rna_restraints.def in data/sequence)? *}
{+ choice: true false +}
{===>} dnarest_on=${str(dnarest_on).lower()};

{=========================== dihedrals restraints ====================}
{* energy constants *}
{+ table: rows=1 "dihedrals" cols=5 "use?" "hot" "cool1" "cool2" "cool3" +}

{+ choice: true false +}
{===>} dihedrals_on=${str(dihedrals_on).lower()};
{===>} dihedrals_hot=${dihedrals_hot};
{===>} dihedrals_cool1=${dihedrals_cool1};
{===>} dihedrals_cool2=${dihedrals_cool2};
{===>} dihedrals_cool3=${dihedrals_cool3};

{* Automatically define backbone dihedral angle restraints from structure? *}
{+ choice: none all alpha alphabeta +}
{===>} ssdihed=${ssdihed};
{===>} error_dih=${error_dih};

{=========================== residual dipolar couplings ======================}

{* Parameters *}
{+ table: rows=5 "class1" "class2" "class3" "class4" "class5"
          cols=25 "type" "firstIt" "lastIt" "Ksani<br>(hot)" "Ksani<br>(cool1)" "Ksani<br>(cool2)" "Ksani<br>(cool3)" "R" "D"
 "Kvean<br>(ini_bor_hot)" "Kvean<br>(fin_bor_hot)"
 "Kvean<br>(ini_bor_cool1)" "Kvean<br>(fin_bor_cool1)"
 "Kvean<br>(ini_bor_cool2)" "Kvean<br>(fin_bor_cool2)"
 "Kvean<br>(ini_bor_cool3)" "Kvean<br>(fin_bor_cool3)"
 "Kvean<br>(ini_cen_hot)" "Kvean<br>(fin_cen_hot)"
 "Kvean<br>(ini_cen_cool1)" "Kvean<br>(fin_cen_cool1)"
 "Kvean<br>(ini_cen_cool2)" "Kvean<br>(fin_cen_cool2)"
 "Kvean<br>(ini_cen_cool3)" "Kvean<br>(fin_cen_cool3)"+}
{+ choice: "NO" "SANI" "XRDC" "VANGLE" +}
{===>} numrdc = ${len(rdc)};
 <% c = 1 %>
% for r in rdc:
{===>} rdc_choice_${c}="${r['choice']}";
{===>} rdc_firstIt_${c}=${r['firstit']};
{===>} rdc_lastIt_${c}=${r['lastit']};
{===>} rdc_hot_${c}=${r['hot']};
{===>} rdc_cool1_${c}=${r['cool1']};
{===>} rdc_cool2_${c}=${r['cool2']};
{===>} rdc_cool3_${c}=${r['cool3']};
{===>} rdc_r_${c}=${r['r']};
{===>} rdc_d_${c}=${r['d']};
{===>} ini_bor_hot_${c}=${r['ini_bor_hot']};
{===>} fin_bor_hot_${c}=${r['fin_bor_hot']};
{===>} ini_bor_cool1_${c}=${r['ini_bor_cool1']};
{===>} fin_bor_cool1_${c}=${r['fin_bor_cool1']};
{===>} ini_bor_cool2_${c}=${r['ini_bor_cool2']};
{===>} fin_bor_cool2_${c}=${r['fin_bor_cool2']};
{===>} ini_bor_cool3_${c}=${r['ini_bor_cool3']};
{===>} fin_bor_cool3_${c}=${r['fin_bor_cool3']};
{===>} ini_cen_hot_${c}=${r['ini_cen_hot']};
{===>} fin_cen_hot_${c}=${r['fin_cen_hot']};
{===>} ini_cen_cool1_${c}=${r['ini_cen_cool1']};
{===>} fin_cen_cool1_${c}=${r['fin_cen_cool1']};
{===>} ini_cen_cool2_${c}=${r['ini_cen_cool2']};
{===>} fin_cen_cool2_${c}=${r['fin_cen_cool2']};
{===>} ini_cen_cool3_${c}=${r['ini_cen_cool3']};
{===>} fin_cen_cool3_${c}=${r['fin_cen_cool3']};\
    <% c += 1 %>
% endfor

{=========================== pseudo contact shifts ===========================}

{* Parameters *}
{+ table: rows=10 "class1" "class2" "class3" "class4" "class5" "class6" "class7" "class8" "class9" "class10"
          cols=9 "type" "firstIt" "lastIt" "Kpcs<br>(hot)" "Kpcs<br>(cool1)" "Kpcs<br>(cool2)" "Kpcs<br>(cool3)" "R" "D" +}
{+ choice: "NO" "XPCS" +}
{===>} numpcs = ${len(pcs)};
<% c = 1 %>
% for p in pcs:
{===>} pcs_choice_${c}="XPCS";
{===>} pcs_firstIt_${c}=${p['firstit']};
{===>} pcs_lastIt_${c}=${p['lastit']};
{===>} pcs_hot_${c}=${p['hot']};
{===>} pcs_cool1_${c}=${p['cool1']};
{===>} pcs_cool2_${c}=${p['cool2']};
{===>} pcs_cool3_${c}=${p['cool3']};
{===>} pcs_r_${c}=${p['r']};
{===>} pcs_d_${c}=${p['d']};\
    <% c += 1 %>
% endfor

{=========================== relaxation data ======================}
{* Parameters *}
{+ table: rows=5 "class1" "class2" "class3" "class4" "class5"
          cols=12 "type" "firstIt" "lastIt" "Kdani(hot)" "Kdani(cool1)" "Kdani(cool2)" "Kdani(cool3)" "Correlation time" "D" "R" "H frequency" "N frequency" +}
{+ choice: "NO" "DANI" +}
{===>} numdani = ${len(dan)};
<% c = 1 %>
% for d in dan:
{===>} dan_choice_${c}="DANI";
{===>} dan_firstIt_${c}=${d['firstit']};
{===>} dan_lastIt_${c}=${d['lastit']};
{===>} dan_hot_${c}=${d['hot']};
{===>} dan_cool1_${c}=${d['cool1']};
{===>} dan_cool2_${c}=${d['cool2']};
{===>} dan_cool3_${c}=${d['cool3']};
{===>} dan_tc_${c}=${d['tc']};
{===>} dan_anis_${c}=${d['anis']};
{===>} dan_r_${c}=${d['r']};
{===>} dan_wh_${c}=${d['wh']};
{===>} dan_wn_${c}=${d['wn']};\
    <% c += 1 %>
% endfor

{========================== Cryo-EM parameters ============================}

{* Centroid definitions *}
{+ choice: true false +}
{===>} centroid_rest=${str(centroid_rest).lower()};
{===>} centroid_kscale=${centroid_kscale};

% if centroid_rest:
    % for n in range(1, ncomp+1):
{* Placement of centroids in absolute coordinates *}
{===>} xcom_${n}=${centroids[n]['coor'][0]};
{===>} ycom_${n}=${centroids[n]['coor'][1]};
{===>} zcom_${n}=${centroids[n]['coor'][2]};
{* Are the centroid retraints ambiguous *}
{+ choice: true false +}
{===>} ambi_${n}=${centroids[n]['ambig']};
    % endfor
% endif

{* Density/XREF restraints *}
{+ choice: true false +}
{===>} em_rest=${em_rest};
{===>} em_kscale=${em_kscale};
{+ choice: true false +}
{===>} em_it0=${em_it0};
{+ choice: true false +}
{===>} em_it1=${em_it1};
{+ choice: true false +}
{===>} em_itw=${em_itw};

{* Resolution of data in angstrom *}
{===>} em_resolution=${em_resolution};

{* Density parameters *}
{* Number of voxels in each dimension *}
{===>} nx=${voxels_num[0]};
{===>} ny=${voxels_num[1]};
{===>} nz=${voxels_num[2]};

{* Length of each dimension in angstrom *}
{===>} xlength=${voxels_dim[0]};
{===>} ylength=${voxels_dim[1]};
{===>} zlength=${voxels_dim[2]};

{* Cryo-EM scoring weights *}
{===>} w_lcc_0=${weights['lcc'][0]};
{===>} w_lcc_1=${weights['lcc'][1]};
{===>} w_lcc_2=${weights['lcc'][2]};

{===================== topology and parameter files ======================}

% for n in range(1, ncomp+1):
{* topology file for molecule  ${n} *}
{===>} prot_top_mol${n}="${partners[n]['top_file']}";
{* linkage file for molecule  ${n} *}
{===>} prot_link_mol${n}="${partners[n]['link_file']}";
{* energy parameter file for molecule  ${n} *}
{===>} prot_par_mol${n}="${partners[n]['par_file']}";
% endfor

{* type of non-bonded parameters *}
{* specify the type of non-bonded interaction *}
{+ choice: "PROLSQ" "PARMALLH6" "PARALLHDG" "OPLSX" +}
{===>} par_nonbonded="${par_nonbonded}";

{============coarse graining topology and parameter files ==================}

% for n in range(1, ncomp+1):
{* topology file for molecule  ${n} *}
{===>} prot_cg_top_mol${n}="protein-CG-Martini-2-2.top";
{* linkage file for molecule  ${n} *}
{===>} prot_cg_link_mol${n}="protein-CG-Martini-2-2.link";
{* energy parameter file for molecule  ${n} *}
{===>} prot_cg_par_mol${n}="protein-CG-Martini-2-2.param";
% endfor

{===================== energy and interaction parameters ==================}

{ Do you want to include dihedral angle energy terms? }
{ choice: true false }
dihedflag=true;

{* Do you want to include the electrostatic energy term for docking? *}
{* Note that it will be automatically included in the solvent refinement *}

{* Include electrostatic during rigid body docking (it0)? *}
{+ choice: true false +}
{===>} elecflag_0=${elecflag_0};
{* Include electrostatic during semi-flexible SA (it1)? *}
{+ choice: true false +}
{===>} elecflag_1=${elecflag_1};

{* Give the epsilon constant for the electrostatic energy term in it0 *}
{* Note that for explicit solvent refinement cdie with epsilon=1 is used *}
{===>} epsilon_0=${epsilon_0};

{* Give the epsilon constant for the electrostatic energy term in it1 *}
{* Note that for explicit solvent refinement cdie with epsilon=1 is used *}
{===>} epsilon_1=${epsilon_1};

{* Use constant (cdie) or distance-dependent (rdie) dielectric in it0? *}
{+ choice: cdie rdie +}
{===>} dielec_0=${dielec_0};

{* Use constant (cdie) or distance-dependent (rdie) dielectric in it1? *}
{+ choice: cdie rdie +}
{===>} dielec_1=${dielec_1};

{* Scaling of intermolecular interactions for rigid body EM*}
{===>} inter_rigid=${inter_rigid};

{* Scaling of intermolecular interactions for semi-flexible SA*}
{+ table: rows=3 "Rigid body dynamic " "SA with flexible side-chains (cool2)" "SA with flexible backbone and side-chains (cool3)"
          cols=2 "Init value" "Final value" +}
{===>} init_rigid=${init_rigid};
{===>} fin_rigid=${fin_rigid};
{===>} init_cool2=${init_cool2};
{===>} fin_cool2=${fin_cool2};
{===>} init_cool3=${init_cool3};
{===>} fin_cool3=${fin_cool3};

{* Interaction matrix for non-bonded interactions*}
{+ table: rows=6 "Mol 1" "Mol 2" "Mol 3" "Mol 4" "Mol 5" "Mol 6" "Mol 7" "Mol 8" "Mol 9" Mol 10" "Mol 11" Mol 12" "Mol 13" "Mol 14" "Mol 15" "Mol 16" "Mol 17" "Mol 18" "Mol 19" "Mol 20"
          cols=6 "Mol 1" "Mol 2" "Mol 3" "Mol 4" "Mol 5" "Mol 6" "Mol 7" "Mol 8" "Mol 9" Mol 10" "Mol 11" Mol 12" "Mol 13" "Mol 14" "Mol 15" "Mol 16" "Mol 17" "Mol 18" "Mol 19" "Mol 20" +}
% for i in range(1, nb_mol_max+1):
    % for j in range(1, nb_mol_max+1):
        % if j >= i:
{===>} int_${i}_${j}=${inter_mat[i-1][j-1]};
        % else:
{===>} int_${i}_${j}="N.A.";
        % endif
    % endfor
% endfor

{===================== Number of structures to dock =======================}
{* Setting for the rigid-body (it0) and semi-flexible refinement (it1) *}

{* number of structures for rigid body docking *}
{===>} structures_0=${structures_0};
       keepstruct_0=&structures_0;
{* number of structures for refinement *}
{===>} structures_1=${structures_1};
       keepstruct_1=&structures_1;
       keepstruct_2=&structures_1;
{* number of structures to be analysed*}
{===>} anastruc_1=${anastruc_1};
       anastruc_0=&anastruc_1;
       anastruc_2=&anastruc_1;

{* - *}

{* Sampling of symmetry related solutions                       *}

{* Sample 180 degrees rotated solutions during rigid body EM?   *}
{+ choice: true false +}
{===>} rotate180_it0=${rotate180_0};

{* Sample 180 degrees rotated solutions during semi-flexible SA?*}
{+ choice: true false +}
{===>} rotate180_it1=${rotate180_1};


{=========================== DOCKING protocol =============================}
{* Cross-dock all combinations in the ensembles of starting structures? *}
{* Turn off this option if you only want to dock structure 1 of ensemble A *}
{*   to structure 1 of ensemble B, structure 2 to structure 2, etc. *}
{+ choice: true false +}
{===>} crossdock=${crossdock};

{* Randomize starting orientations? *}
{+ choice: true false +}
{===>} randorien=${randorien};

{* Expand starting orientations? *}
{+ choice: true false +}
{===>} expand=${expand};

{* Expansion percentage *}
{===>} expansion=${expansion};

{* Random rotation angle *}
{===>} randangle=${randangle};

{* Rebuild missing atoms in the context of the complex? (refinement mode) *}
{+ choice: true false +}
{===>} rebuildcplx=false;

{* Perform initial rigid body minimisation? *}
{+ choice: true false +}
{===>} rigidmini=${rigidmini};

{* Allow translation in rigid body minimisation? *}
{+ choice: true false +}
{===>} rigidtrans=${rigidtrans};

{* Number of trials for rigid body minimisation? *}
{===>} ntrials=${ntrials};

{* initial seed for random number generator *}
{* change to get different initial velocities *}
{===>} iniseed=${iniseed};

{* temperature for rigid body high temperature TAD *}
{===>} tadhigh_t=${tadhigh_t};

{* initial temperature for rigid body first TAD cooling step *}
{===>} tadinit1_t=${tadinit1_t};

{* final temperature after first cooling step *}
{===>} tadfinal1_t=${tadfinal1_t};

{* initial temperature for second TAD cooling step with flexible side-chain at the inferface *}
{===>} tadinit2_t=${tadinit2_t};

{* finale temperature after second cooling step *}
{===>} tadfinal2_t=${tadfinal2_t};

{* initial temperature for third TAD cooling step with fully flexible interface *}
{===>} tadinit3_t=${tadinit3_t};

{* finale temperature after third cooling step *}
{===>} tadfinal3_t=${tadfinal3_t};

{* time step *}
{===>} timestep=${timestep};
{* factor for timestep in TAD *}
{===>} tadfactor=${tadfactor};

{* Number of EM steps for translational minimisation? *}
{===>} emstepstrans=${emstepstrans};

{* number of MD steps for rigid body high temperature TAD *}
{===>} initiosteps=${initiosteps};

{* number of MD steps during first rigid body cooling stage *}
{===>} cool1_steps=${cool1_steps};

{* number of MD steps during second cooling stage with flexible side-chains at interface *}
{===>} cool2_steps=${cool2_steps};

{* number of MD steps during third cooling stage with fully flexible interface *}
{===>} cool3_steps=${cool3_steps};


{======================= Solvated rigid body docking=======================}
{* perform solvated docking ? *}
{+ choice: true false +}
{===>} waterdock=${waterdock};

{* which method to use for solvating? *}
{* db: database-based (recommended), restraints: for restrained solvating to amino-acid most often forming
water mediated contacts and blank (""): for uniform waterlayer *}
{+ choice: "db" "restraints" "" +}
{===>} solvate_method="${solvate_method}";

{* which propensity database to use? *}
{* statistical: based on an analysis of water-mediated contacts in the PDB, kyte-doolittle: based on the Kyte-Doolittle hydrophobicity scalte *}
{+ choice: "statistical" "kytedoolittle" +}
{===>} db_method="${db_method}";

{* initial cutoff for restraints solvating method *}
{* all waters further away from a highly occuring water solvated residue will be removed in the generation
of the initial solvation shell *}
{===>} water_restraint_initial=${water_restraint_initial};

{* cutoff for restraints solvating method *}
{* upper distance limit for defining distance restraints between water and amino-acids often found to be
involved in water-mediated contacts *}
{===>} water_restraint_cutoff=${water_restraint_cutoff};

{* force constant for restrainted solvating method *}
{===>} water_restraint_scale=${water_restraint_scale};

{* fraction of water to keep *}
{* this is the fraction of all interface water after the initial rigid body docking that will be kept
(note that more waters might be removed if the interaction energy is unfavorable  *}
{===>} water_tokeep=${water_tokeep};

{* fraction of water around DNA to keep *}
{* this is the fraction of interface water involving DNA phoshpates after the initial rigid body docking that will be kept
(note that more waters might be removed if the interaction energy is unfavorable  *}
{===>} dnap_water_tokeep=${dnap_water_tokeep};

{* random fraction to be added to the fraction of water to keep *}
{===>} water_randfrac=${water_randfrac};

{* water-protein surface-cutoff *}
{* waters further away than this cutoff distance from any component of the complex will be removed *}
{===>} water_surfcutoff=${water_surfcutoff};

{* do some water analysis *}
{+ choice: true false +}
{===>} water_analysis=${water_analysis};

{* allows translation of water molecules during rigid-body docking, true or false: *}
{+ choice: true false +}
{===>} transwater=${transwater};

{* number of different initial solvation shells to generate *}
{===>} waterensemble=${waterensemble};


{==================== final explicit solvent refinement  ==================}
{* Do you want to refine your docking models in explicit solvent? *}
{+ choice: "yes" "no" +}
{===>} firstwater="${firstwater}";

{* Build explicit solvent shell? (Can be turned off the large molecules or when morphing CG to AA models) *}
{* Only EM will then be performed                                                                         *}
{+ choice: true false +}
{===>} solvshell=${solvshell};

{* Which solvent do you want to use? *}
{+ choice: "water" "dmso" +}
{===>} solvent="${solvent}";

{* number of structures for the explicit solvent refinement *}
{* the n best structures will be refined                    *}
{===>} waterrefine=${waterrefine};
       structures_2=&waterrefine;

{* number of steps for heating phase (100, 200, 300K)?      *}
{===>} waterheatsteps=${waterheatsteps};

{* number of steps for 300K sampling phase?                 *}
{===>} watersteps=${watersteps};

{* number of steps for cooling phase (300, 200, 100K)?      *}
{===>} watercoolsteps=${watercoolsteps};

{* write additional PDB files including solvent ?           *}
{+ choice: true false +}
{===>} keepwater=${keepwater};

{================================ Scoring =================================}
{* Settings for the scoring of the docking solutions *}

{* Define the weights for the various terms for the sorting of structures (scoring) *}
{+ table: rows=14 "Evdw" "Eelec" "Eair" "Erg" "Esani" "Exrdc" "Expcs" "Edani" "Evean" "Ecdih" "Esym" "BSA" "dEint" "Edesolv"
          cols=3 "Rigid body EM" "semi-flexible SA" "Water refinement" +}
% for w in ('vdw', 'elec', 'dist', 'rg', 'sani', 'xrdc', 'xpcs', 'dani', 'vean', 'cdih', 'sym', 'zres', 'bsa', 'deint', 'desolv'):
{===>} w_${w}_0=${weights[w][0]};
{===>} w_${w}_1=${weights[w][1]};
{===>} w_${w}_2=${weights[w][2]};
% endfor

{* It is possible to skip structures in the selection of structure in it0 *}
{* Give for this the number of structures to skip: *}
{===>} skip_struc=${skip_struc};

{======================= analysis and clustering ==========================}
{* Cutoff distance (proton-acceptor) to define an hydrogen bond? *}
{===>} dist_hb=${dist_hb};

{* Cutoff distance (carbon-carbon) to define an hydrophobic contact? *}
{===>} dist_nb=${dist_nb};

{* Clustering method (RMSD or Fraction of Common Contacts (FCC)) *}
{+ choice: "RMSD" "FCC" +}
{===>} clust_meth="${clust_meth}";

{* RMSD cutoff for clustering? (Recommended values: RMSD 7.5, FCC 0.75) *}
{===>} clust_cutoff=${clust_cutoff};

{* Minimum cluster size? *}
{===>} clust_size=${clust_size};

{* Chain-Agnostic Algorithm (used for FCC clustering in symmetrical complexes) *}
{+ choice: "true" "false" +}
{===>} fcc_ignc=${fcc_ignc};

{* Full or limited analysis of results? *}
{+ choice: "full" "cluster" "none" +}
{===>} runana=${runana};


{======================= final clean-up ===================================}
{* Clean up the run directory after completion (only files for struct #1 are kept) ? *}
{+ choice: true false +}
{===>} cleanup=true;

{============================ parallel jobs ===============================}
{* How many nodes do you want to use in parallel? *}
{* leave unused fields blank, make sure that the queues are actually running *}
{+ table: rows=10 "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"
 cols=3 "queue command" "cns executable" "number of jobs" +}
<% c = 1 %>
% for q in queues:
{===>} queue_${c}="${q['queue']}";
{===>} cns_exe_${c}="${q['cns_exe']}";
{===>} cpunumber_${c}=${q['cpunumber']};
    <% c += 1 %>
% endfor

{===========================================================================}
{        things below this line do not normally need to be changed          }
{===========================================================================}

) {- end block parameter definition -}

!for global parameters (local variables (suffix ) => global variables):
evaluate (&saprotocol.crossdock=&crossdock)
evaluate (&saprotocol.randorien=&randorien)
evaluate (&saprotocol.rebuildcplx=&rebuildcplx)
evaluate (&saprotocol.rigidmini=&rigidmini)
evaluate (&saprotocol.rigidtrans=&rigidtrans)
evaluate (&saprotocol.expand=&expand)
evaluate (&saprotocol.expansion=&expansion)
evaluate (&saprotocol.randangle=&randangle)

if (&saprotocol.expand eq true) then
  evaluate (&saprotocol.randorien=false)
  evaluate (&saprotocol.rigidmini=false)
end if
evaluate (&saprotocol.ntrials=&ntrials)
evaluate (&saprotocol.iniseed=&iniseed)
evaluate (&saprotocol.tadhigh_t=&tadhigh_t)
evaluate (&saprotocol.t1_init=&tadinit1_t)
evaluate (&saprotocol.t2_init=&tadinit2_t)
evaluate (&saprotocol.t3_init=&tadinit3_t)
evaluate (&saprotocol.t1_final=&tadfinal1_t)
evaluate (&saprotocol.t2_final=&tadfinal2_t)
evaluate (&saprotocol.t3_final=&tadfinal3_t)
evaluate (&saprotocol.inter_rigid=&inter_rigid)
evaluate (&saprotocol.inter_init_rigid=&init_rigid)
evaluate (&saprotocol.inter_fin_rigid=&fin_rigid)
evaluate (&saprotocol.inter_init_cool2=&init_cool2)
evaluate (&saprotocol.inter_fin_cool2=&fin_cool2)
evaluate (&saprotocol.inter_init_cool3=&init_cool3)
evaluate (&saprotocol.inter_fin_cool3=&fin_cool3)
evaluate (&saprotocol.rotate180_it0=&rotate180_it0)
evaluate (&saprotocol.rotate180_it1=&rotate180_it1)
evaluate (&saprotocol.tempstep=50)
evaluate (&saprotocol.timestep=&timestep)
evaluate (&saprotocol.tadfactor=&tadfactor)
evaluate (&saprotocol.emstepstrans=&emstepstrans)
evaluate (&saprotocol.initiosteps=&initiosteps)
evaluate (&saprotocol.cool1_steps=&cool1_steps)
evaluate (&saprotocol.cool2_steps=&cool2_steps)
evaluate (&saprotocol.cool3_steps=&cool3_steps)
evaluate (&saprotocol.fbeta=100)
evaluate (&saprotocol.mass=100)

evaluate (&filenames.fileroot=&fileroot)
evaluate (&filenames.template=&fileroot + "_1.pdb")

evaluate (&iterations.ini_count    =1)
evaluate (&iterations.structures   =&structures_$iteration)
evaluate (&iterations.keepstruct   =&keepstruct_$iteration)
evaluate (&iterations.w_vdw        =&w_vdw_$iteration)
evaluate (&iterations.w_elec       =&w_elec_$iteration)
evaluate (&iterations.w_dist       =&w_dist_$iteration)
evaluate (&iterations.w_rg         =&w_rg_$iteration)
evaluate (&iterations.w_sani       =&w_sani_$iteration)
evaluate (&iterations.w_xrdc       =&w_xrdc_$iteration)
evaluate (&iterations.w_xpcs       =&w_xpcs_$iteration)
evaluate (&iterations.w_dani       =&w_dani_$iteration)
evaluate (&iterations.w_vean       =&w_vean_$iteration)
evaluate (&iterations.w_cdih       =&w_cdih_$iteration)
evaluate (&iterations.w_sym        =&w_sym_$iteration)
evaluate (&iterations.w_zres       =&w_zres_$iteration)
evaluate (&iterations.w_bsa        =&w_bsa_$iteration)
evaluate (&iterations.w_deint      =&w_deint_$iteration)
evaluate (&iterations.w_desolv     =&w_desolv_$iteration)
evaluate (&iterations.anastruc     =&anastruc_$iteration)
evaluate (&iterations.w_lcc        = &w_lcc_$iteration)


evaluate (&data.ncomponents=&ncomponents)

evaluate ($nmol=1)
while ($nmol <= &data.ncomponents) loop mol

  !aa topology, linkage and parameters files
  evaluate (&toppar.prot_top_$nmol=&prot_top_mol$nmol )
  evaluate (&toppar.prot_link_$nmol=&prot_link_mol$nmol )
  evaluate (&toppar.prot_par_$nmol=&prot_par_mol$nmol )

  !coarse grained topology, linkage and parameters files
  evaluate (&toppar.prot_cg_top_$nmol=&prot_cg_top_mol$nmol )
  evaluate (&toppar.prot_cg_link_$nmol=&prot_cg_link_mol$nmol )
  evaluate (&toppar.prot_cg_par_$nmol=&prot_cg_par_mol$nmol )

  !molecule related (coordinate files, rootname, fix, type, coarse grained, segid)
  evaluate (&toppar.prot_coor_$nmol=&prot_coor_mol$nmol)
  evaluate (&toppar.prot_root_$nmol=&prot_root_mol$nmol)
  evaluate (&toppar.fix_origin_$nmol=&fix_origin_mol$nmol)
  evaluate (&toppar.dna_$nmol=&dna_mol$nmol)
  evaluate (&toppar.cyclicpept_$nmol=&cyclicpept_mol$nmol)
  evaluate (&toppar.shape_$nmol=&shape_mol$nmol)
  evaluate (&toppar.cg_$nmol=&cg_mol$nmol)
  evaluate (&toppar.prot_segid_$nmol=&prot_segid_mol$nmol)

  !semi flexible segments
  evaluate (&toppar.nseg_$nmol=&nseg_$nmol)
  evaluate ($nseg = 1)
  while ($nseg <= &toppar.nseg_$nmol) loop seg
    evaluate (&toppar.start_seg_$nmol_$nseg=&start_seg_$nmol_$nseg)
    evaluate (&toppar.end_seg_$nmol_$nseg=&end_seg_$nmol_$nseg)
    evaluate ($nseg = $nseg + 1)
  end loop seg

  !fully flexible segments
  evaluate (&toppar.nfle_$nmol=&nfle_$nmol)
  evaluate ($nfle = 1)
  while ($nfle <= &toppar.nfle_$nmol) loop fle
    evaluate (&toppar.start_fle_$nmol_$nfle=&start_fle_$nmol_$nfle)
    evaluate (&toppar.end_fle_$nmol_$nfle=&end_fle_$nmol_$nfle)
    evaluate ($nfle = $nfle + 1)
  end loop fle

  !histidine patches
  evaluate (&toppar.autohis=&autohis)
  evaluate (&toppar.nhisd_$nmol=&numhisd_$nmol)
  evaluate ($ncc=1)
  while ($ncc <= &toppar.nhisd_$nmol) loop hisd
    evaluate (&toppar.hisd_resid_$nmol_$ncc=&hisd_$nmol_$ncc)
    evaluate ($ncc = $ncc + 1)
  end loop hisd

  evaluate (&toppar.nhise_$nmol=&numhise_$nmol)
  evaluate ($ncc=1)
  while ($ncc <= &toppar.nhise_$nmol) loop hisd
    evaluate (&toppar.hise_resid_$nmol_$ncc=&hise_$nmol_$ncc)
    evaluate ($ncc = $ncc + 1)
  end loop hisd

  evaluate ($nmol = $nmol + 1)

end loop mol

! non-bonded parameter set to use
evaluate (&toppar.par_nonbonded=&par_nonbonded)

! z-restraining
evaluate (&Data.flags.zres =  &zres_on)
evaluate (&data.numzres=&numzres)
evaluate ($ncc=1)
while ($ncc <= &numzres) loop zres
  evaluate (&toppar.zres_sta_$ncc=&zres_sta_$ncc)
  evaluate (&toppar.zres_end_$ncc=&zres_end_$ncc)
  evaluate (&toppar.zres_seg_$ncc=&zres_seg_$ncc)
  evaluate (&toppar.zres_type_$ncc=&zres_type_$ncc)
  evaluate ($ncc = $ncc + 1)
end loop zres

! NCS restraints
evaluate (&data.kncs=&kncs)
evaluate (&Data.flags.ncs  =  &ncs_on)
evaluate (&data.numncs=&numncs)
evaluate ($ncc=1)
while ($ncc <= &numncs) loop ncs
  evaluate (&toppar.ncs_sta1_$ncc=&ncs_sta1_$ncc)
  evaluate (&toppar.ncs_end1_$ncc=&ncs_end1_$ncc)
  evaluate (&toppar.ncs_seg1_$ncc=&ncs_seg1_$ncc)
  evaluate (&toppar.ncs_sta2_$ncc=&ncs_sta2_$ncc)
  evaluate (&toppar.ncs_end2_$ncc=&ncs_end2_$ncc)
  evaluate (&toppar.ncs_seg2_$ncc=&ncs_seg2_$ncc)
  evaluate ($ncc = $ncc + 1)
end loop ncs

! Symmetry restraints
evaluate (&data.ksym=&ksym)
evaluate (&Data.flags.sym  =  &sym_on)
evaluate (&data.numc2sym=&numc2sym)
evaluate ($nsym=1)
while ($nsym <= &numc2sym) loop sym
  evaluate (&toppar.c2sym_sta1_$nsym=&c2sym_sta1_$nsym)
  evaluate (&toppar.c2sym_end1_$nsym=&c2sym_end1_$nsym)
  evaluate (&toppar.c2sym_seg1_$nsym=&c2sym_seg1_$nsym)
  evaluate (&toppar.c2sym_sta2_$nsym=&c2sym_sta2_$nsym)
  evaluate (&toppar.c2sym_end2_$nsym=&c2sym_end2_$nsym)
  evaluate (&toppar.c2sym_seg2_$nsym=&c2sym_seg2_$nsym)
  evaluate ($nsym = $nsym + 1)
end loop sym

evaluate (&data.numc3sym=&numc3sym)
evaluate ($nsym=1)
while ($nsym <= &numc3sym) loop sym
  evaluate (&toppar.c3sym_sta1_$nsym=&c3sym_sta1_$nsym)
  evaluate (&toppar.c3sym_end1_$nsym=&c3sym_end1_$nsym)
  evaluate (&toppar.c3sym_seg1_$nsym=&c3sym_seg1_$nsym)
  evaluate (&toppar.c3sym_sta2_$nsym=&c3sym_sta2_$nsym)
  evaluate (&toppar.c3sym_end2_$nsym=&c3sym_end2_$nsym)
  evaluate (&toppar.c3sym_seg2_$nsym=&c3sym_seg2_$nsym)
  evaluate (&toppar.c3sym_sta3_$nsym=&c3sym_sta3_$nsym)
  evaluate (&toppar.c3sym_end3_$nsym=&c3sym_end3_$nsym)
  evaluate (&toppar.c3sym_seg3_$nsym=&c3sym_seg3_$nsym)
  evaluate ($nsym = $nsym + 1)
end loop sym

evaluate (&data.nums3sym=&nums3sym)
evaluate ($nsym=1)
while ($nsym <= &nums3sym) loop sym
  evaluate (&toppar.s3sym_sta1_$nsym=&s3sym_sta1_$nsym)
  evaluate (&toppar.s3sym_end1_$nsym=&s3sym_end1_$nsym)
  evaluate (&toppar.s3sym_seg1_$nsym=&s3sym_seg1_$nsym)
  evaluate (&toppar.s3sym_sta2_$nsym=&s3sym_sta2_$nsym)
  evaluate (&toppar.s3sym_end2_$nsym=&s3sym_end2_$nsym)
  evaluate (&toppar.s3sym_seg2_$nsym=&s3sym_seg2_$nsym)
  evaluate (&toppar.s3sym_sta3_$nsym=&s3sym_sta3_$nsym)
  evaluate (&toppar.s3sym_end3_$nsym=&s3sym_end3_$nsym)
  evaluate (&toppar.s3sym_seg3_$nsym=&s3sym_seg3_$nsym)
  evaluate ($nsym = $nsym + 1)
end loop sym

evaluate (&data.numc4sym=&numc4sym)
evaluate ($nsym=1)
while ($nsym <= &numc4sym) loop sym
  evaluate (&toppar.c4sym_sta1_$nsym=&c4sym_sta1_$nsym)
  evaluate (&toppar.c4sym_end1_$nsym=&c4sym_end1_$nsym)
  evaluate (&toppar.c4sym_seg1_$nsym=&c4sym_seg1_$nsym)
  evaluate (&toppar.c4sym_sta2_$nsym=&c4sym_sta2_$nsym)
  evaluate (&toppar.c4sym_end2_$nsym=&c4sym_end2_$nsym)
  evaluate (&toppar.c4sym_seg2_$nsym=&c4sym_seg2_$nsym)
  evaluate (&toppar.c4sym_sta3_$nsym=&c4sym_sta3_$nsym)
  evaluate (&toppar.c4sym_end3_$nsym=&c4sym_end3_$nsym)
  evaluate (&toppar.c4sym_seg3_$nsym=&c4sym_seg3_$nsym)
  evaluate (&toppar.c4sym_sta4_$nsym=&c4sym_sta4_$nsym)
  evaluate (&toppar.c4sym_end4_$nsym=&c4sym_end4_$nsym)
  evaluate (&toppar.c4sym_seg4_$nsym=&c4sym_seg4_$nsym)
  evaluate ($nsym = $nsym + 1)
end loop sym

evaluate (&data.numc5sym=&numc5sym)
evaluate ($nsym=1)
while ($nsym <= &numc5sym) loop sym
  evaluate (&toppar.c5sym_sta1_$nsym=&c5sym_sta1_$nsym)
  evaluate (&toppar.c5sym_end1_$nsym=&c5sym_end1_$nsym)
  evaluate (&toppar.c5sym_seg1_$nsym=&c5sym_seg1_$nsym)
  evaluate (&toppar.c5sym_sta2_$nsym=&c5sym_sta2_$nsym)
  evaluate (&toppar.c5sym_end2_$nsym=&c5sym_end2_$nsym)
  evaluate (&toppar.c5sym_seg2_$nsym=&c5sym_seg2_$nsym)
  evaluate (&toppar.c5sym_sta3_$nsym=&c5sym_sta3_$nsym)
  evaluate (&toppar.c5sym_end3_$nsym=&c5sym_end3_$nsym)
  evaluate (&toppar.c5sym_seg3_$nsym=&c5sym_seg3_$nsym)
  evaluate (&toppar.c5sym_sta4_$nsym=&c5sym_sta4_$nsym)
  evaluate (&toppar.c5sym_end4_$nsym=&c5sym_end4_$nsym)
  evaluate (&toppar.c5sym_seg4_$nsym=&c5sym_seg4_$nsym)
  evaluate (&toppar.c5sym_sta5_$nsym=&c5sym_sta5_$nsym)
  evaluate (&toppar.c5sym_end5_$nsym=&c5sym_end5_$nsym)
  evaluate (&toppar.c5sym_seg5_$nsym=&c5sym_seg5_$nsym)
  evaluate ($nsym = $nsym + 1)
end loop sym

evaluate (&data.numc6sym=&numc6sym)
evaluate ($nsym=1)
while ($nsym <= &numc6sym) loop sym
  evaluate (&toppar.c6sym_sta1_$nsym=&c6sym_sta1_$nsym)
  evaluate (&toppar.c6sym_end1_$nsym=&c6sym_end1_$nsym)
  evaluate (&toppar.c6sym_seg1_$nsym=&c6sym_seg1_$nsym)
  evaluate (&toppar.c6sym_sta2_$nsym=&c6sym_sta2_$nsym)
  evaluate (&toppar.c6sym_end2_$nsym=&c6sym_end2_$nsym)
  evaluate (&toppar.c6sym_seg2_$nsym=&c6sym_seg2_$nsym)
  evaluate (&toppar.c6sym_sta3_$nsym=&c6sym_sta3_$nsym)
  evaluate (&toppar.c6sym_end3_$nsym=&c6sym_end3_$nsym)
  evaluate (&toppar.c6sym_seg3_$nsym=&c6sym_seg3_$nsym)
  evaluate (&toppar.c6sym_sta4_$nsym=&c6sym_sta4_$nsym)
  evaluate (&toppar.c6sym_end4_$nsym=&c6sym_end4_$nsym)
  evaluate (&toppar.c6sym_seg4_$nsym=&c6sym_seg4_$nsym)
  evaluate (&toppar.c6sym_sta5_$nsym=&c6sym_sta5_$nsym)
  evaluate (&toppar.c6sym_end5_$nsym=&c6sym_end5_$nsym)
  evaluate (&toppar.c6sym_seg5_$nsym=&c6sym_seg5_$nsym)
  evaluate (&toppar.c6sym_sta6_$nsym=&c6sym_sta6_$nsym)
  evaluate (&toppar.c6sym_end6_$nsym=&c6sym_end6_$nsym)
  evaluate (&toppar.c6sym_seg6_$nsym=&c6sym_seg6_$nsym)
  evaluate ($nsym = $nsym + 1)
end loop sym

if ( &data.numc2sym eq 6) then
  evaluate (&saprotocol.rotate180_it0 = false)
  evaluate (&saprotocol.rotate180_it1 = false)
end if
if ( &data.numc3sym ne 0) then
  evaluate (&saprotocol.rotate180_it0 = false)
  evaluate (&saprotocol.rotate180_it1 = false)
end if
if ( &data.numc4sym ne 0) then
  evaluate (&saprotocol.rotate180_it0 = false)
  evaluate (&saprotocol.rotate180_it1 = false)
end if
if ( &data.numc5sym ne 0) then
  evaluate (&saprotocol.rotate180_it0 = false)
  evaluate (&saprotocol.rotate180_it1 = false)
end if
if ( &data.numc6sym ne 0) then
  evaluate (&saprotocol.rotate180_it0 = false)
  evaluate (&saprotocol.rotate180_it1 = false)
end if


!Dihedrals, DNA and distance restraints
evaluate (&Data.dnarest = &dnarest_on)
evaluate (&Data.flags.cdih =  &dihedrals_on)
evaluate (&Data.cdih.on = &dihedrals_on)
evaluate (&Data.ssdihed = &ssdihed)
evaluate (&Data.error_dih = &error_dih)
evaluate (&data.dihedrals.on=&dihedrals_on)
evaluate (&data.dihedrals_hot=&dihedrals_hot)
evaluate (&data.dihedrals_cool1=&dihedrals_cool1)
evaluate (&data.dihedrals_cool2=&dihedrals_cool2)
evaluate (&data.dihedrals_cool3=&dihedrals_cool3)
evaluate (&data.hbonds_on=&hbonds_on)

! RDC restraints
evaluate (&Data.flags.vean =  false)
evaluate (&Data.flags.xrdc =  false)
evaluate (&Data.flags.sani =  false)
evaluate (&data.numrdc=&numrdc)
evaluate ($ncc=1)
while ($ncc <= &data.numrdc) loop rdc
  if (&rdc_choice_$ncc = "VANGLE") then
    evaluate (&Data.flags.vean =  true)
  end if
  if (&rdc_choice_$ncc = "SANI") then
    evaluate (&Data.flags.sani =  true)
  end if
  if (&rdc_choice_$ncc = "XRDC") then
    evaluate (&Data.flags.xrdc =  true)
  end if
  evaluate (&data.rdc_choice_$ncc=&rdc_choice_$ncc)
  evaluate (&data.rdc_firstIt_$ncc=&rdc_firstIt_$ncc)
  evaluate (&data.rdc_lastIt_$ncc=&rdc_lastIt_$ncc)
  evaluate (&data.rdc_hot_$ncc=&rdc_hot_$ncc)
  evaluate (&data.rdc_cool1_$ncc=&rdc_cool1_$ncc)
  evaluate (&data.rdc_cool2_$ncc=&rdc_cool2_$ncc)
  evaluate (&data.rdc_cool3_$ncc=&rdc_cool3_$ncc)
  evaluate (&data.rdc_r_$ncc=&rdc_r_$ncc)
  evaluate (&data.rdc_d_$ncc=&rdc_d_$ncc)
  evaluate (&data.ini_bor_hot_$ncc=&ini_bor_hot_$ncc)
  evaluate (&data.ini_bor_cool1_$ncc=&ini_bor_cool1_$ncc)
  evaluate (&data.ini_bor_cool2_$ncc=&ini_bor_cool2_$ncc)
  evaluate (&data.ini_bor_cool3_$ncc=&ini_bor_cool3_$ncc)
  evaluate (&data.ini_cen_hot_$ncc=&ini_cen_hot_$ncc)
  evaluate (&data.ini_cen_cool1_$ncc=&ini_cen_cool1_$ncc)
  evaluate (&data.ini_cen_cool2_$ncc=&ini_cen_cool2_$ncc)
  evaluate (&data.ini_cen_cool3_$ncc=&ini_cen_cool3_$ncc)
  evaluate (&data.fin_bor_hot_$ncc=&fin_bor_hot_$ncc)
  evaluate (&data.fin_bor_cool1_$ncc=&fin_bor_cool1_$ncc)
  evaluate (&data.fin_bor_cool2_$ncc=&fin_bor_cool2_$ncc)
  evaluate (&data.fin_bor_cool3_$ncc=&fin_bor_cool3_$ncc)
  evaluate (&data.fin_cen_hot_$ncc=&fin_cen_hot_$ncc)
  evaluate (&data.fin_cen_cool1_$ncc=&fin_cen_cool1_$ncc)
  evaluate (&data.fin_cen_cool2_$ncc=&fin_cen_cool2_$ncc)
  evaluate (&data.fin_cen_cool3_$ncc=&fin_cen_cool3_$ncc)
  evaluate ($ncc=$ncc+1)
end loop rdc

! PCS restraints
evaluate (&Data.flags.xpcs =  false)
evaluate (&data.numpcs=&numpcs)
evaluate ($ncc=1)
while ($ncc <= &numpcs) loop pcs
  if (&pcs_choice_$ncc = "XPCS") then
    evaluate (&Data.flags.xpcs =  true)
  end if
  evaluate (&data.pcs_choice_$ncc=&pcs_choice_$ncc)
  evaluate (&data.pcs_firstIt_$ncc=&pcs_firstIt_$ncc)
  evaluate (&data.pcs_lastIt_$ncc=&pcs_lastIt_$ncc)
  evaluate (&data.pcs_hot_$ncc=&pcs_hot_$ncc)
  evaluate (&data.pcs_cool1_$ncc=&pcs_cool1_$ncc)
  evaluate (&data.pcs_cool2_$ncc=&pcs_cool2_$ncc)
  evaluate (&data.pcs_cool3_$ncc=&pcs_cool3_$ncc)
  evaluate (&data.pcs_r_$ncc=&pcs_r_$ncc)
  evaluate (&data.pcs_d_$ncc=&pcs_d_$ncc)
  evaluate ($ncc=$ncc+1)
end loop pcs

! DANI restraints
evaluate (&Data.flags.dani =  false)
evaluate (&data.numdani=&numdani)
evaluate ($ncc=1)
while ($ncc <= &numdani) loop dani
  if (&dan_choice_$ncc = "DANI") then
    evaluate (&Data.flags.dani =  true)
  end if
  evaluate (&data.dan_choice_$ncc=&dan_choice_$ncc)
  evaluate (&data.dan_firstIt_$ncc=&dan_firstIt_$ncc)
  evaluate (&data.dan_lastIt_$ncc=&dan_lastIt_$ncc)
  evaluate (&data.dan_hot_$ncc=&dan_hot_$ncc)
  evaluate (&data.dan_cool1_$ncc=&dan_cool1_$ncc)
  evaluate (&data.dan_cool2_$ncc=&dan_cool2_$ncc)
  evaluate (&data.dan_cool3_$ncc=&dan_cool3_$ncc)
  evaluate (&data.dan_tc_$ncc=&dan_tc_$ncc)
  evaluate (&data.dan_anis_$ncc=&dan_anis_$ncc)
  evaluate (&data.dan_r_$ncc=&dan_r_$ncc)
  evaluate (&data.dan_wh_$ncc=&dan_wh_$ncc)
  evaluate (&data.dan_wn_$ncc=&dan_wn_$ncc)
  evaluate ($ncc=$ncc+1)
end loop dani

! planarity restraints
evaluate (&Data.flags.plan =  false)

! distance restraints
evaluate (&Data.flags.noe  =  true)
evaluate (&data.scaling=&air_scaling)
evaluate (&data.totnoe_unamb=&tot_unamb)
evaluate (&data.unamb_firstit=&unamb_firstit)
evaluate (&data.unamb_lastit=&unamb_lastit)
evaluate (&data.unamb_hot=&unamb_hot)
evaluate (&data.unamb_cool1=&unamb_cool1)
evaluate (&data.unamb_cool2=&unamb_cool2)
evaluate (&data.unamb_cool3=&unamb_cool3)
evaluate (&data.noecv=&noecv)
evaluate (&data.ncvpart=&ncvpart)

evaluate (&data.totnoe_amb=&tot_amb)
evaluate (&data.amb_firstit=&amb_firstit)
evaluate (&data.amb_lastit=&amb_lastit)
evaluate (&data.amb_hot=&amb_hot)
evaluate (&data.amb_cool1=&amb_cool1)
evaluate (&data.amb_cool2=&amb_cool2)
evaluate (&data.amb_cool3=&amb_cool3)

evaluate (&data.hbond_firstit=&hbond_firstit)
evaluate (&data.hbond_lastit=&hbond_lastit)
evaluate (&data.hbond_hot=&hbond_hot)
evaluate (&data.hbond_cool1=&hbond_cool1)
evaluate (&data.hbond_cool2=&hbond_cool2)
evaluate (&data.hbond_cool3=&hbond_cool3)

evaluate (&data.mrswi_hot=&mrswi_hot)
evaluate (&data.mrswi_cool1=&mrswi_cool1)
evaluate (&data.mrswi_cool2=&mrswi_cool2)
evaluate (&data.mrswi_cool3=&mrswi_cool3)

evaluate (&data.rswi_hot=&rswi_hot)
evaluate (&data.rswi_cool1=&rswi_cool1)
evaluate (&data.rswi_cool2=&rswi_cool2)
evaluate (&data.rswi_cool3=&rswi_cool3)

evaluate (&data.masy_hot=&masy_hot)
evaluate (&data.masy_cool1=&masy_cool1)
evaluate (&data.masy_cool2=&masy_cool2)
evaluate (&data.masy_cool3=&masy_cool3)

evaluate (&data.asy_hot=&asy_hot)
evaluate (&data.asy_cool1=&asy_cool1)
evaluate (&data.asy_cool2=&asy_cool2)
evaluate (&data.asy_cool3=&asy_cool3)

evaluate (&data.ranair=&ranair)
if (&data.ranair eq true) then
  evaluate (&data.noecv = false)
end if
evaluate (&data.cmrest=&cmrest)
evaluate (&data.cmtight=&cmtight)
evaluate (&data.kcont=&kcont)
evaluate (&data.surfrest=&surfrest)
evaluate (&data.ksurf=&ksurf)


! radius of gydration restraints
evaluate (&data.flags.rg=&rgrest)
evaluate (&data.rgtarg=&rgtarg)
evaluate (&data.krg_hot=&krg_hot)
evaluate (&data.krg_cool1=&krg_cool1)
evaluate (&data.krg_cool2=&krg_cool2)
evaluate (&data.krg_cool3=&krg_cool3)
evaluate (&data.rgsele=&rgsele)

evaluate (&data.kzres=&kzres)
evaluate (&data.zresmax=&zresmax)
evaluate (&data.zresmin=&zresmin)

! keep or delete non-polar hydrogens
evaluate (&toppar.delenph=&delenph)


!Electrostatics:
evaluate (&Data.flags.dihed = &dihedflag)
evaluate (&Data.flags.elec0 = &elecflag_0)
evaluate (&Data.flags.elec1 = &elecflag_1)
evaluate (&Data.epsilon0 = &epsilon_0)
evaluate (&Data.epsilon1 = &epsilon_1)
evaluate (&Data.dielec0  = &dielec_0)
evaluate (&Data.dielec1  = &dielec_1)

!Check for CG and if present force cdie
evaluate($cg = false)
evaluate($nchain1 = 0)
while ($nchain1 < $data.ncomponents) loop cgmol
  evaluate($nchain1 = $nchain1 + 1)
  if ($Toppar.cg_$nchain1 eq true) then
    evaluate($cg = true)
  end if
end loop cgmol

if ($cg = true) then
  evaluate (&Data.dielec0  = cdie)
  evaluate (&Data.dielec1  = cdie)
  evaluate (&Data.epsilon1 = &Data.epsilon0)
  display "FORCING CDIE FOR ELECTROSTATICS BECAUSE OF COARSE GRAINING"
end if

!Interaction matrix:
evaluate ($nmol1=1)
while ($nmol1 <= &data.ncomponents) loop mol1
  evaluate ($nmol2=$nmol1 + 1)
  evaluate (&toppar.int_$nmol1_$nmol1 = &int_$nmol1_$nmol1)
  while ($nmol2 <= &data.ncomponents) loop mol2
    evaluate (&toppar.int_$nmol1_$nmol2 = &int_$nmol1_$nmol2)
    evaluate (&toppar.int_$nmol2_$nmol1 = &int_$nmol1_$nmol2)
    evaluate ($nmol2=$nmol2 + 1)
  end loop mol2
  evaluate ($nmol1 = $nmol1 + 1)
end loop mol1

!intermolecular contacts analysis
evaluate (&data.hb_dist=&dist_hb)
evaluate (&data.nb_dist=&dist_nb)


!water refinement
evaluate (&refine.firstwater=&firstwater)
evaluate (&refine.solvshell=&solvshell)
evaluate (&refine.keepwater=&keepwater)
evaluate (&refine.waterrefine=min(&structures_1,&waterrefine))
evaluate (&refine.solvent=&solvent)
evaluate (&refine.heatsteps=&waterheatsteps)
evaluate (&refine.steps=&watersteps)
evaluate (&refine.coolsteps=&watercoolsteps)


!for the non-bonded parameters (the section was taken out of
!parallhdg5.0.pro and parallhdg5.1.pro, so be careful!):
if (&toppar.par_nonbonded eq "PROLSQ") then
    evaluate (&toppar.repel_radius = 1.0)
    evaluate (&toppar.repel_rcons = 20)
    evaluate (&toppar.repel_rexpo  = 4)
    evaluate (&toppar.repel_irexp  = 1)
elseif (&toppar.par_nonbonded eq "PARMALLH6") then
    evaluate (&toppar.repel_radius = 0.8)
    evaluate (&toppar.repel_rcons = 5.0)
    evaluate (&toppar.repel_rexpo  = 2)
    evaluate (&toppar.repel_irexp  = 2)
elseif (&toppar.par_nonbonded eq "OPLSX") then
    evaluate (&toppar.repel_radius = 0.0)
else        {...now the standard PARALLHDG parameters}
    evaluate (&toppar.repel_radius = 0.78)
    evaluate (&toppar.repel_rcons = 5.0)
    evaluate (&toppar.repel_rexpo  = 2)
    evaluate (&toppar.repel_irexp  = 2)
end if

! Water in rigid body docking
evaluate (&data.waterdock=&waterdock)
evaluate (&data.db_method=&db_method)
evaluate (&data.water_tokeep=&water_tokeep)
evaluate (&data.dnap_water_tokeep=&dnap_water_tokeep)
evaluate (&data.water_randfrac=&water_randfrac)
evaluate (&data.solvate_method=&solvate_method)
evaluate (&data.water_surfcutoff=&water_surfcutoff)
evaluate (&data.water_analysis=&water_analysis)
evaluate (&data.transwater=&transwater)
evaluate (&data.water_restraint_initial=&water_restraint_initial)
evaluate (&data.water_restraint_cutoff=&water_restraint_cutoff)
evaluate (&data.water_restraint_scale=&water_restraint_scale)
evaluate (&data.waterensemble=&waterensemble)

! Centroid parameters
eval(&data.flags.centroids = &centroid_rest)
eval(&data.centroids.kscale = &centroid_kscale)
eval($nchain = 0)
while ($nchain < &ncomponents) loop nloop1
    eval($nchain = $nchain + 1)
    eval(&data.centroids.xcom_$nchain = &xcom_$nchain)
    eval(&data.centroids.ycom_$nchain = &ycom_$nchain)
    eval(&data.centroids.zcom_$nchain = &zcom_$nchain)
    eval(&data.centroids.ambi_$nchain = &ambi_$nchain)
end loop nloop1

if (&saprotocol.expand eq true) then
  eval(&data.flags.centroids = true)
end if

! Cryo-EM parameters
eval(&data.flags.em = &em_rest)
eval(&data.em.kscale = &em_kscale)
eval(&data.em.it0 = &em_it0)
eval(&data.em.it1 = &em_it1)
eval(&data.em.itw = &em_itw)
eval(&data.em.resolution = &em_resolution)
eval(&data.em.nx = &nx)
eval(&data.em.ny = &ny)
eval(&data.em.nz = &nz)
eval(&data.em.xlength = &xlength)
eval(&data.em.ylength = &ylength)
eval(&data.em.zlength = &zlength)

if (&data.flags.em eq true) then
  eval(&data.waterdock = false)
  display "EM restraints and solvated docking are incompatible - turning solvated docking OFF"
end if

if (&data.waterdock eq true) then
  evaluate (&saprotocol.rotate180_it0 = false)
  evaluate (&saprotocol.rotate180_it1 = false)
  evaluate (&SaProtocol.initiosteps = 0)
  evaluate (&SaProtocol.cool1_steps = 0)
  evaluate (&refine.keepwater = true)
  display SOLVATED DOCKING TURNED ON: initiosteps and cool1_steps set to 0, rotate180 set to false
end if

! Flexible region cutoff
evaluate (&data.flcut_nb = &flcut_nb)
