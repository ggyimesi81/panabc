SHELL := /bin/bash
MSMS := /opt/msms/msms.x86_64Linux2.2.6.1.staticgcc
PDB_TO_XYZR := /opt/msms/pdb_to_xyzr
MGLTOOLS := /opt/mgltools_x86_64Linux2_1.5.7

ALL_SYSTEMS := $(patsubst ../%.receptor.pdb,%,$(wildcard ../*.receptor.pdb))

# these have either zero length files or too large structure
BLACKLISTED := 5umd_YMZ_TA 7r2g_MIX_H 7r2g_MIX_I 7r2g_MIX_J 7r2g_MIX_K 7r2g_MIX_L 7r2g_MIX_M 7r2g_MIX_N 7r2g_MIX_O 7r2g_MIX_P 7r2g_MIX_Q 8ogj_YMZ_BC 9eok_TA1_RD 9g1z_YMZ_DC
# these are CA-only structures
BLACKLISTED += 4aqv_TA1_G 4aqw_TA1_G 4ck5_TA1_F 4ck6_TA1_G 4ck7_TA1_G 4lz9_RHQ_B 4uxo_TA1_H 4uxp_TA1_H 4uxr_TA1_H 4uxs_TA1_H 4uxt_TA1_H 4uxy_TA1_H 4uy0_TA1_H 5m5i_TA1_G 5m5l_TA1_G 5m5m_TA1_G 5m5n_TA1_G 5m5o_TA1_G

ALL_SYSTEMS := $(filter-out $(BLACKLISTED),$(ALL_SYSTEMS))

all:	all-areas all-receptor-pdbqt all-surfarea2 all-ligand-areas

all-areas:	$(addsuffix .area,$(ALL_SYSTEMS))

.PHONY:	all all-areas all-receptor-pdbqt all-surfarea2 all-ligand-areas

.SECONDARY:

.DELETE_ON_ERROR:

%.receptor.H.pdb: ../%.receptor.pdb
	-reduce -NOFLIP $< >$@ || echo "reduce error!"

%.receptor.xyzr:	../%.receptor.pdb
	$(PDB_TO_XYZR) $< >$@

all-ligand-areas:	$(addsuffix _ligand.area,$(ALL_SYSTEMS))

%.ligand.pdb:	../%.complex.pdb
	python extract-ligand.py --renumber $* $@

%.ligand.xyzr:	%.ligand.pdb
	$(PDB_TO_XYZR) $< >$@

%_ligand.area %.ligand.surface.vert %.ligand.surface.face:	%.ligand.xyzr
	$(MSMS) -if $< -of $*.ligand.surface -af $*_ligand -probe_radius 1.4

%.area %.surface.vert %.surface.face:	%.receptor.xyzr
	$(MSMS) -if $< -of $*.surface -af $* -probe_radius 1.4

all-receptor-pdbqt:	$(addsuffix .receptor.pdbqt,$(ALL_SYSTEMS))

# PDBQT needs (polar) hydrogens because they will become the H-bond donors (atom type "HD")
%.receptor.pdbqt:	%.receptor.H.pdb
	#obabel -ipdb $< -opdbqt -O $@ -xr -xh --partialcharge gasteiger
	LD_LIBRARY_PATH=$(MGLTOOLS)/lib:"$$LD_LIBRARY_PATH" $(MGLTOOLS)/bin/python2.7 $(MGLTOOLS)/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r $< -o $@ -U nphs

all-surfarea2:	$(addsuffix .surfarea2.tsv,$(ALL_SYSTEMS))

%.surfarea2.tsv:	%.area %.receptor.pdbqt ../%.receptor.pdb
	python surfarea2.py $* >$@
	
