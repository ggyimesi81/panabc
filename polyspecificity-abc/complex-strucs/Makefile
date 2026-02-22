ALL_SYSTEMS := $(patsubst %.complex.cif,%,$(wildcard *.complex.cif))

.PHONY:	all all-receptor

all:	all-receptor

all-receptor:	$(addsuffix .receptor.pdb,$(ALL_SYSTEMS))

%.receptor.pdb:	%.complex.pdb
	grep -v '^HETATM' $< >$@

