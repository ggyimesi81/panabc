ALL_PDBIDS := $(patsubst mmcif/%.cif.gz,%,$(wildcard mmcif/*/*.cif.gz))

.PHONY:	all
all:	all-mmjsons

.PHONY:	all-mmjsons
all-mmjsons:	$(foreach pdb_id,$(ALL_PDBIDS),mmjson/$(pdb_id).json.gz)

ALL_ASSEMBLIES := $(patsubst biounit-mmcif/%.cif.gz,%,$(wildcard biounit-mmcif/*/*.cif.gz))

.PHONY:	all-biounit-mmjsons
all-biounit-mmjsons:	$(foreach assy_id,$(ALL_ASSEMBLIES),biounit-mmjson/$(assy_id).json.gz)

# gemmi sometimes generates floats like -.1234, which python json will refuse to read
# running the file through jq corrects this, we also make it more compact
mmjson/%.json.gz:	mmcif/%.cif.gz
	mkdir -p $$(dirname $@)
	gemmi cif2json -m --bare-tags $< - | jq -c . | gzip -9 >$@

biounit-mmjson/%.json.gz:	biounit-mmcif/%.cif.gz
	mkdir -p $$(dirname $@)
	gemmi cif2json -m --bare-tags $< - | jq -c . | gzip -9 >$@

