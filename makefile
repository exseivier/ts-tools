binary="/Users/javier/binarios/"
source="/Users/javier/tScan/"
PYTHON := $(shell which python 2>/dev/null)
PERL := $(shell which perl 2>/dev/null)
R := $(shell which R 2>/dev/null)


all: installation install_info clear

check_python:
ifdef PYTHON
	@echo "Python was found"
	@touch "python.exist"
else
	@echo "Python was not found"
	@exit 1
endif

check_perl:
ifdef PERL
	@echo "Perl was found"
	@touch "perl.exist"
else
	@echo "Perl was not found"
endif

check_R:
ifdef R
	@echo "R was found"
	@touch "R.exist"
else
	@echo "R was not found"
endif

install_info:
	@echo "Installing in $(binary) and source files are in $(source)" 
	[[ -d $(binary) ]] || (echo "$(binary) was not found" && mkdir $(binary))
	[[ -d $(source) ]] || (echo "$(source) was not found" && mkdir $(source))

installation: check_python python.exist check_perl perl.exist check_R R.exist scripts/*.py scripts/*.sh scripts/*.R scripts/*.pl
	@echo "Procede to install..."
	cp scripts/*.py $(source)
	cp scripts/*.sh $(source)
	cp scripts/*.pl $(source)
	cp scripts/*.R $(source)
	cp ts-tools* $(source)
	ln -s $(source)/fa2ints.py $(binary)fa2ints
	ln -s $(source)/miRNAfamily-build.py $(binary)/miRNAfamily-build
	ln -s $(source)/targetscan_70.pl $(binary)/targetScan
	ln -s $(source)/target2context.pl $(binary)/targetScan-ctx++
	ln -s $(source)/separateTSresults-2.py $(binary)/targetScan-final
	ln -s $(source)/targetScan-pp.sh $(binary)/targetScan-pp
	ln -s $(source)/targetScan-stats-2.R $(binary)/targetScan-stats

clear:
	@rm python.exist perl.exist R.exist
