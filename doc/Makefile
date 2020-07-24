
all: source UserDoc

# Generate doxygen documentation.
source:
	@echo "Start extracting comments from source files and making the html/pdf output."
	@echo " Output of doxygen/latex is written to corresponding log/err file."
	@echo "---------------------------------------------------------------------------"
	doxygen 1> doxygen.log 2> doxygen.err
	make -C latex/ 1> latex.log 2> latex.err
	cp latex/refman.pdf ./
	make -C latex/ clean
	@echo "---------------------------------------------------------------------------"
	@echo "...done."

UserDoc:
	make -C UserDoc
	cp UserDoc/main.pdf ./UserDoc.pdf

# Clean up.
clean:
	rm -f {.,etc,chease}/*.pdf
	rm -f manual/{doc,doc-diff}.pdf
	rm -f {.,manual,etc,chease}/*.{ps,log,aux,nlo,toc,dvi,out}
	rm -rf latex html

########################################################################
### Low-level targets, should usually not be called directly.

%.dvi : %.tex
	@echo make $@ with $< in $(<D) ;\
	cd $(<D); \
	$(LATEX) $(<F) ; \
	while [ -n "$$(grep $(RERUN_LATEX_LINE) *.log)" ]; do \
	     $(LATEX) $(<F); \
	done

%.ps : %.dvi
	@echo make $@ with $< in $(<D) ;\
	cd $(<D) ; pwd ;\
	$(DVIPS) $(<F) -o $(@F)

%.pdf : %.ps
	@echo make $@ with $< in $(<D);\
	cd $(<D) ; pwd ;\
	$(PS2PDF) $(<F) $(@F)

.PHONY: all UserDoc clean
