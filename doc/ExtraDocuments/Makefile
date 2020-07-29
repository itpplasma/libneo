all: coding_style

coding_style: coding_style.pdf

coding_style.pdf: clean coding_style.tex
	pdflatex coding_style
#~ 	makeindex coding_style
	pdflatex coding_style
	latex_count=8 ; \
	while egrep -s 'Rerun (LaTeX|to get cross-references right)' coding_style.log && [ $$latex_count -gt 0 ] ;\
	    do \
	      echo "Rerunning latex...." ;\
	      pdflatex coding_style ;\
	      latex_count=`expr $$latex_count - 1` ;\
	    done
#~ 	makeindex main.idx
	pdflatex coding_style

clean:
	rm -f *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *.brf *.blg *.bbl *.pdf
