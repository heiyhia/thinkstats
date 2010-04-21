LATEX = latex

DVIPS = dvips

PDFFLAGS = -dCompatibilityLevel=1.4 -dPDFSETTINGS=/prepress \
           -dCompressPages=true -dUseFlateCompression=true  \
           -dEmbedAllFonts=true -dSubsetFonts=true -dMaxSubsetPct=100


%.dvi: %.tex
	$(LATEX) $<

%.ps: %.dvi
	$(DVIPS) -o $@ $<

%.pdf: %.ps
	ps2pdf $(PDFFLAGS) $<

all:	book.tex
	latex book
	makeindex book
	latex book
	dvips -t letter -Ppdf -o thinkstats.ps book
#	dvips -T 6.75in,9.25in -Ppdf -o thinkstats.ps book
#	dvips -t b5 -Ppdf -o thinkstats.ps book
#	dvips -T 7in,10in -Ppdf -o thinkstats.ps book
	gv thinkstats.ps

html:	book.tex header.html footer.html
	rm -rf html
	mkdir html
	hevea -O -e latexonly htmlonly book
# the following line is a kludge to prevent imagen from seeing
# the definitions in latexonly
	#grep -v latexonly book.image.tex > a; mv a book.image.tex
	#imagen -png book
	hacha book.html
	cp up.png next.png back.png html
	mv index.html book.css book*.html *motif.gif html

DEST = /home/downey/public_html/greent/compprobstat

distrib:
	ps2pdf $(PDFFLAGS) thinkstats.ps
	rm -rf dist
	mkdir dist dist/tex dist/tex/figs
	rsync -a thinkstats.pdf thinkstats.ps html cover_nolines.png dist
	rsync -a Makefile book.tex latexonly htmlonly dist/tex
	# rsync -a figs/*.fig figs/*.eps dist/tex/figs
	cd dist; zip -r thinkstats.tex.zip tex
	cd dist; zip -r thinkstats.html.zip html
	rsync -a dist/* $(DEST)
	chmod -R o+r $(DEST)/*

clean:
	rm -f *~ *.aux *.log *.dvi *.idx *.ilg *.ind *.toc



