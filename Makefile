SHELL=/bin/bash
PDFLATEX=pdflatex -halt-on-error -interaction nonstopmode -file-line-error

default : enkidu.pdf
.PHONY : default

figures=intro vectorops vectorproj intersect lineparab
figures_py=$(foreach f,$(figures),$(f).py)
figures_src=$(foreach f,$(figures),$(f)-src.tex)
figures_tex=$(foreach f,$(figures),$(f)-fig.tex)
figures_eps=$(foreach f,$(figures),$(f)-fig.eps)
figures_pdf=$(foreach f,$(figures),$(f)-fig.pdf)

tarrable=$(figures_py) enkidu.py enkidu.tex enkidudoc.cls Makefile meta/author meta/date meta/url meta/version py2src enkidu.pdf README
tardir:=enkidu-$(shell cat meta/version)

tar : $(tardir).tar.bz2
.PHONY : tar

$(tardir).tar.bz2 $(tardir).zip : $(tarrable)
	mkdir $(tardir)
	tar cf - $(tarrable) | tar xf - -C $(tardir)
	tar jcf $@ $(tardir)
	zip -r $(tardir).zip $(tardir)
	rm -rf $(tardir)

enkidu.pdf : enkidu.tex meta.tex enkidudoc.cls $(figures_pdf) $(figures_src)
	rm -f enkidu.toc
	$(PDFLATEX) $<
	while egrep '(Rerun|No file enkidu.(aux|toc))' enkidu.log >/dev/null; \
	    do $(PDFLATEX) $< ; done

meta.tex : meta/title meta/author meta/date meta/url
	(echo "\title{`cat meta/title`}"; \
	echo "\pdfinfo{/Title (`cat meta/title`)}"; \
	echo "\author{`cat meta/author`}"; \
	echo "\pdfinfo{/Author (`cat meta/author`)}"; \
	echo "\date{`date -d $$(cat meta/date) '+%Y %B %-d'`}"; \
	echo "\websiteurl{`cat meta/url`}"; \
	echo "\websitetext{`sed 's,~,\\\\textasciitilde ,g' meta/url`}") \
	>$@
meta/title : meta/version
	echo "Documentation for Enkidu v"`cat meta/version` >$@

%-src.tex : %.py
	awk -f py2src $< >$@

# epstopdf (from teTeX) gets the bounding box right, but if allowed
# to run ghostscript itself, gets the status code wrong; ps2pdf (from
# ghostscript) gets the status code right but the bounding box wrong.
%-fig.pdf : %-fig.eps
	epstopdf --nogs $< |ps2pdf - $@
%-fig.eps %-fig.tex : %.py enkidu.py
	python $< $*-fig pdf

clean :
	rm -f enkidu.{aux,log,out,pdf,toc} meta.tex $(figures_tex) $(figures_eps) $(figures_pdf) $(figures_src) enkidu.pyc meta/title
.PHONY : clean

