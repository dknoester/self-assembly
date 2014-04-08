%.pdf: %.eps; ps2pdf14 -dPDFSETTINGS=/prepress -dEPSCrop $< $@
%.eps: %.dot; neato -Tps2 -o $@ $<
%.pdftex: %.fig; fig2dev -L pdftex $< > $@
%.pdftex_t: %.fig; fig2dev -L pdftex_t -p $(patsubst %.fig,%.pdftex,$<) $< > $@
%.pdf: %.fig; fig2dev -L pdf $< > $@

REMOTE_USER:=dk
REMOTE_HOST:=arctic.cse.msu.edu
REMOTE_PATH:=/research/mckinley/Knoester/draft
REMOTE:=${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_PATH}

PAPER:=$(notdir ${CURDIR})
FIGURES:=$(patsubst %.eps,%.pdf,$(wildcard figures/*.eps)) \
	$(patsubst %.dot,%.pdf,$(wildcard figures/*.dot)) \
	$(patsubst %.fig,%.pdftex,$(wildcard figures/*.fig)) \
	$(patsubst %.fig,%.pdftex_t,$(wildcard figures/*.fig))
PDFFIGURES:=$(patsubst %.fig,%.pdf,$(wildcard figures/*.fig))
TARBALL:=${PAPER}.tar.gz
TAR_EXCLUDES:=--exclude *.tar --exclude *.gz \
	--exclude .git --exclude *.sparseimage --exclude src --exclude var \
	--exclude .hg --exclude data* --exclude *.mov
RSYNC_REMOTE:=rsync -avz -e ssh
RSYNC_LOCAL:=rsync -WrLptgo --exclude=*.sparseimage 

all: clean paper

crc: distclean paper

figures: ${FIGURES}

pdffigures: ${PDFFIGURES}

quick: figures
	pdflatex main.tex

paper: figures quick
	-bibtex main.aux
	pdflatex main.tex
	pdflatex main.tex

tarball: clean paper
	tar -cvhf ${PAPER}.tar ${TAR_EXCLUDES} ../${PAPER}
	gzip ${PAPER}.tar

localconfig:
	updmap --cnffile pdflatex.cfg --outputdir .

publish: tarball
	${RSYNC_REMOTE} ${TARBALL} ${REMOTE}/${TARBALL}
	ssh ${REMOTE_USER}@${REMOTE_HOST} \
		'chmod u=rw,g=r,o= ${REMOTE_PATH}/${TARBALL}'	
	@echo ${REMOTE_TARBALL}

clean:
	-rm *.aux *.log *.bbl *.blg *.tar.gz

distclean: clean
	-rm *.pdf diff-* tmp-*
	-rm figures/*.pdf figures/*.pdftex*

checkpoint:
	@echo Creating a disk image...
	hdiutil create -srcfolder ../${PAPER} -fs HFS+ -volname ${PAPER} \
		/tmp/${PAPER}
	mv /tmp/${PAPER}.dmg ./
