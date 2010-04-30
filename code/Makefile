FILES = survey.py first_baby.py Cdf.py \
2002FemPreg.dat.gz  2002FemResp.dat.gz  2002Male.dat.gz

DEST = /home/downey/public_html/greent/thinkstats

all:
	cp $(FILES) $(DEST)

WORKSPACE = /home/downey/workspace.thinkstats/ThinkStats

grab:
	rsync -av --progress $(WORKSPACE)/*.py .