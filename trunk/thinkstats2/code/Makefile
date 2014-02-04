FILES = thinkstats2.py thinkstats2_test.py thinkplot.py \
survey.py first.py risk.py  \
descriptive.py birthdays.py \
populations.py populations_cdf.py brfss.py cumulative.py hypothesis.py \
rankit.py conditional.py class_size.py \
decay.py   \
estimate.py gini.py irs.py brfss_scatter.py \
brfss_corr.py agemodel.py age_lm.py hypothesis_analytic.py \
score_example.py relay.py relay_soln.py relay_cdf.py relay_normal.py \
brfss_figs.py chi.py correlation.py social.py pareto_world.py \
khan.py blinky.py continuous.py bayes_height.py pdf_example.py

DATA = 08in11si.csv 

thinkstats2:
	svn cp $(FILES) ../../thinkstats2/code 
	svn cp $(DATA) ../../thinkstats2/code 

DOCS = thinkstats2.html thinkplot.html

DOCPY = thinkstats2.py thinkplot.py

DEST = /home/downey/public_html/greent/thinkstats2

FIGS = actual_speeds.pdf    length_deltas_cdf.pdf     observed_speeds.pdf \
age_scatter.pdf             locomotive.pdf            pareto_cdf.pdf \
bq_gap1.pdf                 marathon_cdf.pdf          pareto_height.pdf \
bq_gap2.pdf                 marathon_pmf.pdf          populations_loglog.pdf \
brfss_weight_lognormal.pdf  normal_cdf.pdf            populations_logx.pdf \
brfss_weight_log.pdf        nsfg_age_cdf.pdf          populations.pdf \
brfss_weight_model.pdf      nsfg_birthwgt_cdf.pdf     populations_rankit.pdf \
brfss_weight_normal.pdf     nsfg_birthwgt_model.pdf   race_predictor4.pdf \
example_cdf.pdf             nsfg_birthwgt_normal.pdf  scatter1.pdf \
expo_cdf.pdf                nsfg_birthwgt_pmf.pdf     scatter2.pdf \
interarrivals_logy.pdf      nsfg_diffs.pdf            scatter3.pdf \
interarrivals_model.pdf     nsfg_hist.pdf             scatter4.pdf \
interarrivals.pdf           nsfg_pmf.pdf              slashdot.logx.pdf \
pdf_example.pdf \

all_figs:
	python descriptive.py
	python cumulative.py
	python babyboom.py
	python continuous.py
	python hypothesis.py
	python locomotive.py
	python brfss_figs.py
	python brfss_scatter.py

FIGDEST = ../figs

%.html: %.py
	pydoc -w $<

code:
	zip -r thinkstats2.code.zip $(FILES)
	rsync -a thinkstats2.code.zip $(FILES) $(DATA) $(DEST)
	rsync -a $(DOCS) $(DEST)
	chmod -R o+r $(DEST)/*
	cd /home/downey/public_html/greent; sh back

figs:
	rsync -a $(FIGS) $(FIGDEST)

.PHONY: docs $(DOCPY)

docs: $(DOCPY)

$(DOCPY):
	pydoc -w ./$@

