FILES = thinkstats2.py thinkstats2_test.py \
survey.py first.py risk.py  \
descriptive.py birthdays.py \
populations.py populations_cdf.py brfss.py cumulative.py hypothesis.py \
rankit.py conditional.py class_size.py \
decay.py   \
estimate.py gini.py irs.py brfss_scatter.py \
brfss_corr.py agemodel.py age_lm.py hypothesis_analytic.py \
score_example.py relay.py relay_soln.py relay_cdf.py relay_normal.py \
brfss_figs.py chi.py correlation.py social.py pareto_world.py \
khan.py blinky.py continuous.py bayes_height.py 

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
\
actual_speeds.png        nsfg_birthwgt_cdf.png \
age_scatter.png                 nsfg_birthwgt_lognormal.png \
bq_cdf1.png                     nsfg_birthwgt_model.png \
bq_cdf2.png                     nsfg_birthwgt_normal.png \
bq_cdf.png                      nsfg_birthwgt_pmf.png \
bq_gap1.png                     nsfg_diffs.png \
bq_gap2.png                     nsfg_hist.png \
brfss_change2.png               nsfg_oz_pmf.png \
brfss_change.png                nsfg_pmf.png \
brfss_scatter.png               nsfg_ratios.png \
brfss_weight_log.png            observed_speeds.png \
brfss_weight_lognormal.png      pareto_cdf.png \
brfss_weight_model.png          pareto_height.png \
brfss_weight_normal.png         populations.png \
conditional.png                 populations_loglog.png \
example_cdf.png                 populations_logx.png \
expo_cdf.png                    populations_rankit.png \
heri1.png                       posteriors.png \
heri2.png                       race_predictor1.png \
heri3.png                       race_predictor2.png \
income_loglog.png               race_predictor3.png \
income_logx.png                 race_predictor4.png \
interarrivals.png               relay_normal.png \
interarrivals_logy.png          resample_cdf.png \
interarrivals_model.png         resampled_delta_cdf.png \
length_analytic_deltas_cdf.png  scatter1.png \
length_deltas_cdf.png           scatter2.png \
locomotive.png                  scatter3.png \
marathon_cdf.png                scatter4.png \
marathon_pmf.png                sigma.png \
nlsy_scatter.png                slashdot.loglog.png \
nnsfg_birthwgt_model.png        slashdot.logx.png \
normal_cdf.png                  weight_deltas_cdf.png \
nsfg_age_cdf.png

all_figs:
	python descriptive.py
	python cumulative.py
	python babyboom.py
	python continuous.py
	python hypothesis.py
	python locomotive.py
	python brfss_figs.py
	python brfss_scatter.py

FIGDEST = ../../figs

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

red:
	rsync -avp red.htm red.nb HTML* $(DEST)

bayes_tut:
	zip -r bayes_tutorial.zip bayes_tutorial
	rsync -a bayes_tutorial.zip $(DEST)

TUT_FILES = dice.py  install_test.py   sat_ranks.csv  \
thinkbayes.py euro.py  euro2.py myplot.py  sat.py  sat_scale.csv \
train.py

bayes_tut13:
	rm -rf bayes_pycon13
	mkdir bayes_pycon13
	cd bayes_tutorial13; cp $(TUT_FILES) ../bayes_pycon13
	zip -r bayes_pycon13.zip bayes_pycon13
	rsync -a bayes_pycon13.zip $(DEST)
	cd /home/downey/public_html/greent; sh back

.PHONY: docs $(DOCPY)

docs: $(DOCPY)

$(DOCPY):
	pydoc -w ./$@

gss-table:
	cat gss.sp_male.html gss.sp_female.html gss.parent.html gss.transition.table.html gss.generation.table.html > gss.all_tables.html

clean:
	rm gss.regress.*.db

