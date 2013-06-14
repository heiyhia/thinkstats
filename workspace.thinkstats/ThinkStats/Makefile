FILES = survey.py first.py risk.py Pmf.py Pmf_test.py \
thinkstats.py thinkstats_test.py descriptive.py birthdays.py \
populations.py populations_cdf.py brfss.py cumulative.py hypothesis.py \
coin.py coin2.py rankit.py conditional.py class_size.py \
decay.py Cdf.py Cdf_test.py locomotive.py myplot.py  \
erf.py estimate.py gini.py 08in11si.csv irs.py brfss_scatter.py \
brfss_corr.py agemodel.py age_lm.py hypothesis_analytic.py \
score_example.py relay.py relay_soln.py relay_cdf.py relay_normal.py \
brfss_figs.py chi.py correlation.py social.py pareto_world.py \
khan.py blinky.py continuous.py bayes_height.py

DATA = CDBRFS08.ASC.gz soc-Slashdot0902.txt.gz babyboom.dat \
sat_ranks.csv sat_scale.csv populations.csv birthdays.csv \
showcases.2011.csv showcases.2012.csv

DOCS = survey.html first.html risk.html Pmf.html thinkstats.html \
descriptive.html myplot.html brfss.html hypothesis.html cumulative.html \
Cdf.html

DOCPY = survey.py first.py risk.py Pmf.py thinkstats.py \
descriptive.py myplot.py brfss.py hypothesis.py cumulative.py \
Cdf.py

DEST = /home/downey/public_html/greent/thinkstats

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
	zip -r thinkstats.code.zip $(FILES)
	rsync -a thinkstats.code.zip $(FILES) $(DATA) $(DEST)
	rsync -a $(DOCS) $(DEST)
	chmod -R o+r $(DEST)/*

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

BAYES_CODE = train2.py train.py dice.py m_and_m.py monty2.py thinkbayes.py \
monty.py cookie2.py cookie.py train3.py dungeons.py hockey.py hockey_data.csv \
sat.py euro*.py kidney.py species.py paintball.py myplot.py columns.py \
thinkstats.py sat_ranks.csv sat_scale.csv correlation.py \
journal.pone.0047712.s001.csv variability.py price.py redline.py \
redline_data.py thinkplot.py

BAYES_DEST = /home/downey/public_html/greent/thinkbayes

BAYES_FIGS = train*.pdf train*.eps \
euro*.pdf euro*.eps \
dungeons*.pdf dungeons*.eps \
hockey*.pdf hockey*.eps \
sat*.eps sat*.pdf \
kidney*.eps kidney*.pdf \
species*.pdf species*.eps \
paintball*.pdf paintball*.eps \
variability*.pdf variability*.eps \
price*.pdf price*.eps \
redline*.pdf redline*.eps \


BAYES_TABLES = kidney_table*.tex

bayes_run:
	python cookie.py
	python cookie2.py
	python monty.py
	python monty2.py
	python m_and_m.py
	python dice.py
	python train.py
	python train3.py
	python dungeons.py
	python price.py
	python hockey.py
	python redline.py
	python variability.py
	python paintball.py
	python euro.py
	python euro3.py
	python sat.py
	python kidney.py
	python species.py

bayes_code:
	rm -rf thinkbayes_code
	mkdir thinkbayes_code
	cp $(BAYES_CODE) thinkbayes_code
	zip -r thinkbayes_code.zip thinkbayes_code
	rsync -a thinkbayes_code.zip $(BAYES_DEST)
	pydoc -w thinkbayes
	rsync -a thinkbayes.html $(BAYES_DEST)
	rsync -a $(BAYES_CODE) $(BAYES_DEST)
	chmod -R o+r $(BAYES_DEST)/*
	cd /home/downey/public_html/greent; sh back

bayes_figs:
	rsync -a $(BAYES_FIGS) ../../thinkbayes/figs

bayes_tables:
	rsync -a $(BAYES_TABLES) ../../thinkbayes
