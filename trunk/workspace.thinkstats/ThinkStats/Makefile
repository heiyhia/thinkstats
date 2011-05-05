FILES = survey.py first.py risk.py Pmf.py Pmf_test.py \
thinkstats.py thinkstats_test.py descriptive.py birthdays.py \
populations.py populations_cdf.py brfss.py cumulative.py hypothesis.py \
coin.py coin2.py rankit.py conditional.py class_size.py \
decay.py Cdf.py Cdf_test.py locomotive.py myplot.py CDBRFS08.ASC.gz \
erf.py estimate.py gini.py 08in11si.csv irs.py brfss_scatter.py \
brfss_corr.py agemodel.py age_lm.py hypothesis_analytic.py \
score_example.py relay.py relay_soln.py relay_cdf.py relay_normal.py \
brfss_figs.py chi.py correlation.py babyboom.dat

DOCS = survey.html first.html risk.html Pmf.html thinkstats.html \
descriptive.html myplot.html brfss.html hypothesis.html cumulative.html \
Cdf.html

DOCPY = survey.py first.py risk.py Pmf.py thinkstats.py \
descriptive.py myplot.py brfss.py hypothesis.py cumulative.py \
Cdf.py

DEST = /home/downey/public_html/greent/thinkstats

FIGS = actual_speeds.eps        nsfg_birthwgt_cdf.eps \
age_scatter.eps                 nsfg_birthwgt_lognormal.eps \
bq_cdf1.eps                     nsfg_birthwgt_model.eps \
bq_cdf2.eps                     nsfg_birthwgt_normal.eps \
bq_cdf.eps                      nsfg_birthwgt_pmf.eps \
bq_gap1.eps                     nsfg_diffs.eps \
bq_gap2.eps                     nsfg_hist.eps \
brfss_change2.eps               nsfg_oz_pmf.eps \
brfss_change.eps                nsfg_pmf.eps \
brfss_scatter.eps               nsfg_ratios.eps \
brfss_weight_log.eps            observed_speeds.eps \
brfss_weight_lognormal.eps      pareto_cdf.eps \
brfss_weight_model.eps          pareto_height.eps \
brfss_weight_normal.eps         populations.eps \
conditional.eps                 populations_loglog.eps \
example_cdf.eps                 populations_logx.eps \
expo_cdf.eps                    populations_rankit.eps \
heri1.eps                       posteriors.eps \
heri2.eps                       race_predictor1.eps \
heri3.eps                       race_predictor2.eps \
income_loglog.eps               race_predictor3.eps \
income_logx.eps                 race_predictor4.eps \
interarrivals.eps               relay_normal.eps \
interarrivals_logy.eps          resample_cdf.eps \
interarrivals_model.eps         resampled_delta_cdf.eps \
length_analytic_deltas_cdf.eps  scatter1.eps \
length_deltas_cdf.eps           scatter2.eps \
locomotive.eps                  scatter3.eps \
marathon_cdf.eps                scatter4.eps \
marathon_pmf.eps                sigma.eps \
nlsy_scatter.eps                slashdot.loglog.eps \
nnsfg_birthwgt_model.eps        slashdot.logx.eps \
normal_cdf.eps                  weight_deltas_cdf.eps \
nsfg_age_cdf.eps \
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

FIGDEST = ../../figs

%.html: %.py
	pydoc -w $<

code:
	rsync -a $(FILES) $(DEST)
	rsync -a $(DOCS) $(DEST)
	chmod -R o+r $(DEST)/*

figs:
	rsync -a $(FIGS) $(FIGDEST)

.PHONY: docs $(DOCPY)

docs: $(DOCPY)

$(DOCPY):
	pydoc -w ./$@