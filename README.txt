The following steps will reproduce all panels in all GALA II figures in "Subtyping asthma unravels heterogeneous genetic etiologies", Dahl et al 2024+. (With the exception of Fig 3C, which was produced by the Seibold lab.)

Step 0: set up data, global paramters, utility functions:
-global parameters:							setup.R
-convert data matrix to Rdata:	parse_gala_data/read_gala_data.R
-restructure data to Y/Yb/G/X:	parse_gala_data/finalize.R
-read in candidate SNPs:				parse_gala_data/read_asthma_gwas_snps.R 

Step 1: choose_K.R assesses cross-validated likelihood for different K (= # of subtypes) and plot with SFig1.R 

Step 2: find asthma subtypes by running RGWAS on asthma cases (run_mfmr.R) and plot subtype-specific phenotype dist'ns with Fig2.R 
-Validation: Fig3ab.R tests and plots held-out T2-relevant variables: blood eos (Fig3a) and FeNO (Fig3b) 

Step 3: test for subtype-specific effects of clinical covariates on asthma (run_droptest_mlr.R) and 2ndary traits (run_droptest.R)
-For plotting (SFig3_STab1.R), fit subtype-specific effects post-hoc (esp bc center confounds ethnicity) for asthma (run_mlrx_betas.R) and 2ndary traits (run_mfmrx_betas.R) 

Step 4: test for subtype-specific SNP effects:
-Asthma case/control GWAS (gwas.R)
-Subtype-specific GWAS on asthma with Multinomial Logistic Regression (gwas_mlr.R)
-Subtype-specific GWAS on case-only 2ndary traits with standard interactions (gwas_caseonly.R)
-qqplots.R plots the GWAS ouputs (uses functions.R)
-Fig4AB.R draws forest plots for significant subtype-specific SNP-trait pairs

Step 5: subtype-specific PRS effects: first build PRS (build_prs.sh); then test subtype-specificity with MLR (prs_mlr.R) and plot (Fig4C.R) 

Step 6: Estimate subtype-specific heritability for 2ndary features (gxemm.R) and plot with Fig4D.R

Alternative Step 2: use k-means instead of RGWAS (run_kmeans.R), then repeat GWAS (Step 4), then compare QQ plots across k=2,3,4 (kmeans_plots.R)


Additionally, code to reproduce the simulations in the Supplement lives in the folder rgwas_mlr_sims/
