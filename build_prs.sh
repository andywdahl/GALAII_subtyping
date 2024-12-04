echo "0.000000000001 0 0.000000000001" > range_list
echo "0.0000000001 0 0.0000000001" >> range_list
echo "0.00000001 0 0.00000001" >> range_list
echo "0.000001 0 0.000001" >> range_list
echo "0.00001 0 0.00001" >> range_list
echo "0.0001 0 0.0001" >> range_list
echo "0.001 0 0.001" >> range_list
echo "0.01 0 0.01" >> range_list
echo "0.1 0 0.1" >> range_list
echo "1.0 0 1.0" >> range_list

for phen in disease_ALLERGY_ECZEMA_DIAGNOSED blood_EOSINOPHIL_COUNT disease_ASTHMA_DIAGNOSED mental_NEUROTICISM
do

	cut -f1 parse_gala_data/sumstats/${phen}.sumstats | uniq -d > rs_sumstats		 ## filter multi-allelic
	awk '{print $1,$10}' parse_gala_data/sumstats/${phen}.sumstats  | uniq -u -f 1 > SNP.pvalue

	~/misc/plink \
		--bfile parse_gala_data/data/GALAII_all_freeze_041112 \
		--score parse_gala_data/sumstats/${phen}.sumstats 1 4 8 header \
		--q-score-range range_list SNP.pvalue \
		--exclude rs_sumstats \
		--out parse_gala_data/parsed_data/${phen}

done
