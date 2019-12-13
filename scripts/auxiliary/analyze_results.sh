### NOTE: This isn't a very clean script, and you'll have to change the 
### output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/ 
### directory to wherever your colocalization results are stored.


# How many suggestive loci per study?
cut -f1,2,4 output/test-snps/brain-gwas_all-gwas_gwas-pval1e-05_gwas-window500000_snps-considered.txt | sort -k3,3 -k1,1 -k2,2n | cut -f3 | sort | uniq -c

# How many total tests run across all GTEx brain tissues and all GWAS?
cat output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*_clpp_status.txt | wc -l

# How many after filtering?
cat output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*_clpp_status.txt | awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5)) print $0}' | wc -l

###########################
# Individual studies
###########################

for s in `ls output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*clpp*`; do

	study=`basename $s`
	echo $study

	# How many total tests were run?
	cat output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | wc -l

	# How many after filtering?
	awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | wc -l

	# How many are actually colocalized?
	awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | wc -l

	# How many genes colocalized?
	awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f4 | sort | uniq | wc -l

	# How many loci with colocalization?
	awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f1 | sort | uniq | wc -l

	# Venn diagram stuff
	# pairs
	for s2 in `ls output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*clpp*`; do
		study2=`basename $s2`

		echo -e "$study\n$study2"
		# How many genes found in multiple studies?
		join <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f4 | sort | uniq) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study2 | cut -f4 | sort | uniq) | wc -l

		# How many loci found in multiple studies?
		join <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f1 | sort | uniq) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study2 | cut -f1 | sort | uniq) | wc -l

		# triples
		for s3 in `ls output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*clpp*`; do
	                study3=`basename $s3`
			echo -e "$study\n$study2\n$study3"

			join <(join <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f4 | sort | uniq) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study2 | cut -f4 | sort | uniq)) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study3 | cut -f4 | sort | uniq) | wc -l

			join <(join <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f1 | sort | uniq) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study2 | cut -f1 | sort | uniq)) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study3 | cut -f1 | sort | uniq) | wc -l
		done
	done

 	# Tissue specificity - which tissues is each colocalized gene found in?
	awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | sort -k4,4 -k6,6r | cat

done


###########################
# Individual studies
# Exact same thing, except with CLPP-mod
###########################

echo -e "Beginning of CLPP mod results"

for s in `ls output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*clpp*`; do

	study=`basename $s`
	echo $study

	# How many total tests were run?
	cat output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | wc -l

	# How many after filtering?
	awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | wc -l

	# How many are actually colocalized?
	awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | wc -l

	# How many genes colocalized?
	awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f4 | sort | uniq | wc -l

	# How many loci with colocalization?
	awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f1 | sort | uniq | wc -l

	# Venn diagram stuff
	# pairs
	for s2 in `ls output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*clpp*`; do
		study2=`basename $s2`

		echo -e "$study\n$study2"
		# How many genes found in multiple studies?
		join <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f4 | sort | uniq) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study2 | cut -f4 | sort | uniq) | wc -l

		# How many loci found in multiple studies?
		join <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f1 | sort | uniq) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study2 | cut -f1 | sort | uniq) | wc -l

		# triples
		for s3 in `ls output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*clpp*`; do
	                study3=`basename $s3`
			echo -e "$study\n$study2\n$study3"

			join <(join <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f4 | sort | uniq) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study2 | cut -f4 | sort | uniq)) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study3 | cut -f4 | sort | uniq) | wc -l

			join <(join <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | cut -f1 | sort | uniq) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study2 | cut -f1 | sort | uniq)) <(awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study3 | cut -f1 | sort | uniq) | wc -l
		done
	done

 	# Tissue specificity - which tissues is each colocalized gene found in?
	awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/$study | sort -k4,4 -k6,6r | cat

done



for f in `ls output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*clpp*`; do echo $f; awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' $f | cut -f4 | sort | uniq; done

for f in `ls output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*clpp*`; do echo $f; awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($6 >= 0.01)) print $0}' $f | cut -f1 | sort | uniq; done

for f in `ls output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*clpp*`; do echo $f; awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' $f | cut -f4 | sort | uniq; done

for f in `ls output/colocalization/2019-04-25_13-55-25.207812_toms_brain_study_new_ad/*clpp*`; do echo $f; awk '{if (($5 >= 50) && ($7 >= 5) && ($8 >= 5) && ($10 >= 0.3)) print $0}' $f | cut -f1 | sort | uniq; done
