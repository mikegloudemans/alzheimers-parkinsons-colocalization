# Munge 23andme data
bash scripts/munge_gwas/setup_23andme_data.sh

# Munge all other files
python bin/gwas-download/munge/custom_munge.py scripts/munge_gwas/munge_ad.config
