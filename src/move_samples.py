
out_dir = 'pypeaks_fdr5_negip_local/'
if not os.path.exists(out_dir):
    os.system('mkdir ' + out_dir)
os.system('cp fbf1/peaks/combined_fbf1/4 {od}/combined_fbf1.txt\n'.format(od=out_dir))
os.system('cp fbf2/peaks/combined_fbf2/4 {od}/combined_fbf2.txt\n'.format(od=out_dir))
os.system('cp fbf1/peaks/runx_fbf_aacc_20mapq/null_hyp_4.txt {od}/fbf1_aacc.txt\n'.format(od=out_dir))
os.system('cp fbf1/peaks/runx_fbf_gcca_20mapq/null_hyp_4.txt {od}/fbf1_gcca.txt\n'.format(od=out_dir))
os.system('cp fbf1/peaks/runx_fbf_tccg_20mapq/null_hyp_4.txt {od}/fbf1_tccg.txt\n'.format(od=out_dir))

os.system('cp fbf2/peaks/run813_fbf_aacc_20mapq/null_hyp_4.txt {od}/fbf2_aacc.txt\n'.format(od=out_dir))
os.system('cp fbf2/peaks/run813_fbf_gcca_20mapq/null_hyp_4.txt {od}/fbf2_gcca.txt\n'.format(od=out_dir))
os.system('cp fbf2/peaks/run813_fbf_tccg_20mapq/null_hyp_4.txt {od}/fbf2_tccg.txt\n'.format(od=out_dir))


fbf1_dir = 'fbf1/peaks/'
fbf2_dir = 'fbf1/peaks/'

