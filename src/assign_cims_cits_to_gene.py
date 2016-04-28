from get_fasta_of_dCLIP_run_DREME import *

def assign_cims_to_gene(cims_file, gtf_df, label='',
                        out_dir='cims_out_with_genes/',
                        cits=False):
    peaks = pandas.read_csv(cims_file, sep='\t')
    #peaks = peaks[peaks['p_value_local']<1e-5]
    #peaks = peaks[peaks['p_value_neg']<1e-5]
    split_and_assign_to_gene(peaks, gtf_df)
    if not cits:
        peaks.sort(['FDR', 'tagNumber_k'], ascending=[0,0], inplace=True)
    else:
        peaks.sort('score', ascending=0, inplace=True)
    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)
    peaks.to_csv('%s/%s' % (out_dir, label), sep='\t')
    return peaks

def assign_all_in_dir(top_dir, cits=False, out_dir='cims_out_with_genes'):
    gtf_df = pandas.read_csv('../clip/lib/gtf_with_names_column.txt', sep='\t')
    for cims_file in glob.glob(top_dir + '/*'):
        print "Assigning %s to genes..." % cims_file
        assign_cims_to_gene(
            cims_file, gtf_df, fdr=0.001,
            out_dir=out_dir,
            label=os.path.basename(cims_file), cits=cits)

def assign_table(peaks, gtf_df=False, given_gtf=False):
    if not given_gtf:
        gtf_df = pandas.read_csv('../clip/lib/gtf_with_names_column.txt', sep='\t')
    split_and_assign_to_gene(peaks, gtf_df)
    
if __name__ == '__main__':
    gtf_df = pandas.read_csv('../clip/lib/gtf_with_names_column.txt', sep='\t')
    for cims_file in glob.glob('cims_out/*'):
        print "Assigning %s to genes..." % cims_file
        assign_cims_to_gene(
            cims_file, gtf_df, fdr=0.001,
            label=os.path.basename(cims_file))
