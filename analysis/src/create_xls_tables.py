import pandas
import sys
import os

def get_ripchip_targets(get_sam_rank=False):
    ripchip_filename = 'lib/ripchip/sd01.txt'
    ripchip = pandas.read_csv(
        ripchip_filename, sep='\t')
    top_1350 = []
    by_sam_rank = {}
    for index, row in ripchip.iterrows():
        if row['Gene public name'] not in top_1350:
            top_1350.append(row['Gene public name'])
            by_sam_rank[row['Gene public name']] = row['SAM rank']
        if len(list(set(top_1350))) >= 1350:
            break
    #top_1350_rows = ripchip[0:1349]  # There is a row zero.
    #rip_targets = list(top_1350['Gene public name'])
    print "Top targets: {t}.\nLast targets: {b}.".format(
        t=str(top_1350[0:2]), b=str(top_1350[-2:]))
    if get_sam_rank:
        return set(top_1350), by_sam_rank
    return set(top_1350)


def name_cols(df):
    df['Rank'] = range(1, len(df.index)+1)
    df['Chrm'] = df['chrm']
    df['Nucleotide location left end of peak'] = df['left']
    df['Nucleotide location right end of peak'] = df['right']
    df['Gene name'] = df['gene_name']
    df['Peak height (reads/million)'] = df['height']
    df['Strand'] = df['strand']
    df['Transcript'] = df['transcript_name']
    df['Biotype'] = df['biotype']
    df['Location'] = df['location']
    df['Has UGUNNNAU (FBE)?'] = df['has_fbe']
    df['Has -1C FBE?'] = df['minus_one_c']
    df['Has -2C FBE?'] = df['minus_two_c']
    df['Has -1 or -2 C FBE?'] = df['minus_one_or_two_c']
    df['FBF-1 max coverage in peak (reads/million)'] = df['fbf1_reads']
    df['FBF-2 max coverage in peak (reads/million)'] = df['fbf2_reads']
    df['N2 control for FBF-1 max coverage in peak (reads/million)'] = df['fbf1_n2_reads']
    df['N2 control for FBF-2 max coverage in peak (reads/million)'] = df['fbf2_n2_reads']
    df['RNA-seq modencode (Acc. 4594) max coverage in peak (reads/million)'] = df['rna_seq_modencode']
    df['RNA-seq Ortiz et al., oogenic germlines max coverage in peak (reads/million)'] = df['rna_seq_oo']
    df['RNA-seq Ortiz et al., oogenic germlines (RPKM)'] = df['RNA abundance']
    df['Sequence under peak'] = df['seq']
    return df

def _mk(t):
    if os.path.isfile(t):
        t = os.path.dirname(t)
    if not os.path.exists(t):
        os.system('mkdir ' + t)


def normalize_to_rna_seq(df):
    ortiz_min = min([x for x in df.rna_seq_oo if x > 0])
    rpkm_min = min([x for x in df['RNA abundance'].tolist() if x > 0])
    modencode_min = min([x for x in df.rna_seq_modencode if x > 0])
    ortiz_norm = [float(t[0])/float(max([t[1], ortiz_min])) for t in zip(df.height, df.rna_seq_oo)]
    ortiz_rpkm_norm = [float(t[0])/float(max([t[1], rpkm_min])) for t in zip(df.height, df['RNA abundance'].tolist())]
    modencode_norm = [float(t[0])/float(max([t[1], modencode_min])) for t in zip(df.height, df.rna_seq_modencode)]
    df['Peak height/RNA abundance (RPKM for oogenic germlines, Ortiz et al.)'] = ortiz_rpkm_norm
    df['Peak height/RNA abundance (Max RNA coverage for oogenic germlines, Ortiz et al.)'] = ortiz_norm
    df['Peak height/RNA abundance (Max RNA coverage for whole worms, modencode (Acc. 4594))'] = modencode_norm

fbf1 = pandas.read_csv('filtered/combined_fbf1.txt', sep='\t')
fbf1_sep = pandas.read_csv('filtered_separate/combined_fbf1.txt', sep='\t')
fbf2 = pandas.read_csv('filtered/combined_fbf2.txt', sep='\t')
both = pandas.read_csv('filtered_five_reps/combined_fbf.txt', sep='\t')
fbf1.sort(columns='fbf1_reads', inplace=True, ascending=False)
fbf1_sep.sort(columns='fbf1_reads', inplace=True, ascending=False)
fbf2.sort(columns='fbf2_reads', inplace=True, ascending=False)
both.sort(columns='height', inplace=True, ascending=False)

targset = get_ripchip_targets()

def isin(z, aset):
    if z in aset: return 1
    else: return 0

for x in [fbf1, fbf2, fbf1_sep, both]:
    x['is_ripchip'] = [isin(y, targset) for y in x.gene_name]
    x['Is a target by RIP-chip (Kershner et al.)?'] = x['is_ripchip']
    name_cols(x)
    normalize_to_rna_seq(x)

columns = [
    'Rank', 'Chrm', 'Nucleotide location left end of peak',
    'Nucleotide location right end of peak', 'Gene name',
    'Peak height (reads/million)',
    'Strand', 'Transcript', 'Biotype', 'Location',
    'Has UGUNNNAU (FBE)?', 'Has -1C FBE?',
    'Has -1 or -2 C FBE?', 'Has -2C FBE?',
    'FBF-1 max coverage in peak (reads/million)',
    'FBF-2 max coverage in peak (reads/million)',
    'N2 control for FBF-1 max coverage in peak (reads/million)',
    'N2 control for FBF-2 max coverage in peak (reads/million)',
    'RNA-seq modencode (Acc. 4594) max coverage in peak (reads/million)',
    'RNA-seq Ortiz et al., oogenic germlines max coverage in peak (reads/million)',
    'RNA-seq Ortiz et al., oogenic germlines (RPKM)',
    'Peak height/RNA abundance (RPKM for oogenic germlines, Ortiz et al.)',
    'Peak height/RNA abundance (Max RNA coverage for oogenic germlines, Ortiz et al.)',
    'Peak height/RNA abundance (Max RNA coverage for whole worms, modencode (Acc. 4594))',
    'Is a target by RIP-chip (Kershner et al.)?',
    'Sequence under peak']
raw_columns = [
    'chrm', 'left', 'right', 'gene_name',
    'strand', 'ratio', 'biotype', 'location',
    'height', 'seq', 'control', 'fbf1_n2',
    'fbf2_reads', 'fbf1_n2_reads', 'exp_fbf1_GGTT_reads', 'iv',
    'control_reads', 'minus_one_or_two_c', 'fbf1_rep_2', 'fbf1_rep_1',
    'exp_fbf1_TGGC_reads', 'fbf1_reads', 'fbf2_rep_3', 'exp_fbf1_CGGA_reads',
    'fbf2_rep_1', 'exp', 'n2_from_fbf2_reads', 'minus_two_c',  'transcript_id',
    'fbf2_n2_reads', 'exp_fbf1_CGGA', 'number_of_fbes_fbe', 'fbf1_rep_3',
    'has_fbe', 'exp_fbf1_TGGC','minus_one_c', 'exp_fbf1_GGTT', 'Rank',
    'fbf2_rep_2',
    'n2_from_fbf1_reads', 'transcript_name', 'fbf2_n2', 'RNA abundance',
    'rna_seq_modencode', 'rna_seq_oo'
    ]

_mk('tables/')
writer = pandas.ExcelWriter('tables/Table S2.xls')
# Tab 1.
fbf1.to_excel(
    writer,
    sheet_name='FBF-1 controlled to FBF-1 N2',
    index=False,
    columns=columns,
    )
# Tab 2.
fbf1.to_excel(
    writer,
    sheet_name='FBF-1 controlled to FBF-2 N2',
    index=False,
    columns=columns,
    )
# Tab 3.
fbf2.to_excel(
    writer,
    sheet_name='FBF-2 controlled to FBF-2 N2',
    index=False,
    columns=columns,
    )
# Tab 4.
wo_ncrna = both[both['biotype']=='protein_coding']
wo_ncrna.to_excel(
    writer,
    sheet_name='FBF peaks without ncRNA',
    index=False,
    columns=columns,
    )
# Tab 5.
both.to_excel(
    writer,
    sheet_name='FBF peaks with ncRNA',
    index=False,
    columns=columns,
    )
is_rip = both[both['is_ripchip']==1]
is_not_rip = both[both['is_ripchip']==0]
# Tab 6.
is_rip.to_excel(
    writer,
    sheet_name='FBF peaks & also RIP-chip',
    index=False,
    columns=columns,
    )
# Tab 7.
is_not_rip.to_excel(
    writer,
    sheet_name='FBF peaks & not RIP-chip',
    index=False,
    columns=columns,
    )
writer.save()
