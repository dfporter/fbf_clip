import pandas
import sys
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

#df = pandas.read_excel(sys.argv[1])
df = pandas.read_csv('filtered_five_reps/combined_fbf.txt',
                     sep='\t')
name_cols(df)

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
def get_max_heights(df):
    r = df.to_dict('records')
    import collections
    x = collections.defaultdict(float)
    for row in r:
        x[row['Gene name']] = max([
            row['Peak height (reads/million)'], x[row['Gene name']]])
    return x

def get_max_cds_heights(df):
    x = df[df['Location']=='CDS']
    return get_max_heights(x)

def get_max_three_heights(df):
    x = df[df['Location']=="3'UTR"]
    return get_max_heights(x)

def get_max_five_heights(df):
    x = df[df['Location']=="5'UTR"]
    return get_max_heights(x)

print df
max_heights = get_max_heights(df)
max_three_heights = get_max_three_heights(df)
max_five_heights = get_max_five_heights(df)
max_cds_heights = get_max_cds_heights(df)

q ={}
q['five'] = df[df['Location']=="5'UTR"]
q['cds'] = df[df['Location']=="CDS"]
q['three'] = df[df['Location']=="3'UTR"]

has = {}
has['five'] = set(q['five']['Gene name'].tolist())
has['cds'] = set(q['cds']['Gene name'].tolist())
has['three'] = set(q['three']['Gene name'].tolist())


writer = pandas.ExcelWriter('tables/Table S4 revised.xls')
# Tab: peaks in 5'UTRs
tab1_columns = columns
q['five'].sort(
    ['Peak height (reads/million)'], ascending=False, inplace=True)
q['five'].to_excel(
    writer,
    sheet_name="tab1 Targets with 5'UTR peaks",
    index=False,
    columns=columns)
q['cds'].sort(
    ['Peak height (reads/million)'], ascending=False, inplace=True)
q['cds'].to_excel(
    writer,
    sheet_name="tab2 Targets with CDS peaks",
    index=False,
    columns=columns)
q['3and5'] = pandas.DataFrame(
    [
        dict([('Gene name', x),
              ('Maximum peak height (reads/million)', max_heights[x]),
              ("Maximum 5'UTR peak height (reads/million)", max_five_heights[x]),
              ("Maximum 3'UTR peak height (reads/million)", max_three_heights[x])
              ]) \
        for x in set(has['five'] & has['three'])
    ])
q['3and5'].sort(
    ['Maximum peak height (reads/million)'], ascending=False, inplace=True)
q['3and5'].to_excel(
    writer,
    sheet_name="tab3 Trgts with 3'UTR and 5'UTR",
    index=False,
    columns=['Gene name', 'Maximum peak height (reads/million)',
             "Maximum 5'UTR peak height (reads/million)",
             "Maximum 3'UTR peak height (reads/million)"])
#
q['3andCDS'] = pandas.DataFrame(
    [
        dict([('Gene name', x),
              ('Maximum peak height (reads/million)', max_heights[x]),
              ("Maximum CDS peak height (reads/million)", max_cds_heights[x]),
              ("Maximum 3'UTR peak height (reads/million)", max_three_heights[x])
              ]) \
        for x in set(has['three'] & has['cds'])
    ])
q['3andCDS'].sort(
    ['Maximum peak height (reads/million)'], ascending=False, inplace=True)
q['3andCDS'].to_excel(
    writer,
    sheet_name="tab4 Trgts with 3'UTR and CDS",
    index=False,
    columns=['Gene name','Maximum peak height (reads/million)',
             "Maximum CDS peak height (reads/million)",
             "Maximum 3'UTR peak height (reads/million)",])
#
q['5andCDSand3'] = pandas.DataFrame(
    [
        dict([('Gene name', x),
              ('Maximum peak height (reads/million)', max_heights[x]),
              ("Maximum 5'UTR peak height (reads/million)", max_five_heights[x]),
              ("Maximum CDS peak height (reads/million)", max_cds_heights[x]),
              ("Maximum 3'UTR peak height (reads/million)", max_three_heights[x])
              ]) \
        for x in set(has['three'] & has['cds'] & has['five'])
    ])
q['5andCDSand3'].sort(
    ['Maximum peak height (reads/million)'], ascending=False, inplace=True)
q['5andCDSand3'].to_excel(
    writer,
    sheet_name="tab5 Trgts with 5', 3' and CDS",
    index=False,
    columns=['Gene name','Maximum peak height (reads/million)',
             "Maximum 5'UTR peak height (reads/million)",
             "Maximum CDS peak height (reads/million)",
             "Maximum 3'UTR peak height (reads/million)",])
#
writer.save()
