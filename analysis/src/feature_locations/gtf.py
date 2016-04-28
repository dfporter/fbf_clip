
class gtf(object):

    def __init__(self, filename):
        self.file = filename
        self.data = pandas.read_csv(filename, sep='\t')
        self.order = [
            '0', '1', '2', '3', '4', '5', '6', '7', '8',
            'gene_name', 'transcript_id', 'transcript_name',
            'exon_number', 'gene_id', 'biotype']
        self.data_list = convert_to_list(self.data)
        self.utrs = [x for x in self.data_list if x[2]=='UTR']
        self.exon = [x for x in self.data_list if x[2]=='exon']
        self.cds = [x for x in self.data_list if x[2]=='cds']

    def get_exons(self):
        data = self.data  # Concision.
        txpt_to_gene = dict(zip(data.transcript_id, data.gene_name))
        extents = zip(data.transcript_id, data['2'].tolist(),
                    data['4'].tolist(), data['5'].tolist())

    def convert_to_list(self, df):
        _list = df.to_dict('records')
        _list = [[x[t] for t in self.order] for x in _list]
        return _list

    def get_lists(self):
        return (self.utr, self.exon, self.cds, [])

    def flip_minus_strand_features(self, chr_len):
        
#        utrs = [tup for tup in extents if tup[1]=='UTR']
#        exon = [tup for tup in extents if tup[1]=='exon']
#        cds = [tup for tup in extents if tup[1]=='CDS']
#        txpt_to_exon = zip(data.transcript_id, d
