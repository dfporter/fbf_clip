w_fbe = ''
wo_fbe = ''
import re
pat = re.compile('TGT\w\w\wAT', re.IGNORECASE)
with open(infile, 'r') as f:
    while True:
        name = f.readline()
        seq = f.readline()
        if not name: break
        if pat.search(seq):
            w_fbe += name + seq
        else:
            wo_fbe += name + seq
import os
open(os.path.basename(infile).partition('fa')[0] + '_fbes_removed_by_hand.fa,
     'w').write(wo_fbe)
open(os.path.basename(infile).partition('fa')[0] + '_fbes_included_by_hand.fa,
     'w').write(w_fbe)
