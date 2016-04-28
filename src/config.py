import ConfigParser
#import sys
import os
import re

def config(filepath=None):
    Config = ConfigParser.ConfigParser()
    src_dir = os.path.dirname(os.path.realpath(__file__))
    if filepath is None:
        print "src_dir for config.py is {a}".format(a=src_dir)
        Config.read(src_dir + "/config.ini")
    else:
        Config.read(filepath)
    lib = ConfigSectionMap('library', Config)
    lib['clip_replicate'] = [lib[x] for x in lib.keys() if\
                                 re.match('clip_replicate.*', x)]
    lib['clip_replicate_bed'] = [lib[x] for x in lib.keys() if\
                                 re.match('exp_bed.*', x)]
    if 'positive_control_genes' in lib:
        lib['positive_control_genes'] = lib['positive_control_genes'].split(
            ',')
    if 'experiment_name' not in lib:
        lib['experiment_name'] = os.path.basename(__file__)
    return lib

 
def ConfigSectionMap(section, Config):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1
