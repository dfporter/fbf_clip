import numpy as np

def info(arr):
    if type(arr) == type([]):
        arr = np.array(arr)
    if type(arr) != type(np.array([])):
        print "arrinfo.info: type was not array or list."
        return ''
    print "Len {l} Min {a} Max {b} Median {med} Mean {mean} Std {s}".format(
        l=len(arr), a=np.amin(arr), b=np.amax(arr), med=np.median(arr),
        mean=np.mean(arr), s=np.std(arr))
