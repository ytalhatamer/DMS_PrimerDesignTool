import pandas as pd
from globalvars import nucsequence, protsequence, \
    tm_target, genename, nucseqcount
from functions import namer, tmcalculator, checktm, gccontent, reversecomplement,\
                    trimleft, trimright, trimboth, addleft, addboth, addright,\
                    logic, justtrim, score


df = pd.DataFrame([], columns=['Residue #', 'Primer Name', 'Primer Seq(Fwd)',
                               'Primer Sequence (Rev-Comp)', 'Tm', 'GC Content(%)',
                               'Primer Length', 'GC Score', 'Tm Score',
                               'Length Score', 'Wing Score', 'Total Score'])

primers = dict(right=dict(), left=dict(), both=dict())

aprimers = dict (right=dict(), left=dict(), both=dict())

