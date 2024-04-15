#/usr/bin/python

import sys
import os

def flatten(t):
    return [item for sublist in t for item in sublist]

###############################################################################

#all strings should start with space, for concatenation

myexec = " mixed_ll"

# Vegas
vegasopt = " -vegasNc0 8000000 -vegasNc1 16000000 -vegasIt0 10 -vegasIt1 5"

# MPI
mpicores = "60"

cmd = "mpirun -n "+mpicores+myexec
cmd = cmd+vegasopt

allsect=[]

allcorr = ['lo','nloqcd','nloewk','nnlo_rr','nnlo_rv','nnlo_vv','nnlo_s']

allch = {
    "lo":["ns","aa"],
    "nloqcd":["ns","gq","qg"],
    "nloewk":["ns","aq","qa","aa"],
    "nnlo_rr":["ns_ga","ns_qqb","ns_qq","gq","qg","aq","qa","ag","ga"],
    "nnlo_rv":["TBD"],
    "nnlo_vv":["TBD"],
    "nnlo_s" :["TBD"]
}

allsec = {
    "lo":{"ns":["na"],"aa":["na"]},
    #
    "nloqcd":{"ns":["r_is","s","v"],
              "qg":["r_is","s"],
              "gq":["r_is","s"]},
    #
    "nloewk":{"ns":["r_is","r_fs_53","r_fs_54","s","v"],
              "qa":["r_is","s"],
              "aq":["r_is","s"],
              "aa":["r_fs_53","r_fs_54","s","v"]},
    #
    "nnlo_rr":{"ns_ga":flatten([["rr_5161"+i for i in ["a","c"]], ["rr_5262"+i for i in ["a","c"]], ["rr_516"+i for i in ["2","3","4"]],
                                ["rr_526"+i for i in ["1","3","4"]]]),
               "ns_qqb":flatten([["rr_5161"+i for i in ["a","c"]],["rr_5262"+i for i in ["a","c"]]]),
               "ns_qq":flatten([["rr_5161"+i for i in ["a","c"]],["rr_5262"+i for i in ["a","c"]]]),
               #
               "gq":flatten([ ["rr_5161"+i for i in ["a","b","c","d"]], ["rr_516"+i for i in ["2","3","4"]] ]),
               "qg":flatten([ ["rr_5262"+i for i in ["a","b","c","d"]], ["rr_526"+i for i in ["1","3","4"]] ]),
               "aq":flatten([ ["rr_5161"+i for i in ["a","b","c","d"]], ["rr_5262"+i for i in ["a","b","c","d"]], ["rr_5162", "rr_5261"] ]),
               "qa":flatten([ ["rr_5161"+i for i in ["a","b","c","d"]], ["rr_5262"+i for i in ["a","b","c","d"]], ["rr_5162", "rr_5261"] ]),
               "ag":flatten([["rr_5161"+i for i in ["a","c"]], ["rr_5262"+i for i in ["a","c"]], ["rr_5162", "rr_5261"] ]),
               "ga":flatten([["rr_5161"+i for i in ["a","c"]], ["rr_5262"+i for i in ["a","c"]], ["rr_5162", "rr_5261"] ])}
    }

mycorr = {"lo":"lo",
          "nloqcd":"nloqcd",
          "nloewk":"nloewk",
          "nnlo_rr":"nnlo",
          "nnlo_rv":"nnlo",
          "nnlo_vv":"nnlo",
          "nnlo_s" :"nnlo"}

### now user-defined overwrite
allcorr = ['lo','nloqcd','nloewk','nnlo_rr']

print "echo '' > logfile.log"
for corr in allcorr:
    mychs = allch[corr]
    for ch in mychs:
        for sec in allsec[corr][ch]:
            print cmd,"-corr",mycorr[corr],"-ch",ch,"-sec",sec,">> logfile.log"
            print "echo '' >> logfile.log"

