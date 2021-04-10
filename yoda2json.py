
"""
python yoda2json.py \
  fused.yoda \
  30000 \
  '/MCTTBBDecayed/onelep_eq5j_eq2jb/njl,/MCTTBBDecayed/onelep_eq5j_eq3jb/njl,/MCTTBBDecayed/onelep_eq5j_ge4jb/njl,/MCTTBBDecayed/onelep_ge6j_eq2jb/njl,/MCTTBBDecayed/onelep_ge6j_eq3jb/njl,/MCTTBBDecayed/onelep_ge6j_ge4jb/njl' \
  'MUR2_MUF1_PDF261000,btagrate,ctagrate,jes'
"""

# TODO
# up/down variation pairs

import yoda
import json
from sys import argv, stdout

myargs = argv[1:]

infile = myargs.pop(0)
lumi = float(myargs.pop(0))
histnames = myargs.pop(0).split(",")
varnames = myargs.pop(0).split(",")


def lmap(f, l): return list(map(f, l))

def dmap(f, d): return { k : f(k, v) for k, v in d.items() }

def vmap(f, d): return { k : f(v) for k, v in d.items() }

def fmap(f, x): return (f(x[0]), vmap(f, x[1]))

def fnom(f, x): return (f(x[0]), x[1])

def fvars(f, x): return (x[0], f(x[1]))

def id(x): return x

def nom(x): return x[0]

def vars(x): return x[1]

def scale(x):
  def f(h):
    return lmap(lambda y : x*y, h)

  return f

def fold(f, z, h): return reduce(f, h, z)

def integral(h): return fold(float.__add__, 0.0, h)

def normalize(n, h):
  nstart = integral(h)
  return scale(n / nstart)(h)

def merge(c, d):
  d1 = d.copy()
  d1.update(c)
  return d1

def joinWith(f, c, d):
  # make sure we get all the keys from either c or d.
  # will throw an error if the keys aren't the same.
  n = f(nom(c), nom(d))
  vc = vars(c)
  vd = vars(d)
  z = vc.copy()
  z.update(vd)
  return (n , { k : f(vc[k], vd[k]) for k in z.keys() })



def concat(cs): return reduce(list.__add__, cs)

def toList(h):
  return map(yoda.HistoBin1D.area, h.bins())


def readHists(d, names):
  return concat([toList(d[name]) for name in names])

def readHistsWithVars(d, names, varnames):
  return \
    ( readHists(d, names)
    , { varname : readHists(d, [name + "[%s]" % varname for name in names]) for varname in varnames }
    )



aos = yoda.read(infile)

histvariations = readHistsWithVars(aos, histnames, varnames)


def setMin(c):
  def f(h):
    return lmap(lambda x: max(c, x), h)

  return f

def setKey(k, v, d):
  d[k] = v
  return d


hists = fmap(setMin(1e-5), histvariations)
norm = integral(nom(hists))
hists = fmap(lambda h: normalize(norm, h), hists)
normsys = scale(1.1)(nom(hists))

hists = fvars(lambda vs: setKey("ttbarnorm", normsys, vs), hists)

output = {}
output["Data"] = lmap(round, scale(lumi)(nom(hists)))
output["Nominal"] = \
  { "Bkgs" : {"ttjets" : nom(hists) }
  , "Sig" : [0]
  , "Mig" : [scale(0)(nom(hists))]
  , "Lumi" : lumi
  }

def tovar(h):
  d = \
    { "InitialValue" : 0
    , "Prior" : { "Normal" : { "Mu" : 0, "Sigma" : 1 } } 
    }

  d["Variation"] = { "Bkgs" : {"ttjets" : h } }

  return d


output["ModelVars"] = vmap(tovar, vars(hists))

json.dump(output, stdout, indent=2, sort_keys=True)
