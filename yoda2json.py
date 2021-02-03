
"""
python yoda2json.py \
  sherpa_ttbb_finalscales_decOff_hps8.MUR1_MUF1_MUQ1_PDF260400.yoda \
  30000 \
  '/hxswg_ttbjets_stable_v3_25/0_logPT_J1_(25)' \
  '/hxswg_ttbjets_stable_v3_25/1_logPT_J1_(25)' \
  '/hxswg_ttbjets_stable_v3_25/2_logPT_J1_(25)' \
  'MUR2_MUF1_PDF260400'
"""

# TODO
# up/down variation pairs

import yoda
import json
from sys import argv, stdout

myargs = argv[1:]

infile = myargs.pop(0)
lumi = float(myargs.pop(0))
ttjetsname = myargs.pop(0)
ttbname = myargs.pop(0)
ttbbname = myargs.pop(0)
varnames = myargs


aos = yoda.read(infile)

def variations(basename, varnames):
  return \
    ( aos[basename]
    , { varname : aos[basename + "[%s]" % varname] for varname in varnames }
    )

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


def concat(c, d): return joinWith(list.__add__, c, d)

def toList(h):
  return map(yoda.HistoBin1D.area, h.bins())

ttjets = fmap(toList, variations(ttjetsname, varnames))
ttb = fmap(toList, variations(ttbname, varnames))
ttbb = fmap(toList, variations(ttbbname, varnames))


def setMin(c):
  def f(h):
    return lmap(lambda x: max(c, x), h)

  return f

def setKey(k, v, d):
  d[k] = v
  return d


hists = fmap(setMin(1e-5), concat(concat(ttjets, ttb), ttbb))
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
