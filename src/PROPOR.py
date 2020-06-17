"""Confidence intervals for proportions"""
#/***********************************************************************
# * Licensed Materials - Property of IBM
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 1989, 2020
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp.
# ************************************************************************/"""Partial Least Squares Regression Module"""

from extension import Template, Syntax, setnegativedefaults
import spss, spssaux
from math import sqrt
import random, sys

__version__ = "0.9.1"
__author__ = "JKP"

# History
# 19-Feb-2008 Initial version
# 27-Feb-2008 Handle p=0 case

# Help text
helptext = """Confidence Intervals for Proportions and Differences in Proportions.

PROPOR /HELP displays this help and does nothing else.

Syntax:

PROPOR NUM=list DENOM=list [ID=varname] 
[/DATASET NAME=dsname]
[/LEVEL ALPHA=value]
[/HELP]

Example:
PROPOR NUM= 55 DENOM=100.
displays Binomial and Poisson confidence intervals for p=55/100.

The list value can be one or more numbers or a variable name.  If a name,
the variable values are used.  

There must be the same number of values for NUM and DENOM
except that a single number will be considered to be repeated as many times as needed.

ID can name a variable to identifiy the proportions.  It can only be used when list is a variable.

The following two forms are equivalent.

propor num=52 55 denom=100 100.

DATA LIST FREE /n d.
begin data.
52 100
55 100
end data.
propor num=n denom=d.

CI's will be produced for p=52/100, p=55/100 and for the difference in the two proportions.

Proportion difference CI's are for the first proportion(p0) compared individually to each other proportion.

alpha defaults to .05.

By default, the active dataset is used for variable values.  /DATASET NAME can specify a different dataset to be used.

This procedure does not aggregate the data.  If the data are not in appropriate form with numerator and denominator
portions, AGGREGATE can be used to construct the appropriate input.

Example,
DATASET DECLARE edbyregion.
AGGREGATE  /OUTFILE='edbyregion'  /BREAK=region
  /posths_sum=SUM(posths)  /N=N.
PROPOR num=posths_sum denom=n id=region /dataset name=edbyregion.
"""

def Run(args):
    """Execute the PROPOR command"""
    
    debug = False
    if debug:
        print(args)   #debug
    args = args[list(args.keys())[0]]
    # Note that the keys of args are the names of the subcommands that were given.
    if debug:
        print(args)
    
    # define the syntax
    oobj = Syntax([
        Template("NUM", subc="",  ktype="str", islist=True),
        Template("DENOM", subc="",  ktype="str", islist=True),
        Template("ID", subc="", ktype="existingvarlist", islist=False),
        Template("HELP", subc="", ktype="bool"),
        
        Template("NAME", subc="DATASET", var="dsname", ktype="varname"),
        Template("ALPHA", subc="LEVEL",  ktype="float", vallist=(.0000000001, .99999999999)),
    ])
    

   # A HELP subcommand overrides all else
    if "HELP" in args:
        print(helptext)
    else:
        try:
            # parse and execute the command
            oobj.parsecmd(args, vardict = spssaux.VariableDict())
            ###print oobj.parsedparams
            dopropor(**oobj.parsedparams)
        except:
            # Exception messages are printed here, but the exception is not propagated, and tracebacks are suppressed,
            # because as an Extension command, the Python handling should be suppressed (unless debug mode)
            if debug:
                raise
            else:
                print(sys.exc_info()[1])
                sys.exc_clear()
    
def dopropor(num=None, denom=None, id=None, dsname="*", alpha=.05, adjust='bonferroni'):
    
    if num is None or denom is None:
        raise ValueError("Error: NUM and DENOM keywords are required")
    if spss.PyInvokeSpss.IsUTF8mode():
        unistr = str
    else:
        unistr = str

    
    currentds = spss.ActiveDataset()
    if currentds=="*":
        currentds = "S"+ str(random.uniform(0, 1))
        spss.Submit("DATASET NAME %s" % currentds)
        dsnamed = True
    else:
        dsnamed = False

    numvec, denomvec, idvec = getvalues(num, denom, id, dsname)
    # clean data, discard missing
    droplist = []
    for i in range(len(numvec)):
        droplist.append(numvec[i] is not None and  denomvec[i] is not None)  #missing data
        if (droplist[i] and (numvec[i] > denomvec[i] or denomvec[i] <= 0)):
            raise ValueError("Error: NUM value greater than DENOM value or zero denominator: %s, %s" % (numvec[i], denomvec[i]))
    for lis in numvec, denomvec, idvec:
        lis = [x for f,x in zip(droplist, lis) if f]  #prune missing values
    if len(numvec) == 0:
        raise ValueError("Error: No valid proportions were found to analyze")

    alphalow = alpha/2
    alphahigh = 1-alphalow
    dotest = len(numvec) > 1
    try:
        spss.StartDataStep()  #TODO: pending transformations
    except:
        spss.Submit("EXECUTE")
        spss.StartDataStep()

    
    # calculate ci's via SPSS IDFs
    
    ds = spss.Dataset(name=None)
    spss.SetActive(ds)
    ds.varlist.append("p",0)
    ds.varlist.append("num",0)
    ds.varlist.append("denom",0)
    
    p0 = numvec[0]/denomvec[0]
    sdvec = []
    for i in range(len(numvec)):
        p1 = numvec[i] / denomvec[i]
        sdvec.append(sqrt(p0*(1-p0)/denomvec[0] + p1*(1-p1)/denomvec[i]))
        #p = (numvec[i] + numvec[0]) / (denomvec[i] + denomvec[0])
        #z = (p1 - p0)/sqrt(p * (1 - p)*(1/denomvec[0] + 1/denomvec[i]))
        
        ds.cases.append([p1, numvec[i], denomvec[i]])
    spss.EndDataStep()
    
    cmd =r"""COMPUTE PLOWBI = IDF.BETA(%(alphalow)s, num + .5, denom-num + .5).
    COMPUTE PHIGHBI = IDF.BETA(%(alphahigh)s, num + .5,  denom - num + .5).
    DO IF num > 0.
    COMPUTE PLOWPOIS = (IDF.CHISQ(%(alphalow)s, 2*num)/2)/denom.
    ELSE.
    COMPUTE PLOWPOIS = 0.
    END IF.
    COMPUTE PHIGHPOIS = (IDF.CHISQ(%(alphahigh)s, 2*(num+1))/2) / denom.
    COMPUTE ZTAIL = IDF.NORMAL(%(alphahigh)s, 0,1).
    EXECUTE."""\
    % {"alphalow": alphalow, "alphahigh": alphahigh}
    
    spss.Submit(cmd)
    plowbi = []
    phighbi = []
    plowpois = []
    phighpois = []
    spss.StartDataStep()
    ds = spss.Dataset(name="*")
    for case in ds.cases:
        i = 3
        for v in plowbi, phighbi, plowpois, phighpois:
            v.append(case[i])
            i += 1
    zalpha2 = case[-1]
    try:
        closeafter = False
        spss.SetActive(spss.Dataset(name=currentds))
    except:
        closeafter = True
    ds.close()
    spss.EndDataStep()

    from spss import CellText
    spss.StartProcedure("Proportions")
    table = spss.BasePivotTable("Proportion Confidence Intervals", "Proportions")
    titlefootnote = "Alpha = %.3f" % alpha
    if 0. in numvec:
        titlefootnote += " (One-sided %.3f when p = 0)" % (alpha/2.)
    table.TitleFootnotes(titlefootnote)
    rowdim = table.Append(spss.Dimension.Place.row, "Proportions")
    coldim =  table.Append(spss.Dimension.Place.column, "Statistics")
    cols = ["p", "Binomial\nLower CI", "Binomial\nUpper CI", "Poisson\nLower CI", "Poisson\nUpper CI", "Difference\nfrom p0", 
        "Difference from p0\nLower CI", "Difference from p0\nUpper CI"]
    table.SetCategories(coldim, [CellText.String(v) for v in cols])
    idvec = [not v is None and unistr(v) or unistr(i+1) for i,v in enumerate(idvec)]
    table.SetCategories(rowdim, [CellText.String(v) for v in idvec])
    for i in range(len(numvec)):
        p1 = numvec[i]/denomvec[i]
        if i > 0:
            zdifflow = p1 - p0 - sdvec[i] * zalpha2
            zdiffhigh = p1 - p0 + sdvec[i] * zalpha2
        else:
            zdifflow = zdiffhigh = 0.
        table.SetCellsByRow(CellText.String(idvec[i]), [CellText.Number(v) for v in (numvec[i]/denomvec[i], plowbi[i], phighbi[i],
            plowpois[i], phighpois[i], p1-p0, zdifflow, zdiffhigh)])
        if i == 0:
            table[(CellText.String(idvec[0]), CellText.String(cols[-3]))] = CellText.String("-")
            table[(CellText.String(idvec[0]), CellText.String(cols[-2]))] = CellText.String("-")
            table[(CellText.String(idvec[0]), CellText.String(cols[-1]))] = CellText.String("-")
    spss.EndProcedure()
    if closeafter:
        spss.Submit(r"""NEW FILE.
        DATASET NAME %s.""" %"S"+ str(random.uniform(0, 1)))
    
    
    
    
def getvalues(num, denom, id, dsname):
    """return vectors of num.  denom, and id values from constants in syntax or variable values"""
    
    if isname(num[0]) or isname(denom[0]) or isname(id):
        spss.StartDataStep()
        ds = spss.Dataset(dsname)
    else:
        ds = None
    id = [id]
    try:
        vallist = []
        if ds:
            vl = [v.name.lower() for v in ds.varlist]  # variables in the dataset
        for v in num, denom, id:
            try:
                vallist.append([float(val) for val in v])
            except:   #variable name argument or None
                if v[0] is None:  # can only happen with id variable
                    vallist.append([None])   # null label in case no id variable
                else:
                    if len(v) > 1:
                        raise ValueError("Error: Only one variable may be named on each of NUM, DENOM, and ID, and a variable may not be combined with a value: " + " ".join(v))
                    try:
                        vindex = vl.index(v[0].lower())
                        vallist.append([val[vindex] for val in ds.cases])
                    except:
                        raise ValueError("Error: An undefined variable name was specified in NUM, DENOM, or ID: " + " ".join(v))
    finally:
        spss.EndDataStep()
        
    # check and fix value list lengths
    maxlen = max([len(vl) for vl in vallist])
    
    for i in range(len(vallist)):
        if len(vallist[i]) == 1:
            vallist[i] = maxlen * vallist[i]
        if len(vallist[i]) != maxlen:
            raise ValueError("Error: NUM, DENOM and optional ID do not all have the same number of items")
    return vallist
        
def isname(v):
    if v is None:
        return False
    try:
        float(v)
        return False
    except:
        return True