"""Interface to rpy2.glm

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import rpy2.robjects as robjects
r = robjects.r


def linear_model(model, print_flag=True):
    """Submits model to r.lm and returns the result."""
    model = r(model)
    res = r.lm(model)
    if print_flag:
        print_summary(res)
    return res


def run_model(model, family=robjects.r.binomial(), print_flag=True):
    """Submits model to r.lm and returns the result."""
    model = r(model)
    res = r.glm(model, family=family)
    if print_flag:
        print_summary(res)
    return res


def print_summary(res):
    """Prints results from r.lm (just the parts we want)."""
    flag = False
    lines = r.summary(res)
    lines = str(lines)

    for line in lines.split('\n'):
        # skip everything until we get to coefficients
        if line.startswith('Coefficients'):
            flag = True
        if line.startswith('Signif'):
            continue
        if flag:
            print line
    print


def get_coeffs(res):
    """Gets just the lines that contain the estimates.

    res: R glm result object

    Returns: map from coefficient name to (estimate, error, z-value) tuple
    """
    flag = False
    lines = r.summary(res)
    lines = str(lines)

    res = {}
    for line in lines.split('\n'):
        line = line.strip()
        if line.startswith('---'):
            break
        if flag:
            t = line.split()
            var = t[0]
            est = float(t[1])
            error = float(t[2])
            z = float(t[3])
            res[var] = est, error, z
        if line.startswith('Estimate'):
            flag = True
    return res

def inject_col_dict(col_dict, prefix=''):
    """Copies data columns into the R global environment.
    
    col_dict: map from attribute name to column of data
    prefix: string prepended to the attribute names
    """
    for name, col in col_dict.iteritems():
        robjects.globalEnv[prefix+name] = robjects.FloatVector(col)


