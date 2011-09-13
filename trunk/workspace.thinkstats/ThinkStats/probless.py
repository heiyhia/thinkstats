

def CdfProbLess(cdf1, cdf2):
    """Probability that a value from cdf1 is less than one from cdf2.
    For continuous distributions F and G, the chance that a sample
    from F is less than a sample from G is
    Integral(x): pdf_F(x) * (1 - cdf_G(x))
    This function computes an approximation of this Integral using
    discrete CDFs.

    Args:
        cdf1: CDF object
        cdf2: CDF object

    Returns:
        float probability
    """
    total = 0.0
    i = 0
    j = 0
    x = float('-Inf')
    while True:
        # sweep through cdf1 and compute p, the marginal prob of v2
        unused_x1, p1 = cdf1.data[i]
        x2, p2 = cdf1.data[i+1]
        p = p2 - p1
    
        # incr through cdf2 to find Prob{x < x2}
        while x <= x2:
            x, y = cdf2.data[j]
            if j == len(cdf2.data)-1:
                break
            else:
                j += 1
    
        # add up the integral
        total += p * (1 - y)
        i += 1
        if i == len(cdf1.data)-1:
            break

    return total
