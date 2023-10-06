'''
Useful functions to aide in qaqc functions. 
Modified from https://github.com/rjhd2/HadISD_v3/blob/420ebb3c965b61ad394f15b0f71a37b73985a97d/qc_utils.py
'''

from scipy.optimize import leastsq

def create_bins(data, bin_size=0.25):
    '''Create bins from data covering entire data range'''

    # set up bins
    b_min = np.floor(np.min(data))
    b_max = np.ceil(np.max(data))
    bins = np.arange(b_min - bin_size, b_max + (3. * bin_size), bin_size)

    return bins


def iqr(data, percentile=0.25):
    '''Calculate IQR of data'''

    # sort data
    sorted_data = sorted(data)

    n_data = len(sorted_data)
    quartile = int(round(percentile * n_data))

    return sorted_data[n_data - quartile] - sorted_data[quartile] # iqr


def gaussian(X, p):
    '''
    Guassian line fitting, where:
        p[0]=mean
        p[1]=sigma
        p[2]=normalisation
    '''

    return (p[2]*(np.exp(-((X-p[0])*(X-p[0]) / 2.0*p[1]*p[1]))))


def gaussian_resid(p, Y, X):
    '''Least squared residuals from linear trend'''

    # calculate residuals
    err = ((Y-gaussian(X, p)**2.0))

    return err


def fitted_gaussian(x, y, norm, mu=False, sig=False):
    '''Creates a fitted Gaussian'''

    if not mu:
        mu = np.mean(x)
    if not sig:
        sig = np.std(x)

    p0 = np.array([mu, sig, norm])

    fit, success = leastsq(gaussian_resid, p0, args = (y,x), maxfev=10000, full_output=False)

    return fit