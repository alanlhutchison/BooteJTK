import numpy as np
from scipy.stats import kendalltau as kt
from scipy.stats import circmean as sscircmean
from scipy.stats import circstd as sscircstd

def get_stat_probs(dorder,new_header,periods,phases,widths,int size):
    cdef double period,phase,width
    cdef double tau,p
    cdef double m_tau,s_tau,m_per,s_per,m_ph,s_ph,m_na,s_na
    cdef double[6] out1
    cdef double[2] out2
    
    waveform = 'cosine'
    d_taugene,d_pergene,d_phgene,d_nagene = {},{},{},{}
    totals = np.array([complex(0,0),complex(0,0),complex(0,0),complex(0,0),complex(0,0)])

    ### kkey is a sequence of ranks to have eJTK performed on them
    ### dorder[kkey] is a probability of observing that sequence in the bootstraps
    rs = []
    for kkey in dorder:
        res = []
        for period in periods:
            for phase in phases:
                for width in widths:
                    serie = kkey
                    reference = generate_base_reference(new_header,waveform,period,phase,width)
                    tau,p = kt(reference,serie)#generate_mod_series(reference,serie)
                    p = p/2.0
                    tau = farctanh(tau)
                    maxloc = new_header[serie.index(max(serie))]
                    minloc = new_header[serie.index(min(serie))]  
                    
                    res.append([tau,p,period,phase,periodic(phase+width),maxloc,minloc])
        r = pick_best_match(res)
        r = list(r)

        d_taugene.setdefault(r[0],0)
        d_taugene[r[0]]+=dorder[kkey]
        
        d_pergene.setdefault(r[2],0)
        d_pergene[r[2]]+=dorder[kkey]

        d_phgene.setdefault(r[3],0)
        d_phgene[r[3]]+=dorder[kkey]

        d_nagene.setdefault(r[4],0)
        d_nagene[r[4]]+=dorder[kkey]

        for _ in xrange(int(np.round(size*dorder[kkey]))):
            rs.append(r)
    rs = np.array(rs)
    #print rs

    m_tau = np.mean(rs[:,0])
    s_tau = np.std(rs[:,0])
    m_per = np.mean(rs[:,2])
    s_per = np.std(rs[:,2])
    m_ph = sscircmean(rs[:,3],high=24,low=0)
    s_ph = sscircstd(rs[:,3],high=24,low=0)    
    m_na = sscircmean(rs[:,4],high=24,low=0)
    s_na = sscircstd(rs[:,4],high=24,low=0)    
    out1,out2 = [m_per,s_per,m_ph,s_ph,m_na,s_na],[m_tau,s_tau]
    return out1,out2,d_taugene,d_pergene,d_phgene,d_nagene

def generate_base_reference(header,waveform="cosine",double period=24.,double phase=0.,double width=12.):
    """
    This will generate a waveform with a given phase and period based on the header, 
    """
    cdef double coef,w,tpoint
    #print header,phase
    ZTs = np.array(header,dtype=float)
    coef = 2.0 * np.pi / period
    w = width * coef
    tpoints = (ZTs - phase) *coef

    def trough(double x,double w):
        x = x % (2*np.pi)
        w = w % (2*np.pi)
        if x <= w:
            y = 1 + -x/w
        elif x > w:
            y = (x-w)/(2*np.pi - w)
        return y

    def cosine(double x,double w):
        cdef double y
        x = x % (2*np.pi)
        w = w % (2*np.pi)
        if x <= w:
            y = np.cos(x/(w/np.pi))
        elif x > w:
            y = np.cos( (x+2.*(np.pi-w))*np.pi/ (2*np.pi - w) )
        return y
    d_f = {'cosine':cosine,'trough':trough}
    reference = [d_f[waveform](tpoint,w) for tpoint in tpoints]
    return reference

def farctanh(double x):
    if x>0.99:
        return np.arctanh(0.99)
    elif x<-0.99:
        return np.arctanh(-0.99)
    else:
        return np.arctanh(x)

def periodic(double x):
    x = float(x)
    while x>12:
        x=x-24.
    while x<=-12:
        x= x+24.
    return x

def pick_best_match(res):
    cdef double maxtau,minphasediff,mindiff
    cdef int ind
    
    taus = [r[0] for r in res]
    maxtau = max(taus)
    tau_mask = np.array([maxtau==tau for tau in taus])
    if np.sum(tau_mask)==1:
        ind = list(tau_mask).index(True)
        return res[ind]

    res = np.array(res)[tau_mask]
    phases = [np.abs(r[3]-r[5]) for r in res]
    minphasediff = min(phases)
    phasemask = np.array([minphasediff==phase for phase in phases])
    if np.sum(phasemask)==1:
        ind = list(phasemask).index(True)
        return res[ind]

    res = np.array(res)[phasemask]
    diffs = [np.abs(r[4]-r[6]) for r in res]
    mindiff = min(diffs)
    diffmask = np.array([mindiff==diff for diff in diffs])
    if np.sum(diffmask)==1:
        ind = list(diffmask).index(True)
        return res[ind]

    ### If we've gotten down here everything has failed
    print 'Ties remain...',res
    return res[np.random.randint(len(res))]
