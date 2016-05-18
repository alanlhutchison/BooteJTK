import numpy as np
#cimport numpy as np
from scipy.stats import kendalltau as kt
from scipy.stats import circmean as sscircmean
from scipy.stats import circstd as sscircstd
#from multiprocessing import Pool


def get_stat_probs(dorder,new_header,triples,dref,int size):
    #periods,phases,widths,int size):
    cdef double period,phase,width,nadir
    cdef double tau,p
    cdef double m_tau,s_tau,m_per,s_per,m_ph,s_ph,m_na,s_na
    cdef double[6] out1
    cdef double[2] out2
    cdef double[2] pair
    #cdef cnp.ndarray pairs = np.zeros((len(phases)*len(widths),2))
    cdef str waveform = 'cosine'
    cdef int i,j
    cdef dict d_taugene,d_pergene,d_phgene,d_nagene 
    d_taugene,d_pergene,d_phgene,d_nagene = {},{},{},{}
    #waveform = 'cosine'

    totals = np.array([complex(0,0),complex(0,0),complex(0,0),complex(0,0),complex(0,0)])

    ### kkey is a sequence of ranks to have eJTK performed on them
    ### dorder[kkey] is a probability of observing that sequence in the bootstraps

    rs = []
    for kkey in dorder:
        lamb_get_matches = lambda triple: get_matches(kkey,triple,dref,new_header)
        res = np.array(map(lamb_get_matches,triples))
        r = pick_best_match(res)
        #r = list(r)
          
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
    if waveform=='cosine':

        def cosine(double x,double w):
            cdef double y
            x = x % (2*np.pi)
            w = w % (2*np.pi)
            if x <= w:
                y = np.cos(x/(w/np.pi))
            elif x > w:
                y = np.cos( (x+2.*(np.pi-w))*np.pi/ (2*np.pi - w) )
            return y
        f_wav = cosine
    elif waveform=='trough':
        def trough(double x,double w):
            cdef double y
            x = x % (2*np.pi)
            w = w % (2*np.pi)
            if x <= w:
                y = 1 + -x/w
            elif x > w:
                y = (x-w)/(2*np.pi - w)
            return y
        f_wav = trough
        
    reference = [f_wav(tpoint,w) for tpoint in tpoints]
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
    #cdef np.ndarray taus,phases,diffs
    #cdef np.ndarray tau_mask,phasemask,diffmask
    cdef int ind
    
    res = np.array(res)
    taus = res[:,0]
    #maxtau = max(taus)
    tau_mask = (max(taus)==taus)
    if np.sum(tau_mask)==1:
        ind = list(tau_mask).index(True)
        return res[ind]

    res = res[tau_mask]
    phases = np.abs(res[:,3]-res[:,5])
    #minphasediff = min(phases)
    phasemask = (min(phases)==phases)
    if np.sum(phasemask)==1:
        ind = list(phasemask).index(True)
        return res[ind]

    res = res[phasemask]
    diffs = np.abs(res[:,4]-res[:,6])
    diffmask = (min(diffs)==diffs)
    if np.sum(diffmask)==1:
        ind = list(diffmask).index(True)
        return res[ind]

    ### If we've gotten down here everything has failed
    print 'Ties remain...',res
    return res[np.random.randint(len(res))]


def get_waveform_list(periods,phases,widths):
    cdef int lper = len(periods)
    cdef int lpha = len(phases)
    cdef int lwid = len(widths)
    #cdef np.ndarray
    triples = np.zeros((lper*lpha*lwid/2,3))
    cdef int i,j
    cdef int period,phase,width,nadir

    
    for i,period in enumerate(periods):
        j = 0
        pairs = [[0,0]]*(lpha*lwid/2)
        for phase in phases:
            for width in widths:
                nadir = (phase+width)%period
                #print
                #print phase,nadir
                pair = [nadir,phase]
                #print pair, 'not in', pairs,pair not in pairs
                if pair not in pairs:
                    #print 'Adding', phase,nadir
                    pairs[j] = [phase,nadir]
                    triples[i*lper+j] = np.array([period,phase,width])
                    j+=1
        #for pair in pairs:
        #    print pair
    triples = np.array(triples,dtype=float)
    #for trip in triples:
    #    print trip
    return triples


def make_references(new_header,triples,waveform='cosine'):#,period,phase,width):
    cdef double period,phase,width
    cdef dict dref ={}
    for triple in triples:
        period,phase,width = triple
        reference = generate_base_reference(new_header,waveform,period,phase,width)
        dref[(period,phase,width)] = reference
    return dref


def get_matches(kkey,triple,d_ref,new_header):
    cdef double period,phase,width,nadir,tau,p
    
    reference = d_ref[tuple(triple)]
    period,phase,width = triple
    nadir = (phase+width)%period
    serie = kkey            
    tau,p = kt(reference,serie)#generate_mod_series(reference,serie)
    p = p/2.0
    tau = farctanh(tau)
    maxloc = new_header[serie.index(max(serie))]
    minloc = new_header[serie.index(min(serie))]
    r =  [tau,p,period,phase,nadir,maxloc,minloc]
    if tau < 0:
        r = [np.abs(tau),p,period,nadir,phase,maxloc,minloc]
    
    return r

