#!/usr/bin/env python
"""
Created on Nov 1, 2015
@author: Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron R. Dinner Group, University of Chicago

This script is a bootstrapped expansion of the eJTK method described in

Hutchison AL, Maienschein-Cline M, and Chiang AH et al. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data, PLoS Computational Biology 2015 11(3): e 1004094. doi:10.1371/journal.pcbi.1004094

This script bootstraps time series and provides phase and tau distributions from those bootstraps to allow for measurement of the variance on phase and tau estimates.


Please use ./BooteJTK -h to see the help screen for further instructions on running this script.

"""
VERSION="0.1"

import cmath
import scipy.stats as ss
import numpy as np
from scipy.stats import kendalltau as kt
from scipy.stats import multivariate_normal as mn
from scipy.stats import rankdata
from scipy.stats import norm
from scipy.special import polygamma

import pickle
from operator import itemgetter
import sys
import argparse
import time
import os.path

def main(args):
    fn = args.filename
    prefix = args.prefix
    fn_waveform = args.waveform
    fn_period = args.period
    fn_phase = args.phase
    fn_width = args.width
    fn_out = args.output
    fn_out_pkl = args.pickle # This is the output file which could be read in early
    fn_list = args.id_list # This is the the list of ids to go through
    fn_null_list = args.null_list # These are geneIDs to be used to estimate the SD
    size = int(args.size)
    ### If no list file set id_list to empty
    if fn_list!='DEFAULT':
        id_list = read_in_list(fn_list)
    else:
        id_list = []

    ### If no pkl out file, modify the option variable
    if os.path.isfile(fn_out_pkl):
        opt = 'premade'
    else:
        opt = ''

    ### The number of bootstraps to be done
    # This is now a input parameter
    #size = 100

    ### If not given a new name, name fn_out after the fn file
    if fn_out == "DEFAULT":
        if ".txt" in fn:
            fn_out=fn.replace(".txt","_"+prefix+"_bootejtk.txt")
        else:
            fn_out = fn+"_" +prefix + "_bootejtk.txt"
        
    ### If not given a new name, name fn_out_pkl based on the fn file
    if fn_out_pkl == 'DEFAULT':
        fn_out_pkl = fn.replace('.txt','_boot{}_order_probs.pkl'.format(int(np.log10(size))))

    ### Name vars file based on pkl out file
    fn_out_pkl_vars = fn_out_pkl.replace('.pkl','_vars.pkl')

 
    assert fn!='DEFAULT' or fn_out_pkl!='DEFAULT'   
    ### If we already have the PKL file, we just need a place to put the header information
    if fn=='DEFAULT' and fn_out_pkl!='DEFAULT':
        #original = False # This is currently unused
        series = [[key,0] for key in d_data_master.keys()]
        fn_out= fn_out_pkl.replace('.pkl','_{0}-bootejtk.txt'.format(int(np.log10(size))))
        if fn_out_pkl[:-4]=='2.pkl':
            new_header = [0,4,8,12,16,20,0,4,8,12,16,20]
        else:
            new_header = [0,4,8,12,16,20]

    ### If we have the initial data we can get it 
    elif fn!='DEFAULT':
        header,data = read_in(fn)
        header,series = organize_data(header,data)
        d_data_master,new_header = get_data(header,data)
        print new_header

    
    if 'premade' not in opt:
        null_ids = read_in_list(fn_null_list)
        D_null = get_series_data(data,null_ids)
        d_data_master = eBayes(d_data_master,D_null)
    elif 'premade' in opt:
        d_data_master,d_order_probs = open_pickle_append2(fn_out_pkl)
        # RUN THIS IN THE EM CASE
        #d_data_master,d_order_probs = pickle.load(open(fn_out_pkl,'rb'))#open_pickle_append2(fn_out_pkl)        

    add_on = 1
    while os.path.isfile(fn_out):
        print fn_out, "already exists, take evasive action!!!"
        if add_on==1:
            origstr = '.txt'
        else:
            origstr = '_{0}.txt'.format(add_on-1)
        fn_out = fn_out.replace(origstr,'_{0}.txt'.format(add_on))
        add_on = add_on + 1

    ### This is for the Hogen case ###
    #new_header = [0,2,4,6,8,10,12,14,16,18,20,22]

    
    #for key in d_data_master:
    #d_order_probs = get_order_prob2(d_data_master,size)
    #pickle.dump([d_data_master,d_order_probs], open(fn_out_pkl,'wb'))
    ### Read in lists of search parameters
    waveforms = read_in_list(fn_waveform)
    periods = read_in_list(fn_period)
    phases = read_in_list(fn_phase)
    widths = read_in_list(fn_width)
    
    # If BooteJTK2
    #new_header = list(new_header)+list(new_header)


    
    print 'Pickled Orders'
    waveform = []
    Ps = []
    print 'Begin eJTK'
    #with open(fn_out,'w') as g:
    d_tau = {}
    d_ph = {}
    #if 1:
    done = []
    remaining = []
    """Time Limit for Code to Run (in hours)"""
    time_limit = 34
    time_limit_sec = float(60*60*time_limit)
    g = open(fn_out,'a')
    g.write("ID\tWaveform\tPeriodMean\tPeriodStdDev\tPhaseMean\tPhaseStdDev\tNadirMean\tNadirStdDev\tMean\tStd_Dev\tMax\tMin\tMax_Amp\tFC\tIQR_FC\tNumBoots\tTauMean\tTauStdDev\n")
    g.close()
    time_original = time.time()
    for serie in series:
        ### If we have an ID list, we only want to deal with data from it.
        if id_list!=[]:
            if serie[0] not in id_list:
                #print 'Passing'
                pass
        ### We have time limits here so we don't blow out the Midway allocation
        time_diff = time.time() - time_original
        if time_diff < time_limit_sec:
            if fn=='DEFAULT':
                mmax,mmin,MAX_AMP = np.nan,np.nan,np.nan
                sIQR_FC = np.nan
                smean = np.nan
                sstd = np.nan
                sFC = np.nan
            else:
                mmax,mmin,MAX_AMP=series_char(serie)
                sIQR_FC=IQR_FC(serie)
                smean = series_mean(serie)
                sstd = series_std(serie)
                sFC = FC(serie)

            #local_ps = []

            s_stats = [smean,sstd,mmax,mmin,MAX_AMP,sFC,sIQR_FC,size]

            ### The first line of the data is the geneID
            geneID = serie[0] 
            ### Need to make this file if it doesn't exist already
            if 'premade' not in opt:
                d_data_sub = {geneID:d_data_master[geneID]}
                d_order_probs = get_order_prob(d_data_sub,size)
                f1 = open(fn_out_pkl,'ab')
                pickle.dump([d_data_sub,d_order_probs],f1)
                f1.close()

            #totals = np.array([complex(0,0),complex(0,0),complex(0,0),complex(0,0),complex(0,0)])
            #d_tau[geneID] = {}
            #d_ph[geneID] = {}
            if geneID in d_order_probs:
                #print geneID,'in d_order_probs'
                out1,out2,d_taugene,d_phgene = get_stat_probs(d_order_probs[geneID],new_header,periods,phases,widths)
                print out1,out2
                out_line = [geneID,waveform]+out1+s_stats+out2

                out_line = [str(l) for l in out_line]
                g = open(fn_out,'a')
                g.write("\t".join(out_line)+"\n")
                g.close()
                print geneID
                done.append(geneID)
                sys.stdout.flush()
        
                pickle.dump([{geneID:d_taugene},{geneID:d_phgene}],open(fn_out_pkl_vars,'ab'))
            #else:
                #pprint 'Gene not in pkl',geneID
        else:
            remaining.append(serie[0])
            #print 'Time is up'
    if len(remaining)>0:
        fn_remaining = fn_out.replace('.txt','_remaining_list.txt')
        with open(fn_remaining,'w') as g:
            for r in remaining:
                g.write(r+'\n')
                
def get_stat_probs(dorder,new_header,periods,phases,widths):
    RealKen = KendallTauP()
    waveform = 'cosine'
    d_taugene,d_phgene = {},{}
    totals = np.array([complex(0,0),complex(0,0),complex(0,0),complex(0,0),complex(0,0)])
    for kkey in dorder:
        res = []
        for period in np.array(periods,dtype=float):
            for phase in np.array(phases,dtype=float):
                for width in np.array(widths,float):
                    serie = kkey
                    reference = generate_base_reference(new_header,waveform,period,phase,width)
                    tau,p = generate_mod_series(reference,serie,RealKen)
                    res.append([tau,p,period,phase,periodic(phase+width)])
        r = sorted(res)[-1]
        if r[0] not in d_taugene:
            d_taugene[r[0]] = 0
        d_taugene[r[0]]+=dorder[kkey]

        if r[3] not in d_phgene:
            d_phgene[r[3]] = 0
        d_phgene[r[3]]+=dorder[kkey]

        r[3] = r[3] * np.pi /12.
        r[4] = r[4] * np.pi /12.
        r=np.array([r[0],r[1],r[2],complex(np.cos(r[3]),np.sin(r[3])),complex(np.cos(r[4]),np.sin(r[4]))])
        sub = r*dorder[kkey]
        totals+= sub
    m = totals
    m = np.array([m[0].real,m[1].real,m[2].real,cmath.phase(m[3]),cmath.phase(m[4])],dtype=float)*np.array([1,1,1,12/np.pi,12/np.pi])
    cos_sum1,sin_sum1,cos_sum2,sin_sum2,tau_var,p_var,per_var = 0.,0.,0.,0.,0.,0.,0.            
    for kkey in dorder:
        res = []
        for period in periods:
            for phase in phases:
                for width in widths:
                    serie = kkey
                    reference = generate_base_reference(new_header,waveform,period,phase,width)
                    tau,p = generate_mod_series(reference,serie,RealKen)
                    res.append([tau,p,period,phase,periodic(phase+width)])
        r = np.array(sorted(res)[-1],dtype=float)
        r[3] = r[3] * np.pi /12.
        r[4] = r[4] * np.pi /12.
        tau_var+=(r[0]-float(m[0]))*(r[0]-float(m[0]))*dorder[kkey]
        p_var+=(r[1]-float(m[1]))*(r[1]-float(m[1]))*dorder[kkey]
        per_var+=(r[2]-float(m[2]))*(r[2]-float(m[2]))*dorder[kkey]
        cos_sum1+=(np.cos(r[3])*dorder[kkey])
        sin_sum1+=(np.sin(r[3])*dorder[kkey])
        cos_sum2+=(np.cos(r[4])*dorder[kkey])
        sin_sum2+=(np.sin(r[4])*dorder[kkey])   

    ph_var=(1.-np.sqrt(cos_sum1**2 + sin_sum1**2))*24.
    nad_var=(1.-np.sqrt(cos_sum2**2 + sin_sum2**2))*24.
    tau_mean = m[0]
    p_mean = m[1]
    per_mean = m[2]
    ph_mean = m[3]
    nad_mean = m[4]

    out1,out2 = [per_mean,np.sqrt(per_var),ph_mean,np.sqrt(ph_var),nad_mean,np.sqrt(nad_var)],[tau_mean,np.sqrt(tau_var)]
    return out1,out2,d_taugene,d_phgene

            
def rename(x):
    return x
def open_pickle_append3(fn):
    d = {}
    objs =[]
    f = open(fn,'rb')
    while 1:
        try:
            objs.append(pickle.load(f))
        except EOFError:
            break
    f.close()
    d1 = {}
    d2 = {}
    d3 = {}
    for aa in objs:
        new = rename(aa[0].keys()[0])
        d1[new]=aa[0][aa[0].keys()[0]]
        d2[new]=aa[1][aa[1].keys()[0]]
        d3[new]=aa[2][aa[2].keys()[0]]

def open_pickle_append2(fn):
    d = {}
    objs =[]
    f = open(fn,'rb')
    while 1:
        try:
            objs.append(pickle.load(f))
        except EOFError:
            break
    f.close()
    d1 = {}
    d2 = {}
    for aa in objs:
        new = rename(aa[0].keys()[0])
        d1[new]=aa[0][aa[0].keys()[0]]
        d2[new]=aa[1][aa[1].keys()[0]]
    return [d1,d2]


def append_out(fn_out,line):
    line = [str(l) for l in line]
    with open(fn_out,'a') as g:
        g.write("\t".join(line)+"\n")

def write_out(fn_out,output):
    with open(fn_out,'w') as g:
        for line in output:
            g.write(str(line)+"\n")

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_in_list(fn):
    with open(fn,'r') as f:
        lines = f.read().splitlines()
    return lines
        
def read_in(fn):
    """Read in data to header and data"""
    with open(fn,'r') as f:
        data=[]
        start_right=0
        for line in f:
            words = line.strip().split()
            words = [word.strip() for word in words]
            if words[0] == "#":
                start_right = 1
                header = words[1:]
            else:
                if start_right == 0:
                    print "Please enter file with header starting with #"
                elif start_right == 1:
                    data.append(words)
    return header, data

def organize_data(header,data):
    """
    Organize list of lists from such that genes with similar time-series holes match (for null distribution calc)
    Return a header ['#','ZTX','ZTY'...] and a list of lists [ lists with similar holes (identical null distribution) , [],[],[]] 
    """
    L = data

    for i in xrange(1,len(header)):
        L=sorted(L, key=itemgetter(i))
    return header,L

def generate_base_reference(header,waveform="cosine",period=24,phase=0,width=12):
    """
    This will generate a waveform with a given phase and period based on the header, 
    """
    #print header,phase
    tpoints = []
    ZTs = header
    coef = 2.0 * np.pi / float(period)
    w = float(width) * coef
    for ZT in ZTs:
        z = ZT
        tpoints.append( (float(z)-float(phase) ) * coef)


    def trough(x,w):
        x = x % (2*np.pi)
        w = w % (2*np.pi)
        if x <= w:
            y = 1 + -x/w
        elif x > w:
            y = (x-w)/(2*np.pi - w)
        return y

    def cosine(x,w):
        x = x % (2*np.pi)
        w = w % (2*np.pi)
        if x <= w:
            y = np.cos(x/(w/np.pi))
        elif x > w:
            y = np.cos( (x+2.*(np.pi-w))*np.pi/ (2*np.pi - w) )
        return y


    if waveform == "cosine":
        reference=[cosine(tpoint,w) for tpoint in tpoints]
    elif waveform == "trough":
        reference=[trough(tpoint,w) for tpoint in tpoints]
    return reference



def IQR_FC(series):
    qlo = __score_at_percentile__(series, 25)
    qhi = __score_at_percentile__(series, 75)
    if (qlo=="NA" or qhi=="NA"):
        return "NA"
    elif (qhi==0):
        return 0
    elif ( qlo==0):
        return "NA"
    else:
        iqr = qhi/qlo
        return iqr

def FC(series):
    series=[float(s) if s!="NA" else 0 for s in series[1:] if s!="NA"  ]
    if series!=[]:
        mmax = max(series)
        mmin = min(series)
        if mmin==0:
            sFC = -10000
        else:
            sFC = mmax / mmin
    else:
        sFC = "NA"
    return sFC


def series_char(series):
    """Uses interquartile range to estimate amplitude of a time series."""
    series=[float(s) for s in series[1:] if s!="NA"]
    if series!=[]:
        mmax = max(series)
        mmin = min(series)
        diff=mmax-mmin
    else:
        mmax = "NA"
        mmin = "NA"
        diff = "NA"
    return mmax,mmin,diff


def series_mean(series):
    """Finds the mean of a timeseries"""
    series = [float(s) for s in series[1:] if s!="NA"]
    return np.mean(series)

def series_std(series):
    """Finds the std dev of a timeseries"""
    series = [float(s) for s in series[1:] if s!="NA"]
    return np.std(series)

def __score_at_percentile__(ser, per):
    ser = [float(se) for se in ser[1:] if se!="NA"]
    if len(ser)<5:
        score ="NA"
        return score
    else: 
        ser = np.sort(ser)
        i = (per/100. * len(ser))
        if (i % 1 == 0):
            score = ser[i]
        else:
            interpolate = lambda a,b,frac: a + (b - a)*frac
            score = interpolate(ser[int(i)], ser[int(i) + 1], i % 1)
        return float(score)

def generate_mod_series(reference,series,RealKen):
    """
    Takes the series from generate_base_null, takes the list from data, and makes a null
    for each gene in data or uses the one previously calculated.
    Then it runs Kendall's Tau on the exp. series against the null
    """
    values = series
    binary = np.array([1.0 if value!="NA" else np.nan for value in values])
    reference = np.array(reference)
    temp = reference*binary
    mod_reference = [value for value in temp if not np.isnan(value)]
    mod_values = [float(value) for value in values if value!='NA']

    if len(mod_values) < 3:
        tau,p = np.nan,np.nan
    elif mod_values.count(np.nan) == len(mod_values):
        tau,p = np.nan,np.nan
    elif mod_values.count(0) == len(mod_values):
        tau,p = np.nan,np.nan
    else:
        tau,p=kt(mod_values,mod_reference)
        if not np.isnan(tau):
            if len(mod_values) < 150:
                pk = RealKen.pval(tau,len(mod_values))
                if pk is not None:
                    p=pk
            else:
                p = p / 2.0
                if tau < 0:
                    p = 1-p
    return tau,p

##################################
### HERE WE INSERT BOOT FUNCTIONS
##################################

def get_data(header,data):
    new_h = [float(h[2:])%24 for h in header]
    print new_h
    seen = []
    dref = {}
    for i,h in enumerate(new_h):
        if h not in seen:
            seen.append(h)
            dref[h]=[i]
        else:
            dref[h].append(i)
    d_data = {}
    for dat in data:
        name=dat[0]
        series = [float(da) for da in dat[1:]]
        out = [[],[],[]]
        for i,s in enumerate(seen):
            points = [series[idx] for idx in dref[s]]
            N = len(points)
            m = np.mean(points)
            std = np.std(points)
            out[0].append(m)
            out[1].append(std)
            out[2].append(N)
        #print name,seen,out
        d_data[name]=out
    #print d_data.keys()
    return d_data,seen


def get_series_data(data,id_list):
    d_data = {}
    for dat in data:
        name = dat[0]
        if name in id_list:
            series = [float(da) for da in dat[1:]]
            d_data[name] = [np.mean(series),np.std(series),len(series)]
    return d_data

def eBayes(d_data,D_null={}):
    """
    This is based on Smyth 2004 Stat. App. in Gen. and Mol. Biol. 3(1)3
    It generates a dictionary with eBayes adjusted variance values.
    Ns have been set to 1 to not complicate downstream analyses.
    """
    def get_d0_s0(d_data,D_null):
        digamma = lambda x: polygamma(0,x)
        trigamma = lambda x: polygamma(1,x)
        tetragamma = lambda x: polygamma(2,x)

        def solve_trigamma(x):
            """To solve trigamma(y)=x, x>0"""
            y0 = 0.5 + 1./x
            d =1000000.
            y = y0
            while -d/y >1e-8:
                d = trigamma(y)*(1-trigamma(y)/x)/tetragamma(y)
                y = y+d
            return y

        if D_null!={}:
            dg = np.hstack(np.array(D_null.values())[:,2])
            s  = np.hstack(np.array(D_null.values())[:,1])            
        else:
            dg = np.hstack(np.array(d_data.values())[:,2])
            s  = np.hstack(np.array(d_data.values())[:,1])

        G = len(d_data)
        s2 = np.array([ss for ss in s if ss!=0])
        dg2 = np.array([dg[i] for i,ss in enumerate(s) if ss!=0])
        G = len(dg2)
        z = 2.*np.log(s2)

        e = z - digamma(dg2/2) + np.log(dg2/2)
        emean = np.mean(e)
        d0 = 2.* solve_trigamma( np.mean( (e-emean)**2 *G/(G-1)-trigamma(dg2/2) )   )
        s0 = np.sqrt(np.exp(emean + digamma(d0/2)- np.log(d0/2.)))

        #print d0,s0
        return d0,s0
    def posterior_s(d0,s0,s,d):
        return np.sqrt( (d0*s0**2 + d*s**2) /(d0+d) )

    d0,s0=get_d0_s0(d_data,D_null)
    
    for key in d_data:
        d_data[key][1]=[posterior_s(d0,s0,d_data[key][1][i],d_data[key][2][i]) for i in xrange(len(d_data[key][1]))  ]
        d_data[key][2]=[1 for i in xrange(len(d_data[key][1]))  ]
    print d0,s0
    return d_data

def get_order_prob(d_data,size):
    d_order_prob = {}
    #print d_data.keys()
    for key in d_data:
        res = d_data[key]
        d_order_prob[key]=dict_of_orders(res[0],res[1],res[2],size)
    return d_order_prob

def get_order_prob2(d_data,size):
    d_order_prob = {}
    #print d_data.keys()
    for key in d_data:
        res = d_data[key]
        res[0] = list(res[0])+list(res[0])
        res[1] = list(res[1])+list(res[1])
        res[2] = list(res[2])+list(res[2])        
        d_order_prob[key]=dict_of_orders(res[0],res[1],res[2],size)
    return d_order_prob



def dict_of_orders(M,SDS,NS,size):
    """
    This produces a dictionary of probabilities
    for the different orders given the data
    """
    index = range(len(M))
    dorder = dict_order_probs(M,SDS,NS,size)
    d = {}
    for idx in dorder:
        d[idx]=dorder[idx]
    SUM = sum(d.values())
    for key in d:
        d[key]=d[key]/SUM
    return d

def dict_order_probs(ms,sds,ns,size=1000000):
    sds = [sds[i]/np.sqrt(ns[i]) for i in xrange(len(sds))]
    cov=np.diag(sds)
    A = mn(ms,cov)
    d = {}
    for s in A.rvs(size):
        r=tuple(map(int,rankdata(s)))
        if r not in d:
            d[r] = 1./size
        else:
            d[r] += 1./size
    return d

def periodic(x):
    x = float(x)
    while x>12:
        x=x-24.
    while x<=-12:
        x= x+24.
    return x




def __create_parser__():
    p = argparse.ArgumentParser(
        description="Python script for running empirical JTK_CYCLE with asymmetry search as described in Hutchison, Maienschein-Cline, and Chiang et al. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data, PLoS Computational Biology 2015 11(3): e1004094. This script was written by Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron R. Dinner Group, University of Chicago.",
        epilog="Please contact the correpsonding author if you have any questions.",
        version=VERSION
        )

                   
    #p.add_argument("-t", "--test",
    #               action='store_true',
    #               default=False,
    #               help="run the Python unittest testing suite")
    p.add_argument("-o", "--output",
                   dest="output",
                   action='store',
                   metavar="filename string",
                   type=str,
                   default = "DEFAULT",
                   help="You want to output something. If you leave this blank, _jtkout.txt will be appended to your filename")

    p.add_argument("-k","--pickle",
                          dest="pickle",
                          metavar="filename string",
                          type=str,
                          action='store',
                          default="DEFAULT",
                          help='Should be a file with phases you wish to search for listed in a single column separated by newlines.\
                          Provided file is "period_24.txt"')


    p.add_argument("-l","--list",
                          dest="id_list",
                          metavar="filename string",
                          type=str,
                          action='store',
                          default="DEFAULT",
                          help='A filename of the ids to be run in this time series. If time is running out on your job, this will be compared to the ids that have already been completed and a file will be created stating what ids remain to be analyzed.')


    p.add_argument("-n","--null",
                          dest="null_list",
                          metavar="filename string",
                          type=str,
                          action='store',
                          default="DEFAULT",
                          help='A filename of the ids upon which to calculate the standard deviation. These ids are non-cycling, so the standard deviation can be taken across the entire time series. Using this argument is useless if the bootstraps have already been performed.')
    

    analysis = p.add_argument_group(title="JTK_CYCLE analysis options")

    analysis.add_argument("-f", "--filename",
                   dest="filename",
                   action='store',
                   metavar="filename string",
                   type=str,
                   help='This is the filename of the data series you wish to analyze.\
                   The data should be tab-spaced. The first row should contain a # sign followed by the time points with either CT or ZT preceding the time point (such as ZT0 or ZT4). Longer or shorter prefixes will not work. The following rows should contain the gene/series ID followed by the values for every time point. Where values are not available NA should be put in it\'s place.')


    analysis.add_argument('-x',"--prefix",
                          dest="prefix",
                          type=str,
                          metavar="string",
                          action='store',
                          default="",
                          help="string to be inserted in the output filename for this run")


    analysis.add_argument('-z',"--size",
                          dest="size",
                          type=str,
                          metavar="int",
                          action='store',
                          default="",
                          help="Number of bootstraps to be performed")
    

    analysis.add_argument('-w',"--waveform",
                          dest="waveform",
                          type=str,
                          metavar="filename string",
                          action='store',
                          default="cosine",
                          #choices=["waveform_cosine.txt","waveform_rampup.txt","waveform_rampdown.txt","waveform_step.txt","waveform_impulse.txt","waveform_trough.txt"],
                          help='Should be a file with waveforms  you wish to search for listed in a single column separated by newlines.\
                          Options include cosine (dflt), trough')

    analysis.add_argument("--width", "-a", "--asymmetry",
                          dest="width",
                          type=str,
                          metavar="filename string",
                          action='store',
                          default="widths_02-22.txt",
                          #choices=["widths_02-22.txt","widths_04-20_by4.txt","widths_04-12-20.txt","widths_08-16.txt","width_12.txt"]
                          help='Should be a file with asymmetries (widths) you wish to search for listed in a single column separated by newlines.\
                          Provided files include files like "widths_02-22.txt","widths_04-20_by4.txt","widths_04-12-20.txt","widths_08-16.txt","width_12.txt"\nasymmetries=widths')
    analysis.add_argument('-s', "-ph", "--phase",
                          dest="phase",
                          metavar="filename string",
                          type=str,
                          default="phases_00-22_by2.txt",
                          help='Should be a file with phases you wish to search for listed in a single column separated by newlines.\
                          Example files include "phases_00-22_by2.txt" or "phases_00-22_by4.txt" or "phases_00-20_by4.txt"')

    analysis.add_argument("-p","--period",
                          dest="period",
                          metavar="filename string",
                          type=str,
                          action='store',
                          default="period_24.txt",
                          help='Should be a file with phases you wish to search for listed in a single column separated by newlines.\
                          Provided file is "period_24.txt"')


    distribution = analysis.add_mutually_exclusive_group(required=False)
    distribution.add_argument("-e", "--exact",
                              dest="harding",
                              action='store_true',
                              default=False,
                              help="use Harding's exact null distribution (dflt)")
    distribution.add_argument("-g", "--gaussian","--normal",
                              dest="normal",
                              action='store_true',
                              default=False,
                              help="use normal approximation to null distribution")
    
    return p


# instantiate class to precalculate distribution
# usage: 
#   K = KendallTauP()
#   pval = K.pval(tau,n,two_tailed=True)
class KendallTauP:
    def __init__(self,N=150):        
        # largest number of samples to precompute
        self.N = N
        Nint = self.N*(self.N-1)/2

        # first allocate freq slots for largest sample array
        # as we fill this in we'll save the results for smaller samples

        # total possible number of inversions is Nint + 1
        freqN = np.zeros(Nint + 1)
        freqN[0] = 1.0

        # save results at each step in freqs array
        self.freqs = [np.array([1.0])]
        for i in xrange(1,self.N):
            last = np.copy(freqN)
            for j in xrange(Nint+1):
                # update each entry by summing over i entries to the left
                freqN[j] += sum(last[max(0,j-i):j])
            # copy current state into freqs array
            # the kth entry of freqs should have 1+k*(k-1)/2 entries
            self.freqs.append(np.copy(freqN[0:(1+(i+1)*i/2)]))
            
        # turn freqs into cdfs
        # distributions still with respect to number of inversions
        self.cdfs = []
        for i in xrange(self.N):
            self.cdfs.append(np.copy(self.freqs[i]))
            # turn into cumulative frequencies
            for j in xrange(1,len(self.freqs[i])):
                self.cdfs[i][j] += self.cdfs[i][j-1]
            # convert freqs to probs
            self.cdfs[i] = self.cdfs[i]/sum(self.freqs[i])
            
    # plot exact distribution compared to normal approx
    def plot(self,nlist):
        colors = cm.Set1(np.linspace(0,1,len(nlist)))

        # for plotting gaussian
        x = np.linspace(-1.2,1.2,300)
        # plot pdfs
        plt.figure()
        for i in xrange(len(nlist)):
            ntot = len(self.freqs[nlist[i]-1])-1
            tauvals = (ntot - 2.0*np.arange(len(self.freqs[nlist[i]-1])))/ntot
            probs = ((ntot+1.0)/2.0)*self.freqs[nlist[i]-1]/sum(self.freqs[nlist[i]-1])
            plt.scatter(tauvals,probs,color=colors[i])
            # now plot gaussian comparison
            var = 2.0*(2.0*nlist[i]+5.0)/(nlist[i]*(nlist[i]-1)*9.0)
            plt.plot(x,norm.pdf(x,0.0,np.sqrt(var)),color=colors[i])
        plt.legend(nlist,loc='best')
        # plt.savefig('pdfs.png')
        plt.show()

        # now plot cdfs
        plt.figure()
        for i in xrange(len(nlist)):
            ntot = len(self.freqs[nlist[i]-1])-1
            tauvals = -1.0*(ntot - 2.0*np.arange(len(self.freqs[nlist[i]-1])))/ntot
            probs = self.cdfs[nlist[i]-1]
            plt.scatter(tauvals,probs,color=colors[i])
            # now plot gaussian comparison
            var = 2.0*(2.0*nlist[i]+5.0)/(nlist[i]*(nlist[i]-1)*9.0)
            plt.plot(x,norm.cdf(x,0.0,np.sqrt(var)),color=colors[i])
        plt.legend(nlist,loc='best')
        # plt.savefig('cdfs.png')
        plt.show()

    # use cdfs to return pval
    # default to return two tailed pval
    def pval(self,tau,n,two_tailed=False):
        # enforce tau is between -1 and 1
        if tau <= -1.000001 or tau >= 1.000001:
            sys.stderr.write(str(type(tau))+"\n")
            sys.stderr.write(str(tau)+"\n")
            sys.stderr.write("invalid tau\n")
            #print 'invalid tau'
            return None
        # enforce n is less than our precomputed quantities
        if n > self.N:
            #print 'n is too large'
            sys.stderr.write("n is too large/n")
            return None

        # convert tau to value in terms of number of inversions
        ntot = n*(n-1)/2
        inv_score = int(round((ntot - tau * ntot)/2.0))
        # I'm a little worried about the precision of this,
        # but probably not enough to be really worried for reasonable n
        # since we really only need precision to resolve ntot points

        # if two tailed, we're getting a tail from a symmetric dist
        min_inv_score = min(inv_score,ntot-inv_score)

        if two_tailed:
            pval = self.cdfs[n-1][min_inv_score]*2.0
        else:
            # if one tailed return prob of getting that or fewer inversions
            pval = self.cdfs[n-1][inv_score]

        # if inv_score is 0, might have larger than 0.5 prob
        return min(pval,1.0)



if __name__=="__main__":
    parser = __create_parser__()
    args = parser.parse_args()
    main(args)
