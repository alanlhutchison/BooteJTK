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

#import cmath
from scipy.stats import circmean as sscircmean
from scipy.stats import circstd as sscircstd
#import scipy.stats as ss
import numpy as np
#from scipy.stats import kendalltau as kt
from scipy.stats import multivariate_normal as mn
from scipy.stats import rankdata
from scipy.stats import norm
from scipy.special import polygamma

import pickle
#from operator import itemgetter
import sys
import argparse
import time
import os.path

from get_stat_probs import get_stat_probs as gsp_get_stat_probs
from get_stat_probs import get_waveform_list as gsp_get_waveform_list
from get_stat_probs import make_references as gsp_make_references
from get_stat_probs import  kt ### this is kendalltau



def main(args):

    fn = args.filename
    fn_means = args.means
    fn_sds = args.sds
    fn_ns = args.ns
    
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
    reps = int(args.reps)

    ### If no list file set id_list to empty
    id_list = read_in_list(fn_list) if fn_list.split('/')[-1]!='DEFAULT' else []
    null_list = read_in_list(fn_null_list) if fn_null_list.split('/')[-1]!='DEFAULT' else []
        
    ### If no pkl out file, modify the option variable
    opt = 'premade' if os.path.isfile(fn_out_pkl) and fn_out_pkl.split('/')[-1]!='DEFAULT' else ''
    ### If not given a new name, name fn_out after the fn file
    if fn_out.split('/')[-1] == "DEFAULT":
        endstr = '_{0}_boot{1}-rep{2}.txt'.format(prefix,int(size),int(reps))
        fn_out = fn.replace('.txt',endstr) if  ".txt" in fn else fn+endstr

    ### If not given a new name, name fn_out_pkl based on the fn file
    if fn_out_pkl.split('/')[-1] == 'DEFAULT':
        endstr = '_{0}_boot{1}-rep{2}_order_probs.pkl'.format(prefix,int(size),int(reps))
        fn_out_pkl = fn.replace('.txt',endstr) if  ".txt" in fn else fn+endstr

    ### Name vars file based on pkl out file
    fn_out_pkl_vars = fn_out_pkl.replace('.pkl','_vars.pkl') 

    assert fn.split('/')[-1]!='DEFAULT' or fn_out_pkl.split('/')[-1]!='DEFAULT'   

    ### If we already have the PKL file, we just need a place to put the header information
    #if fn.split('/')[-1]=='DEFAULT' and fn_out_pkl.split('/')[-1]!='DEFAULT':
    #    d_series = dict(zip([key for key in d_data_master.keys()],[[]*len(d_data_master)]))
    #    fn_out= fn_out_pkl.replace('.pkl','_{0}-bootejtk.txt'.format(int(size)))
    #    new_header = [0,4,8,12,16,20,0,4,8,12,16,20] if fn_out_pkl[:-4]=='2.pkl' else [0,4,8,12,16,20]
    ### If we have the initial data we can get it 
    #elif fn.split('/')[-1]!='DEFAULT':
    ### WE HAVE CHANGED HOW THE DATA GETS PUT IN

    """ Read in the data """
    header,data = read_in(fn)
    d_series = dict_data(data)
    if fn_means!='DEFAULT' and fn_sds!='DEFAULT' and fn_ns!='DEFAULT':
        print 'Taking Limma/noreps input'
        header2,means = read_in(fn_means)
        _,sds = read_in(fn_sds)
        _,ns = read_in(fn_ns)
        d_data_master,new_header = get_data_multi(header,header2,means,sds,ns)
    else:
        d_data_master,new_header = get_data2(header,data)
        
        #new_header = list(new_header)*reps
    
        if 'premade' not in opt:
            print 'Running internal eBayes'
            D_null = get_series_data(d_data_master,null_list) if null_list!=[] else {}
            d_data_master = eBayes(d_data_master,D_null)
        elif 'premade' in opt:
            d_data_master,d_order_probs = pickle.load(open(fn_out_pkl,'rb'))

        
    def add_on_out(outfile):
        add_on = 1
        while os.path.isfile(outfile):
            print outfile, "already exists, take evasive action!!!"
            if '.txt' in outfile:
                end = '.txt'
            elif '.pkl' in outfile:
                end = '.pkl'
            origstr = end if add_on==1 else '_{0}'.format(add_on-1)+end
            outfile = outfile.replace(origstr,'_{0}'.format(add_on))+end
            add_on = add_on + 1
        return outfile

    
    fn_out = add_on_out(fn_out)
    fn_out_pkl = add_on_out(fn_out_pkl)
    fn_out_pkl_vars = add_on_out(fn_out_pkl_vars)    

    periods = np.array(read_in_list(fn_period),dtype=int)
    phases = np.array(read_in_list(fn_phase),dtype=int)
    widths = np.array(read_in_list(fn_width),dtype=int)

    waveform = 'cosine'

    triples = gsp_get_waveform_list(periods,phases,widths)

    dref = gsp_make_references(new_header,triples,waveform)
    
    d_data_master1 = {}
    d_order_probs_master = {}
    d_boots_master = {}

    Ps = []

    d_tau = {}
    d_ph = {}

    done = []
    remaining = []

    id_list = d_series.keys() if id_list==[] else id_list
    out_lines = []
    for geneID in d_data_master:

        ### If we have an ID list, we only want to deal with data from it.
        if geneID in id_list:

            if fn=='DEFAULT':
                mmax,mmin,MAX_AMP = np.nan,np.nan,np.nan
                sIQR_FC = np.nan
                smean = np.nan
                sstd = np.nan
                sFC = np.nan
            else:
                serie = d_series[geneID]                
                mmax,mmin,MAX_AMP=series_char(serie)
                sIQR_FC=IQR_FC(serie)
                smean = series_mean(serie)
                sstd = series_std(serie)
                sFC = FC(serie)

            s_stats = [smean,sstd,mmax,mmin,MAX_AMP,sFC,sIQR_FC,size]

            ### Need to make this file if it doesn't exist already
            if 'premade' not in opt:
                d_data_sub = {geneID:d_data_master[geneID]}
                d_data_master1 = dict(d_data_master1.items()+d_data_sub.items())
                
                d_order_probs,d_boots = get_order_prob(d_data_sub,size,reps)
                d_order_probs_master = dict(d_order_probs_master.items()+d_order_probs.items())
                d_boots_master = dict(d_boots_master.items()+d_boots.items())
                
            if geneID in d_order_probs:
                out1,out2,d_taugene,d_pergene,d_phgene,d_nagene = gsp_get_stat_probs(d_order_probs[geneID],new_header,triples,dref,size)
                out_line = [geneID,waveform]+out1+s_stats+out2
                out_line = [str(l) for l in out_line]
                out_lines.append(out_line)

                done.append(geneID)

                d_tau = dict(d_tau.items()+{geneID:d_taugene}.items())
                d_ph = dict(d_ph.items()+{geneID:d_phgene}.items())

        
        if len(remaining)>0:
            fn_remaining = fn_out.replace('.txt','_remaining_list.txt')
            with open(fn_remaining,'w') as g:
                for r in remaining:
                    g.write(r+'\n')
    pickle.dump([d_tau,d_ph],open(fn_out_pkl_vars,'wb'))                    
    pickle.dump([d_data_master1,d_order_probs_master,d_boots_master],open(fn_out_pkl,'wb'))
    #print out_lines
    taus = [[i,float(out[-2])] for i,out in enumerate(out_lines)]
    taus = sorted(taus,key=lambda x: np.abs(x[1]),reverse=True)
    #print taus
    indexes = np.array([i[0] for i in taus])
    #print indexes
    out_lines = np.array(out_lines)[np.array(indexes)]
    g = open(fn_out,'a')
    g.write("ID\tWaveform\tPeriodMean\tPeriodStdDev\tPhaseMean\tPhaseStdDev\tNadirMean\tNadirStdDev\tMean\tStd_Dev\tMax\tMin\tMax_Amp\tFC\tIQR_FC\tNumBoots\tTauMean\tTauStdDev\n")
    for out_line in out_lines:
        g.write("\t".join(out_line)+"\n")
    g.close()

    return fn_out,fn_out_pkl,header
    

    
def read_in_EMdata(fn):
    """Reads in one of the two EM files """
    WT = {}
    with open(fn,'r') as f:
        for line in f:
            words = line.strip('\n').split()
            if words[0]=='#':
                header = words
            else:
                key = words[0]
                m = [float(w) for w in words[1:7]]
                s = [float(s) for s in words[7:]]
                    
                WT[key] = [m,s,np.ones(6)*3]
    return WT


def get_SD_distr(dseries,reps,size):
    data = np.zeros(size)
    for j in xrange(size):
        ser = []
        for i in xrange(len(dseries[0])):
            ser.append(np.random.normal(dseries[0][i],dseries[1][i],size=reps))
        ser = np.concatenate(ser)
        SD = np.std(ser)
        data[j] = SD
    return data


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
        return np.isfinite(float(s))
    except TypeError and ValueError:
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
            if words[0] == "#" or words[0]=='ID':
                start_right = 1
                header = words[1:]
            else:
                if start_right == 0:
                    print "Please enter file with header starting with #"
                elif start_right == 1:
                    data.append(words)
    return header, data

def dict_data(data):
    d_series = {}
    for dat in data:
        d_series[dat[0]] = dat
    return d_series


def IQR_FC(series):
    qlo = __score_at_percentile__(series, 25)
    qhi = __score_at_percentile__(series, 75)
    if not is_number(qlo) or not is_number(qhi):
        return np.nan
    elif (qhi==0):
        return 0
    elif ( qlo==0):
        return np.nan
    else:
        #print qhi,qlo
        iqr = qhi/qlo
        return iqr

def FC(series):
    series=[float(s) if is_number(s) else 0 for s in series[1:]]
    if series!=[]:
        mmax = max(series)
        mmin = min(series)
        if mmin==0:
            sFC = -10000
        else:
            sFC = mmax / mmin
    else:
        sFC = np.nan
    return sFC


def series_char(series):
    """Uses interquartile range to estimate amplitude of a time series."""
    series=[float(s) for s in series[1:] if is_number(s)]
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
    series = [float(s) for s in series[1:] if is_number(s)]
    return np.mean(series)

def series_std(series):
    """Finds the std dev of a timeseries"""
    series = [float(s) for s in series[1:] if is_number(s)]
    return np.std(series)

def __score_at_percentile__(ser, per):
    ser = [float(se) for se in ser[1:] if is_number(se)]
    if len(ser)<5:
        score ="NA"
        return score
    else: 
        ser = np.sort(ser)
        i = int(per/100. * len(ser))
        if (i % 1 == 0):
            score = ser[i]
        else:
            interpolate = lambda a,b,frac: a + (b - a)*frac
            score = interpolate(ser[int(i)], ser[int(i) + 1], i % 1)
        return float(score)

def generate_mod_series(reference,series):
    """
    Takes the series from generate_base_null, takes the list from data, and makes a null
    for each gene in data or uses the one previously calculated.
    Then it runs Kendall's Tau on the exp. series against the null
    """
    tau,p=kt(series,reference)
    p = p / 2.0
    return tau,p

##################################
### HERE WE INSERT BOOT FUNCTIONS
##################################

def get_data(header,data):
    """This function does not use the header information to set the number of replicates """
    new_h = [float(h[2:])%24 if 'ZT' in h or 'CT' in h else float(h)%24 for h in header]
    #print new_h
    length = len(new_h)
    seen = []
    dref = {}
    for i,h in enumerate(new_h):
        if h not in seen:
            seen.append(h)
            dref[h]=[i]
        else:
            dref[h].append(i)
    d_data = {}
    #print dref
    for dat in data:
        name=dat[0]
        series = [float(da) if is_number(da) else np.nan for da in dat[1:]]
        if len(series)==length:
            out = [[],[],[]]
            for i,s in enumerate(seen):
                points = [series[idx] for idx in dref[s]]
                N = len([p for p in points if not np.isnan(p)])
                m = np.nanmean(points)
                std = np.nanstd(points)
                out[0].append(m)
                out[1].append(std)
                out[2].append(N)
            #print name,seen,out
            d_data[name]=out
    #print d_data.keys()
    return d_data,seen


def get_data2(header,data):
    """ This function uses the header information to set the number of replicates"""
    new_h = [float(h[2:])%24 if 'ZT' in h or 'CT' in h else float(h)%24 for h in header]
    length = len(new_h)
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
        series = [float(da) if is_number(da) else np.nan for da in dat[1:]]

        out = [[],[],[]]
        for i,s in enumerate(new_h):
            points = [series[idx] for idx in dref[s]]
            N = len([p for p in points if not np.isnan(p)])
            m = np.nanmean(points)
            std = np.nanstd(points)
            out[0].append(m)
            out[1].append(std)
            out[2].append(N)
            #print name,seen,out
        d_data[name]=out
    #print d_data.keys()
    return d_data,new_h



def get_data_multi(header,header2,data,sds,ns):
    """ This function takes in several pre-eBayes files to create the d_data dictionary"""
    new_h = [float(h[2:])%24 if 'ZT' in h or 'CT' in h else float(h)%24 for h in header]
    h2 = [float(h[2:])%24 if 'ZT' in h or 'CT' in h else float(h)%24 for h in header2]    
    length = len(new_h)
    seen = []
    dref = {}
    for i,h in enumerate(new_h):
        seen.append(h2.index(h))
    d_data = {}
    for j,dat in enumerate(data):
        g_sds = sds[j]
        g_ns = ns[j]
        name=dat[0]
        g_means = [float(da) if is_number(da) else np.nan for da in dat[1:]]
        g_means = [g_means[i] for i in seen]
        g_sds =   [float(da) if is_number(da) else np.nan for da in g_sds[1:]]
        g_sds = [g_sds[i] for i in seen]
        g_ns =   [float(da) if is_number(da) else np.nan for da in g_ns[1:]]
        g_ns = [g_ns[i] for i in seen]
        out = [g_means,g_sds,g_ns]
        d_data[name]=out

    return d_data,new_h




def get_series_data(d_data_master,id_list):
    d_data = {}
    for key in d_data_master:
        if key in id_list:
            dataset = d_data_master[key]
            #print key,datase
            N = np.sum(dataset[2])
            length = len(dataset[0])
            one = 1./N*np.sum([dataset[2][i]*(dataset[1][i]**2+dataset[0][i]**2) for i in xrange(length)])
            two = (1./N * np.sum([dataset[2][i]*dataset[0][i] for i in xrange(length)]))**2
            std = np.sqrt(one+two)
            m = 1./N *np.sum([dataset[2][i]*dataset[0][i] for i in xrange(length)])
            d_data[key] = [m,std,N]
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
        #print D_null.keys()    
        #print dg,s
        G = len(d_data)
        s2 = np.array([ss for ss in s if ss!=0])
        dg2 = np.array([dg[i] for i,ss in enumerate(s) if ss!=0])
        G = len(dg2)
        z = 2.*np.log(s2)

        e = z - digamma(dg2/2) + np.log(dg2/2)
        emean = np.nanmean(e)
        d0 = 2.* solve_trigamma( np.nanmean( (e-emean)**2 *G/(G-1)-trigamma(dg2/2) )   )
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

def get_order_prob(d_data,size,reps):
    ### WE WANT TO CHANGE HOW REPS GET INCORPORATED...
    d_order_prob = {}
    d_boots = {}
    for key in d_data:
        res = d_data[key]
        d_order_prob[key],d_boots[key]=dict_of_orders(res[0],res[1],res[2],size)        
        #d_order_prob[key],d_boots[key]=dict_of_orders(list(res[0])*reps,list(res[1])*reps,list(res[2])*reps,size)
    return d_order_prob,d_boots

def dict_of_orders(M,SDS,NS,size):
    """
    This produces a dictionary of probabilities
    for the different orders given the data
    """
    index = range(len(M))
    dorder,s2 = dict_order_probs(M,SDS,NS,size)
    d = {}
    for idx in dorder:
        d[idx]=dorder[idx]
    SUM = sum(d.values())
    for key in d:
        d[key]=d[key]/SUM
    return d,s2

def dict_order_probs(ms,sds,ns,size=100):
    """ Check for SD problems """
    sds = [sd if is_number(sd) else np.nanmean(sds) for sd in sds]
    d = {}
    """ Create array for bootstraps of time series """
    s3 = np.zeros((size,len(sds)))
    """ Fill array time point by time point """
    for i in xrange(len(sds)):
        s3[:,i] = np.random.normal(ms[i],sds[i],size=size)
    """ Turn into ranks for storage """
    for s in s3:
        r=tuple(map(int,rankdata(s)))
        d.setdefault(r,0)
        d[r] += 1./size
    return d,s3


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
                    default='DEFAULT',                          
                   type=str,
                   help='This is the filename of the data series you wish to analyze.\
                   The data should be tab-spaced. The first row should contain a # sign followed by the time points with either CT or ZT preceding the time point (such as ZT0 or ZT4). Longer or shorter prefixes will not work. The following rows should contain the gene/series ID followed by the values for every time point. Where values are not available NA should be put in it\'s place.')

    analysis.add_argument("-F", "--means",
                   dest="means",
                   action='store',
                   metavar="filename string",
                    default='DEFAULT',                          
                   type=str,
                   help='This is the filename of the time point means of the data series you wish to analyze. The data should be tab-spaced. The first row should contain a # sign followed by the time points with either CT or ZT preceding the time point (such as ZT0 or ZT4). Longer or shorter prefixes will not work. The following rows should contain the gene/series ID followed by the values for every time point. Where values are not available NA should be put in it\'s place.')

    analysis.add_argument("-S", "--sds",
                   dest="sds",
                   action='store',
                   metavar="filename string",
                    default='DEFAULT',                          
                   type=str,
                   help='This is the filename of the time point standard devations of the data series you wish to analyze. The data should be tab-spaced. The first row should contain a # sign followed by the time points with either CT or ZT preceding the time point (such as ZT0 or ZT4). Longer or shorter prefixes will not work. The following rows should contain the gene/series ID followed by the values for every time point. Where values are not available NA should be put in it\'s place.')

    analysis.add_argument("-N", "--ns",
                   dest="ns",
                   action='store',
                   metavar="filename string",
                    default='DEFAULT',
                   type=str,
                   help='This is the filename of the time point replicate numbers of the data series you wish to analyze. \
                   The data should be tab-spaced. The first row should contain a # sign followed by the time points with either CT or ZT preceding the time point (such as ZT0 or ZT4). Longer or shorter prefixes will not work. The following rows should contain the gene/series ID followed by the values for every time point. Where values are not available NA should be put in it\'s place.')
    
    

    analysis.add_argument('-x',"--prefix",
                          dest="prefix",
                          type=str,
                          metavar="string",
                          action='store',
                          default="",
                          help="string to be inserted in the output filename for this run")


    analysis.add_argument('-r',"--reps",
                          dest="reps",
                          type=int,
                          metavar="int",
                          action='store',
                          default=2,
                          help="# of reps of each time point to bootstrap (1 or 2, generally)")

    analysis.add_argument('-z',"--size",
                          dest="size",
                          type=int,
                          metavar="int",
                          action='store',
                          default=50,
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




if __name__=="__main__":
    parser = __create_parser__()
    args = parser.parse_args()
    main(args)
