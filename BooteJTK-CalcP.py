#!/usr/bin/env python
"""
Created on Nov 1, 2015
@author: Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron R. Dinner Group, University of Chicago

This script is a bootstrapped expansion of the eJTK method described in

Hutchison AL, Maienschein-Cline M, and Chiang AH et al. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data, PLoS Computational Biology 2015 11(3): e 1004094. doi:10.1371/journal.pcbi.1004094

This script bootstraps time series and provides phase and tau distributions from those bootstraps to allow for measurement of the variance on phase and tau estimates.


Please use ./BooteJTK -h to see the help screen for further instructions on running this script.

"""
VERSION="1.0"

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
import pandas as pd
import pickle
#from operator import itemgetter
import sys
import argparse
import time
import os.path

import subprocess
#from get_stat_probs import get_stat_probs as gsp_get_stat_probs
#from get_stat_probs import get_waveform_list as gsp_get_waveform_list
#from get_stat_probs import make_references as gsp_make_references
#from get_stat_probs import  kt ### this is kendalltau

binpath=os.path.join(os.path.dirname(sys.argv[0]),'bin/')
sys.path.insert(1,binpath)

import BooteJTK
import CalcP

def main(args):

    fn = args.filename

    ### INTEGRATING THIS INTO THE R CODE
    #fn_means = args.means
    #fn_sds = args.sds
    #fn_ns = args.ns
    
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


    if args.basic==False:
        args.limma==True
        args.vash ==True
    elif args.basic==True and args.limma==True and args.vash==False:
        args.limma=True
        args.vash =False

        
    #fn_null = args.null
    if args.noreps==True:
        args.prefix = 'NoRepSD_'+args.prefix
        print 'No replicates, skipping Limma procedure'
        print 'Estimating time point variance from arrhythmic genes'
        try:
            df = pd.read_table(fn,index_col='ID')
        except ValueError:
            df = pd.read_table(fn,index_col='#')
        except ValueError:
            print 'Header needs to begin with "ID" or with "#"'

        j = pd.read_table(args.jtk,index_col='ID')
        mean = df.loc[j[j.GammaP>0.8].index].std(axis=1).dropna().mean()

        df_sds = pd.DataFrame(np.ones(df.shape)*mean,index=df.index,columns=df.columns)
        df_ns =  pd.DataFrame(np.ones(df.shape),index=df.index,columns=df.columns)
        fn_sds = fn.replace('.txt','_Sds_noRepsEst.txt')
        df_sds.to_csv(fn_sds,na_rep=np.nan,sep='\t')
        fn_ns = fn.replace('.txt','_Ns_noRepsEst.txt')
        df_sds.to_csv(fn_ns,na_rep=np.nan,sep='\t')
                      
        args.means = fn
        args.sds = fn_sds
        args.ns = fn_ns
        
    elif args.limma==True:

        """Rscript command for Limma"""
        if args.vash==False:
            print 'Running the Limma commands'
            args.prefix = 'Limma_'+args.prefix
            path2script = binpath+'Limma_voom_script.R'
            args.means = fn.replace('.txt','_Means_postLimma.txt')
            args.sds = fn.replace('.txt','_Sds_postLimma.txt')
            args.ns = fn.replace('.txt','_Ns_postLimma.txt')
        else:
            print 'Running the Vash commands'
            args.prefix = 'Vash_'+args.prefix
            path2script = binpath+'Limma_voom_vash_script.R'
            args.means = fn.replace('.txt','_Means_postVash.txt')
            args.sds = fn.replace('.txt','_Sds_postVash.txt')
            args.ns = fn.replace('.txt','_Ns_postVash.txt')

        command = 'Rscript'

        pref=fn.replace('.txt','')
        period='24'
        if args.rnaseq:
            arguments = [fn, pref, period,'rnaseq']
        else:
            arguments = [fn, pref, period]            
        cmd = [command, path2script] + arguments
        ret = subprocess.call(cmd)
        #print 'Subprocess call is ', ret
        if ret != 0:
            if ret < 0:
                print "Killed by signal", -ret
            else:
                print "Command failed with return code", ret
            assert False
    else:
        pass
    fn_out,fn_out_pkl,header = BooteJTK.main(args)

    #args.output = fn_out.replace('boot','NULL1000-boot')
    #args.pickle = fn_out_pkl.replace('boot','NULL1000-boot')

    #print args.pickle
    #print args.output
    fn_null_out = args.output

    fn_null = fn.replace('.txt','_NULL1000.txt')
    
    sims = 1000
    with open(fn_null,'w') as g:
        g.write('\t'.join(['#']+header)+'\n')
        for i in xrange(sims):
            line = ['wnoise_'+str(i)] + map(str,list(np.random.normal(0,1,len(header))))
            g.write('\t'.join(line)+'\n')
            
    args.filename = fn_null
    if args.noreps==True:
        print 'No replicates, skipping Limma procedure'
        print 'Estimating time point variance from arrhythmic genes'
        try:
            df = pd.read_table(fn,index_col='ID')
        except ValueError:
            df = pd.read_table(fn,index_col='#')
        except ValueError:
            print 'Header needs to begin with "ID" or with "#"'

        j = pd.read_table(args.jtk,index_col='ID')
        mean = df.loc[j[j.GammaP>0.8].index].std(axis=1).dropna().mean()

        df_sds = pd.DataFrame(np.ones(df.shape)*mean,index=df.index,columns=df.columns)
        df_ns =  pd.DataFrame(np.ones(df.shape),index=df.index,columns=df.columns)
        fn_sds = fn_null.replace('.txt','_Sds_noRepsEst.txt')
        df_sds.to_csv(fn_sds,na_rep=np.nan,sep='\t')
        fn_ns = fn_null.replace('.txt','_Ns_noRepsEst.txt')
        df_sds.to_csv(fn_ns,na_rep=np.nan,sep='\t')
                      
        args.means = fn_null
        args.sds = fn_sds
        args.ns = fn_ns
    elif args.limma==True:
        """Rscript command for Limma"""
        if args.vash==False:
            path2script = binpath+'Limma_voom_script.R'
            args.means = fn_null.replace('.txt','_Means_postLimma.txt')
            args.sds = fn_null.replace('.txt','_Sds_postLimma.txt')
            args.ns = fn_null.replace('.txt','_Ns_postLimma.txt')
        else:
            path2script = binpath+'Limma_voom_vash_script.R'
            args.means = fn_null.replace('.txt','_Means_postVash.txt')
            args.sds = fn_null.replace('.txt','_Sds_postVash.txt')
            args.ns = fn_null.replace('.txt','_Ns_postVash.txt')

        command = 'Rscript'
        pref=fn_null.replace('.txt','')
        period='24'
        arguments = [fn_null, pref, period]
        cmd = [command, path2script] + arguments
        subprocess.call(cmd)    
    else:
        pass
    
    fn_null_out,_,_ = BooteJTK.main(args)
    args.filename = fn_out
    args.null = fn_null_out
    args.fit = ''    
    CalcP.main(args)


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
    analysis.add_argument("-F", "--means",
                   dest="means",
                   action='store',
                   metavar="filename string",
                    default='DEFAULT',                          
                   type=str,
                   help='This is the filename of the time point means of the data series you wish to analyze.\
                   The data should be tab-spaced. The first row should contain a # sign followed by the time points with either CT or ZT preceding the time point (such as ZT0 or ZT4). Longer or shorter prefixes will not work. The following rows should contain the gene/series ID followed by the values for every time point. Where values are not available NA should be put in it\'s place.')

    analysis.add_argument("-S", "--sds",
                   dest="sds",
                   action='store',
                   metavar="filename string",
                    default='DEFAULT',                          
                   type=str,
                   help='This is the filename of the time point standard devations of the data series you wish to analyze.\
                   The data should be tab-spaced. The first row should contain a # sign followed by the time points with either CT or ZT preceding the time point (such as ZT0 or ZT4). Longer or shorter prefixes will not work. The following rows should contain the gene/series ID followed by the values for every time point. Where values are not available NA should be put in it\'s place.')

    analysis.add_argument("-N", "--ns",
                   dest="ns",
                   action='store',
                   metavar="filename string",
                    default='DEFAULT',
                   type=str,
                   help='This is the filename of the time point replicate numbers of the data series you wish to analyze.                         \
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


    analysis.add_argument("-V","--vash",
                          dest="vash",
                          action='store_true',
                          default=False,
                          help='Determine if you would like to use limma or Vash')

    analysis.add_argument("-U","--noreps","--unique",
                          dest="noreps",
                          action='store_true',
                          default=False,
                          help='Determine if your data has no replicates and therefore the standard deviation should be estimated from the arrhythmic time series')

    analysis.add_argument("-R","--rnaseq",
                          dest="rnaseq",
                          action='store_true',
                          default=False,
                          help='Flag for data that is RNA-Seq and for which voom should be used.')
    
    analysis.add_argument("-L","--limma",
                          dest="limma",
                          action='store_true',
                          default=False,
                          help='Flag for using the limma variance shrinkage methods')
    

    analysis.add_argument("-J","--jtk",
                          dest="jtk",
                          metavar="filename string",
                          type=str,
                          action='store',
                          help="The eJTK file to use if you don't have replicates in in your time series. The standard deviation between points will be estimated based on the arrhythmic time series.")

    analysis.add_argument("-B","--basic",
                          dest="basic",
                          action='store_true',
                          default=False,
                          help='Flag for not using the limma or vash settings')

    
    

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
