#BooteJTK

BooteJTK is an implementation of empirical JTK (eJTK) on parametrically bootstrapped resamplings of time series.

Information on BooteJTK can be found in [Hutchison AL _et al._ (2016), "BooteJTK: Improved Rhythm Detection via Bootstrapping"](), available at bioRxiv.

Information on eJTK can be found in [Hutchison AL, Maienscein-Cline M, Chiang AH, _et al._ “Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data.” PLoS Computational Biology 2015 Mar. Vol. 11, No. 3, pp. e1004094. doi:10.1371/journal.pcbi.1004094](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004094)


##Instructions

###Run this code to compile the external module:

<pre><code>python setup.py build_ext --inplace</code></pre>


###Here is example code to test the method:

This code uses 10 bootstrap resamplings for data that has 2 replicates. It searches for phases at 0h, 2h, 4h, and so forth through 22h, and asymmetries at 2h, 4h, 6h, and so forth through 22h. It looks for a period of 24h, and gives the output files an identifier of _OTHERTEXT_.

<pre><code>./BooteJTK-CalcP.py -f example/TestInput4.txt -p ref_files/period24.txt -s ref_files/phases_00-22_by2.txt -a ref_files/asymmetries_02-22_by2.txt -x OTHERTEXT -r 2 -z 10</code></pre>

1. **example/TestInput4_OTHERTEXT_boot19-rep2_GammaP.txt**

   This is the output of BooteJTK with a Gamma fit of the null dataset used to assign p-values to the Tau values from running BooteJTK.

2. **example/TestInput4_OTHERTEXT_boot10-rep2.txt**
   This is the main output of BooteJTK, it contains the best reference waveform matching each time series. Best is defined as having the highest absolute Tau value. This will become input for CalcP.py, which is done automatically.


3. **example/TestInput4_OTHERTEXT_boot10-rep2_order_probs.pkl**
   This is further output of BooteJTK, a Python pickle file of two dictionaries. The first dictionary has keys which are IDs and contains a list of three lists, matching the number of unique time points in the time series. The first list contains the means, the second list contains the standard deviations, and the third list contains the number of replicates.
   The second dictionary has keys which are IDs and contains dictionaries of dictionaries whose keys are all the rank ordered time series sampled by BooteJTK and whose values are the frequency of those orderings.


4. **example/TestInput4_OTHERTEXT_boot10-rep2_order_probs_vars.pkl**

   This is further ouput of BooteJTK, a Python pickle file of two dictionaries.
   The first dictioanry has keys which are IDs and whose values are dictionaries whose keys are Tau values and whose values are probabilities of those Tau values.
   The second dictioanry has keys which are IDs and whose values are dictionaries whose keys are phase values and whose values are probabilities of those phase values.   

5. **example/TestInput4_NULL1000.txt**

   This is the file of randomly generated time series to be used to generate the null distribution of Tau values from which the Gamma null distribution will be fit.
   The other files that being with this prefix are the BooteJTK output analagous to the output from the original time series.


If you run the above command as is will produce files with a '_1' appended, as these files already exist in the examples folder.




##Version information:

* Python 2.7.11 (default, Dec  5 2015, 14:44:47)
* [GCC 4.2.1 Compatible Apple LLVM 7.0.0 (clang-700.1.76)] on darwin

* cython.__version__ 0.24
* scipy.__version__ 0.15.1
* numpy.__version__ 1.11.0
* statsmodels.__version__ 0.6.1

## License:
This code is released with the MIT License. See the License.txt file for more information.
