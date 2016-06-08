#BooteJTK

BooteJTK is an implementation of empirical JTK (eJTK) on parametrically bootstrapped resamplings of time series.

Information on BooteJTK can be found in [Hutchison AL _et al._ (2016), "BooteJTK: Improved Rhythm Detection via Bootstrapping"](), available at bioRxiv.

Information on eJTK can be found in [Hutchison AL, Maienscein-Cline M, Chiang AH, _et al._ “Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data.” PLoS Computational Biology 2015 Mar. Vol. 11, No. 3, pp. e1004094. doi:10.1371/journal.pcbi.1004094](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004094)


##Instructions

###Run this code to compile the external module:

<pre><code>python setup.py build_ext --inplace</code></pre>


###Here is example code to test the method:

<pre><code>./BooteJTK2-cython.py -f example/TestInput4.txt -p ref_files/period24.txt -s ref_files/phases_00-22_by2.txt -a ref_files/asymmetries_02-22_by2.txt -x cos24_ph00-22_by2_a02-22_by2 -r 2 -z 25</code></pre>







##Version information:

Python 2.7.11 (default, Dec  5 2015, 14:44:47)
[GCC 4.2.1 Compatible Apple LLVM 7.0.0 (clang-700.1.76)] on darwin

cython.__version__ 0.24
scipy.__version__ 0.15.1
numpy.__version__ 1.11.0
