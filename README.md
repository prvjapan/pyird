# PYIRD

Free from messy directory and file management of IRD analysis.

INSTALL
------------------------------------------
python setup.py install

g++ -fPIC -Wall -O3 -ffast-math -msse2 -shared -o process_RNe.so src/process_RNe.cpp  -I /home/kawahara/anaconda3/pkgs/python-2.7.18-h15b4118_1/include/python2.7 -lboost_python-py27 -lboost_numpy-py27

cp process_RNe.so /home/kawahara/anaconda3/envs/py27/lib/python2.7/site-packages/pyird-0.0.0-py2.7.egg/pyird/



Requirements
------------------------------------------

Recommended environment:

```
conda config --add channels http://ssb.stsci.edu/astroconda
conda create -n iraf27 python=2.7 iraf-all pyraf-all stsci
source activate iraf27
mkdir ~iraf
cd iraf
mkiraf
```

Install boost (>1.63)

Compile such as (this is for python 27)

```
g++ -fPIC -Wall -O3 -ffast-math -msse2 -shared -o process_RNe.so process_RNe.cpp  -I /home/kawahara/anaconda3/pkgs/python-2.7.18-h15b4118_1/include/python2.7 -lboost_python-py27 -lboost_numpy-py27
```

python setup.py install

Classes
------------------

- fitsset.FitsSet --  sets of fits files 
- irdstream.Stream2D -- 2D fits stream from raw images


Scripts
------------------------------

- makemask.py (by M. Kuzuhara, modified by H. Kawahara)
- IRD_bias_sube.py (by M. Kuzuhara, modified by H. Kawahara)
- process_RNe.cpp (by T. Hirano, wrapped by H. Kawahara)


Aperture
------------------------------

For n=51 (YJ) and n=21 (H)
Mask 104 (YJ) 
