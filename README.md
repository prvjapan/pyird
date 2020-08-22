# PYIRD

Free from messy directory and file management of IRD analysis.

INSTALL
------------------------------------------
python setup.py install

g++ -fPIC -Wall -O3 -ffast-math -msse2 -shared -o process_RNe.so src/process_RNe.cpp  -I /home/kawahara/anaconda3/pkgs/python-3.6.18-h15b4118_1/include/python3.6 -lboost_python-py36 -lboost_numpy-py36

cp process_RNe.so /home/kawahara/anaconda3/envs/py36/lib/python3.6/site-packages/pyird-0.0.0-py3.6.egg/pyird/



Requirements
------------------------------------------

Recommended environment:

```
conda config --add channels http://ssb.stsci.edu/astroconda
conda create -n iraf36 python=3.6 iraf-all pyraf-all stsci
source activate iraf36
mkdir ~iraf
cd iraf
mkiraf
```




for the use of T.Hirano's code
-----------------------------------

Install boost (>1.63)

Compile such as (this is for python 36)

```
g++ -fPIC -Wall -O3 -ffast-math -msse2 -shared -o process_RNe.so process_RNe.cpp  -I /home/kawahara/anaconda3/pkgs/python-3.6.18-h15b4118_1/include/python3.6 -lboost_python-py36 -lboost_numpy-py36
```

for the use of hdsis_ecf
----------------------

```
cp pyraf/scripts/hdsis_ecf.cl /home/kawahara/.iraf/scripts/
```

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

For n=51 or 52 (YJ) and n=21 (H)
Mask 104 (YJ) 


Help tool for ECIDENTIFY
--------------------------

- calref.py

s option: The number of the orders varies depending on environment. s option gives an offset that defines the order numbering. Try it if you couldn't find good match. In particular, for SMF/YJ, try s = 1.

xorder 4
yorder 3