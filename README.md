# PYIRD


Free from messy directory and file management of IRD analysis.
Currently under heavily construction. Use develop branch for dev.

Requirements
------------------------------------------

Recommended environment:

```
conda config --add channels http://ssb.stsci.edu/astroconda
conda create -n iraf37 python=3.7 iraf-all pyraf-all ds9
source activate iraf37
mkdir ~iraf
cd iraf
mkiraf
```

BASIC INSTALL
------------------------------------------
python setup.py install

Copy a Th-Ar list to the iraf linelist directory:

```
cp pyird/data/thar_ird2.dat /home/USERS/anaconda3/envs/iraf37/iraf/noao/lib/linelists/
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
