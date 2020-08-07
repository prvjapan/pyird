g++ -fPIC -Wall -O3 -ffast-math -msse2 -shared -o pyird/process_RNe.so src/process_RNe.cpp  -I /home/kawahara/anaconda3/pkgs/python-2.7.18-h15b4118_1/include/python2.7 -lboost_python-py27 -lboost_numpy-py27

cp process_RNe.so /home/kawahara/anaconda3/envs/py27/lib/python2.7/site-packages/pyird-0.0.0-py2.7.egg/pyird/

