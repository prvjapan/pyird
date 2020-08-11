import pathlib
from pyird.fitsset import FitsSet
from pyraf import iraf
from pyird import IRD_bias_sube
from pyird import processRN
import tqdm

__all__ = ['Stream1D','Stream2D']

class Stream1D(FitsSet):    
    def __init__(self, streamid, rawdir, anadir, rawtag="IRDA000", extension=""):
        super(Stream1D,self).__init__(rawtag, rawdir, extension="")
        self.streamid = streamid
        self.anadir = anadir
        self.unlock=False
        
class Stream2D(FitsSet):    
    def __init__(self, streamid, rawdir, anadir, rawtag="IRDA000", extension=""):
        super(Stream2D,self).__init__(rawtag, rawdir, extension="")
        self.streamid = streamid
        self.rawdir = rawdir
        self.anadir = anadir
        self.unlock=False
        
    @property
    def fitsid(self):
        return self._fitsid

    @fitsid.setter
    def fitsid(self, fitsid):
        self._fitsid = fitsid
        self.rawpath=self.path(string=False,check=True)

    def fitsid_increment(self):
        for i in range(0,len(self.fitsid)):
            self.fitsid[i]=self.fitsid[i]+1
        self.rawpath=self.path(string=False,check=True)
            
    def fitsid_decrement(self):
        for i in range(0,len(self.fitsid)):
            self.fitsid[i]=self.fitsid[i]-1
        self.rawpath=self.path(string=False,check=True)
            
    def extpath(self,extension,string=False,check=True):
        f=self.fitsdir
        e=self.extension
        self.fitsdir=self.anadir
        self.extension=extension
        path_=self.path(string,check)
        self.fitsdir=f
        self.extension=e
        return path_
        
    def remove_bias(self,rot=None,method = 'reference',hotpix_img = None):
        print("Bias Correction by M. KUZUHARA.")
        if rot=="r":
            print("180 degree rotation applied.")
        IRD_bias_sube.main(self.anadir,self.rawpath, method,rot,hotpix_im = hotpix_img)
        self.fitsdir=self.anadir
        self.extension="_rb"
    
    def rm_readnoise(self,maskfits):
        print("READ NOISE REDUCTION by T. Hirano.")
        rbn=self.extpath("_rbn",string=False,check=False)
        rb=self.extpath("_rb",string=True,check=False)

        rbn_noexist=[]
        rb_noexist=[]
        skip=0
        for i,rbni in enumerate(rbn):
            if not rbni.exists():
                rbn_noexist.append(str(rbni))
                rb_noexist.append(str(rb[i]))
            else:
                skip=skip+1
            
        if skip>1:
            print("Read Noise Correction: Skipped "+str(skip)+" files because they already exists.")

        for i,fitsid in enumerate(tqdm.tqdm(rb_noexist)):
            processRN.wrap_processRN(filen=rb_noexist[i],filemmf=maskfits.path()[0],fitsout=rbn_noexist[i])
        self.fitsdir=self.anadir
        self.extension="_rbn"
