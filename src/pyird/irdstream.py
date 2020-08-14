import pathlib
from pyird.fitsset import FitsSet
from pyraf import iraf
from pyird import IRD_bias_sube
from pyird import processRN
import tqdm
import os 
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

    def flatfielding(self,apflat,apref,extin="rb",extout="rbf",lower=-1,upper=2,badf="none"):
        iraf.task(hdsis_ecf = "home$scripts/hdsis_ecf.cl")
        currentdir=os.getcwd()
        os.chdir(str(self.anadir))

        ####CHECK THESE VALUES##
        plotyn="no"    #plot
        apflat_path=apflat.path(string=False,check=True)[0].name #ref_ap
        apref_path=apref.path(string=False,check=True)[0].name #ref_ap
        #IS3
        lower=str(lower)#"-1"    #st_x
        upper=str(upper)#"2"      #ed_x
        #badf="none"  #badfix
          #badfix
        
        ########################
        iraf.imred()
        iraf.eche()
        ## CHECK EXISTENCE for RB
        ext_noexist, extf_noexist = self.check_existence(extin,extout)
        
        for i,fitsid in enumerate(tqdm.tqdm(ext_noexist)):
            iraf.hdsis_ecf(inimg=ext_noexist[i],outimg=extf_noexist[i],plot=plotyn,st_x=lower,ed_x=upper,flatimg=apflat_path,ref_ap=apref_path,badfix=badf)

        os.chdir(currentdir)

    def extract1D(self,apref,extin="rb",extout="rboned"):
        currentdir=os.getcwd()
        os.chdir(str(self.anadir))
        apref_path=apref.path(string=False,check=True)[0].name #ref_ap

        ## CHECK EXISTENCE
        ext_noexist, extoned_noexist = self.check_existence(extin,extout)
        for i,fitsid in enumerate(tqdm.tqdm(ext_noexist)):
            iraf.imred.echell.apall(input=ext_noexist[i],output=extoned_noexist[i],find="n",recenter="n",resize="n",edit="n",trace="n",fittrace="n",extract="y",references=apref_path,review="n",interactive="n")

        os.chdir(currentdir)


    def check_existence(self,extin,extout):
        extf=self.extpath("_"+extout,string=False,check=False)
        ext=self.extpath("_"+extin,string=True,check=False)
        extf_noexist=[]
        ext_noexist=[]
        skip=0
        for i,extfi in enumerate(extf):
            if not extfi.exists():
                extf_noexist.append(str(extfi))
                ext_noexist.append(str(ext[i]))
            else:
                skip=skip+1
            
        if skip>1:
            print("Skipped "+str(skip)+" files because they already exists.")

        return ext_noexist, extf_noexist
