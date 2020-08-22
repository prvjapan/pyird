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
    
    def extclean(self,extension):
        #Clean i.e. remove fits files if exists.
        import os
        if extension=="":
            print("extclean cannot clean w/o extension fits.")
            return
        
        for i in range(len(self.path())):
            if self.path(check=False)[i].exists():
                os.remove(self.extpath(extension,check=False,string=True)[i])
                print("rm old "+self.extpath(extension,check=False,string=True)[i])

    
    def remove_bias(self,rot=None,method = 'reference',hotpix_img = None):
        print("Bias Correction by M. KUZUHARA.")
        if rot=="r":
            print("180 degree rotation applied.")
        IRD_bias_sube.main(self.anadir,self.rawpath, method,rot,hotpix_im = hotpix_img)
        self.fitsdir=self.anadir
        self.extension="_rb"

    def rm_readnoise(self,maskfits,extin="_rb",extout="_rbo"):
        print("READ NOISE REDUCTION by H. Kawahara (OEGP).")
        currentdir=os.getcwd()
        os.chdir(str(self.anadir))

        ext_noexist, extf_noexist = self.check_existence(extin,extout)
        for i,fitsid in enumerate(tqdm.tqdm(ext_noexist)):
            filemask=maskfits.path()[0]
            processRN.wrap_kawahara_processRN(filen=ext_noexist[i],filemask=filemask,fitsout=extf_noexist[i])
        self.fitsdir=self.anadir
        self.extension=extout

        os.chdir(currentdir)
        
    def rm_readnoise_spline(self,maskfits):
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
            processRN.wrap_hirano_processRN(filen=rb_noexist[i],filemmf=maskfits.path()[0],fitsout=rbn_noexist[i])
        self.fitsdir=self.anadir
        self.extension="_rbn"

    def flatfielding1D(self,apflat,apref,wavref=None,extin="_rb",extout="_rb_f1d",extwout="_rb_f1dw",lower=-1,upper=2,badf="none"):
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
        
        ########################
        iraf.imred()
        iraf.eche()
        ## CHECK EXISTENCE for RB
        ext_noexist, extf_noexist = self.check_existence(extin,extout)
        
        for i,fitsid in enumerate(tqdm.tqdm(ext_noexist)):
            iraf.hdsis_ecf(inimg=ext_noexist[i],outimg=extf_noexist[i],plot=plotyn,st_x=lower,ed_x=upper,flatimg=apflat_path,ref_ap=apref_path,badfix=badf)

        if wavref != None:
            wavref_path=wavref.path(string=False,check=True)[0].name #wav ref
            ext_noexist, extonedw_noexist = self.check_existence(extin,extwout)
            for i,fitsid in enumerate(tqdm.tqdm(ext_noexist)):
                iraf.refs(input=extf_noexist[i],references=wavref_path,sort="mjd",group="mjd")
                iraf.dispcor(input=extf_noexist[i],output=extonedw_noexist[i])
                        
        os.chdir(currentdir)

    def apscat(self,apref,extin="_rb",extout="_rbs",ulimit=8,llimit=-8):
        ulimit=str(ulimit)
        llimit=str(llimit)
        apref_path=apref.path(string=False,check=True)[0].name #ref_ap

        currentdir=os.getcwd()
        os.chdir(str(self.anadir))
        iraf.imred()
        iraf.eche()
        ext_noexist, extf_noexist = self.check_existence(extin,extout)
        for i,fitsid in enumerate(tqdm.tqdm(ext_noexist)):
            iraf.apresize(input=ext_noexist[i],ulimit=ulimit,llimit=llimit)
            iraf.apscatter(input=ext_noexist[i],output=extf_noexist[i],references=apref_path,find="n",recenter="n",resize="n",edit="n",trace="n",fittrace="n",subtrac="y",smooth="y",fitscat="y",fitsmoo="y")
        
        os.chdir(currentdir)
        
    def extract1D(self,apref,wavref=None,extin="_rb",extout="_rb_1d",extwout="_rb_1dw"):
        currentdir=os.getcwd()
        os.chdir(str(self.anadir))
        apref_path=apref.path(string=False,check=True)[0].name #ref_ap
        
        ext_noexist, extoned_noexist = self.check_existence(extin,extout)
        iraf.imred()
        iraf.eche()

        for i,fitsid in enumerate(tqdm.tqdm(ext_noexist)):
            iraf.apall(input=ext_noexist[i],output=extoned_noexist[i],find="n",recenter="n",resize="n",edit="n",trace="n",fittrace="n",extract="y",references=apref_path,review="n",interactive="n")

        if wavref != None:
            wavref_path=wavref.path(string=False,check=True)[0].name #wav ref
            ext_noexist, extonedw_noexist = self.check_existence(extin,extwout)
            for i,fitsid in enumerate(tqdm.tqdm(ext_noexist)):
                iraf.refs(input=extoned_noexist[i],references=wavref_path,select="match")
                iraf.dispcor(input=extoned_noexist[i],output=extonedw_noexist[i],flux="no")
            
        os.chdir(currentdir)


    def check_existence(self,extin,extout):
        extf=self.extpath(extout,string=False,check=False)
        ext=self.extpath(extin,string=False,check=False)
        extf_noexist=[]
        ext_noexist=[]
        skip=0
        for i,extfi in enumerate(extf):
            if not extfi.exists():
                extf_noexist.append(str(extfi.name))
                ext_noexist.append(str(ext[i].name))
            else:
                skip=skip+1
            
        if skip>1:
            print("Skipped "+str(skip)+" files because they already exists.")

        return ext_noexist, extf_noexist
