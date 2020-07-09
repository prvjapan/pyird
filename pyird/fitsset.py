import pathlib
import IRD_bias_sube
import processRN
from pyraf import iraf
import tqdm

__all__ = ['FitsSet']

class FitsSet(object):
    
    def __init__(self, tag, fitsdir, extension=""):
        self.tag = tag
        self.extension = extension
        self.fitsdir = fitsdir
        self.fitsid = None
        
    def path(self,string=False,check=True):
        fitsarray=[]
        if self.fitsid is None:
            d=self.fitsdir/(self.tag+self.extension+".fits")
            if string:
                fitsarray.append(str(d))
            else:                
                fitsarray.append(d)
        else:
            for num in self.fitsid:
                d=self.fitsdir/(self.tag+str(num)+self.extension+".fits")
                if check and not d.exists():
                    print(d, "does not exist.")
                else:
                    if string:
                        fitsarray.append(str(d))
                    else:                
                        fitsarray.append(d)
        return fitsarray


    
class FitsStream(FitsSet):    
    def __init__(self, streamid, rawdir, anadir, rawtag="IRDA0000", extension=""):
        super(FitsStream,self).__init__(rawtag, rawdir, extension="")
        self.streamid = streamid
        self.rawdir = rawdir
        self.anadir = anadir
        self.rbcombined = FitsSet(streamid+"_rbcombined",self.anadir)
        self.rbncombined = FitsSet(streamid+"_rbncombined",self.anadir)
        
    def rawpath(self,string=False,check=True):
        return self.path(string,check)

    def rbpath(self,string=False,check=True):
        #rb files after BIAS CORRECTION
        self.fitsdir=self.anadir
        self.extension="_rb"
        return self.path(string,check)
    
    def rbnpath(self,string=False,check=True):
        #rbn files after BIAS/READNOISE correction
        self.fitsdir=self.anadir
        self.extension="_rbn"
        return self.path(string,check)        
 
    def remove_bias(self,method = 'reference',hotpix_img = None):
        print("Bias Correction by M. KUZUHARA.")
        IRD_bias_sube.main(self.anadir,self.rawpath(), method, hotpix_im = hotpix_img)
                
    def at_list(self,objl):
        listname=self.streamid+".list"
        atlistname=str(self.anadir/listname)
        f=open(str(self.anadir/listname),mode="w")
        for fits in objl:
            f.write(fits+"\n")
        f.close()
        atlistname='@'+atlistname
        return atlistname

    def combine_rb(self,combine="median"):
        iraf.imcombine(input=self.at_list(self.rbpath(string=True)),output=self.rbcombined.path(check=False)[0],combine=combine)
        
    def combine_rbn(self,combine="median"):
        iraf.imcombine(input=self.at_list(self.rbnpath(string=True)),output=self.rbncombined.path(check=False)[0],combine=combine)
        
    def rm_readnoise(self,maskfits):
        print("READ NOISE REDUCTION by T. Hirano.")
        rbn=self.rbnpath(string=True,check=False)
        rb=self.rbpath(string=True,check=False)
        for i,fitsid in enumerate(tqdm.tqdm(self.fitsid)):
            processRN.wrap_processRN(filen=rb[i],filemmf=maskfits.path()[0],fitsout=rbn[i])

    
        
