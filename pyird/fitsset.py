import pathlib
import IRD_bias_sube
from pyraf import iraf

__all__ = ['FitsSet']


class FitsSet(object):
    
    def __init__(self, tag, fitsdir, extension=""):
        self.tag = tag
        self.extension = extension
        self.fitsdir = fitsdir
        self.fitsid = []
        
    def path(self,string=False):
        fitsarray=[]
        for num in self.fitsid:
            d=self.fitsdir/(self.tag+str(num)+self.extension+".fits")
            if not d.exists():
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
        self.combined = self.anadir/(self.streamid+"_rbcombined.fits")
        
    def rawpath(self,string=False):
        return self.path(string)

    def rbpath(self,string=False):
        self.fitsdir=self.anadir
        self.extension="_rb"
        return self.path(string)        
        
    def remove_bias(self,method = 'reference',hotpix_img = None):
        print("Bias Correction")
        IRD_bias_sube.main(self.anadir,self.rawpath(), method, hotpix_im = hotpix_img)
        
        
    def at_rblist(self):
        listname=self.streamid+".list"
        objl=self.rbpath(string=True)
        atlistname=str(self.anadir/listname)
        f=open(str(self.anadir/listname),mode="w")
        for fits in objl:
            f.write(fits+"\n")
        f.close()
        atlistname='@'+atlistname
        return atlistname

    def rb_combine(self,combine="median"):
        iraf.imcombine(input=self.at_rblist(),output=self.combined,combine=combine)


        
