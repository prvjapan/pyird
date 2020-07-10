import pathlib
from pyraf import iraf

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

    def at_list(self,listname="tmp"):
        listname=listname+".list"
        atlistname=str(self.fitsdir/listname)
        f=open(atlistname,mode="w")
        for fits in self.path(string=True):
            f.write(fits+"\n")
        f.close()
        atlistname='@'+atlistname
        return atlistname

    def imcombine(self,tag,combine="median"):
        combined_fitsset=FitsSet(tag,self.anadir)
        iraf.imcombine(input=self.at_list(listname=tag),output=combined_fitsset.path(check=False)[0],combine=combine)
        return combined_fitsset
    

    
        
