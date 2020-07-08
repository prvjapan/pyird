import pathlib

__all__ = ['FitsSet']


class FitsSet(object):
    
    def __init__(self, tag, obsmode = "MMF-MMF"):
        self.tag = tag
        self.datadir = None
        self.anadir = None

        #Ideintifiers
        self.targets = []
        self.flat_mmf1 = [] 
        self.flat_mmf2 = []
        self.comparison = []

    def idall(self):
        num=self.targets["id"]+self.flat_mmf1["id"]+self.flat_mmf2["id"]+self.comparison["id"]
        return num
    
    def rawfitslist(self,string=False):
        fitslist=extract_rawfitsset(self.comparison,string)
        a=extract_rawfitsset(self.targets,string)
        fitslist=fitslist+a
        a=extract_rawfitsset(self.flat_mmf1,string)
        fitslist=fitslist+a
        a=extract_rawfitsset(self.flat_mmf2,string)
        fitslist=fitslist+a
        return fitslist
    
    def rbfitslist(self,string=False):
        fitslist=extract_rbfitsset(self.comparison,self.anadir,string)
        a=extract_rbfitsset(self.targets,self.anadir,string)
        fitslist=fitslist+a
        a=extract_rbfitsset(self.flat_mmf1,self.anadir,string)
        fitslist=fitslist+a
        a=extract_rbfitsset(self.flat_mmf2,self.anadir,string)
        fitslist=fitslist+a
        return fitslist
            
def extract_rbfitsset(obj,odir,string=False):
    fitsarray=[]
    for num in obj["id"]:
        d=odir/(obj["tag"]+str(num)+"_rb.fits")
        if not d.exists():
            print(d, "does not exist.")
        else:
            if string:
                fitsarray.append(str(d))
            else:                
                fitsarray.append(d)
    return fitsarray
    
def extract_rawfitsset(obj):
    fitsarray=[]
    for num in obj["id"]:
        d=obj["dir"]/(obj["tag"]+str(num)+".fits")
        if not d.exists():
            print(d, "does not exist.")
        else:
            if string:
                fitsarray.append(str(d))
            else:                
                fitsarray.append(d)
    return fitsarray

def at_rblist(obj,adir,listname="tmp.list"):
    objl=extract_rbfitsset(obj,adir,string=True)
    atlistname=str(adir/listname)
    f=open(str(adir/listname),mode="w")
    for fits in objl:
        f.write(fits+"\n")
    f.close()
    atlistname='@'+atlistname
    return atlistname
