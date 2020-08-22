import pathlib
from pyraf import iraf
import os
import tqdm
__all__ = ['FitsSet']

class FitsSet(object):
    
    def __init__(self, tag, fitsdir, extension=""):
        self.tag = tag
        self.extension = extension
        self.fitsdir = fitsdir
        self.fitsid = None
        self.unlock=True #unlock for cleanning (i.e. remove)
        
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

    def clean(self):
        #Clean i.e. remove fits files if exists.
        import os
        if self.unlock:
            for i in range(len(self.path())):
                if self.path(check=False)[i].exists():
                    os.remove(self.path(check=False,string=True)[i])
                    print("rm old "+self.path(check=False,string=True)[i])
    
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
        combined_fitsset=FitsSet(tag,self.fitsdir)
        combined_fitsset.clean()
        iraf.imcombine(input=self.at_list(listname=tag),output=combined_fitsset.path(check=False)[0],combine=combine)
        return combined_fitsset
    
    def apall(self,apref,wavref=None,extout="_1d",extwout=""):
        currentdir=os.getcwd()
        os.chdir(str(self.fitsdir))
        apref_path=apref.path(string=False,check=True)[0].name #ref_ap
        iraf.imred()
        iraf.eche()

        for i,eachpath in enumerate(tqdm.tqdm(self.path(string=False,check=True))):
            iraf.apall(input=eachpath.name,output=(eachpath.name).replace(".fits",extout+".fits"),find="n",recenter="n",resize="n",edit="n",trace="n",fittrace="n",extract="y",references=apref_path,review="n",interactive="n")

        if wavref != None:
            wavref_path=wavref.path(string=False,check=True)[0].name #wav ref
            for i,eachpath in enumerate(tqdm.tqdm(self.path(string=False,check=True))):
                out=(eachpath.name).replace(".fits","w.fits")
                print("WAV",wavref_path, eachpath.name, out)
                iraf.refs(input=eachpath.name,references=wavref_path,select="match")
                iraf.dispcor(input=eachpath.name,output=out,flux="no")
        os.chdir(currentdir)
        return
    
    def apnormalize(self,tag,ref):
        currentdir=os.getcwd() 
        os.chdir(anadir)

        apnormalfits=FitsSet(tag,self.fitsdir)
        iraf.imred.echell.apnormaize(input=self.path(check=False,string=True)[0],output=apnormalfits.path(check=False,string=True)[0],find="n",recenter="n",resize="n",edit="y",trace="n",fittrace="n",references=ref,order=35,niterat=1)

        os.chdir(currentdir)
        return apnormalfits

    
