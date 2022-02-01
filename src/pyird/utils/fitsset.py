"""Set of fits files."""

#from pyraf import iraf
import sys
__all__ = ['FitsSet']


class FitsSet(object):

    def __init__(self, tag, fitsdir, extension=''):
        self.tag = tag
        self.extension = extension
        self.fitsdir = fitsdir
        self.fitsid = None
        self.unlock = True  # unlock for cleanning (i.e. remove)
        self.not_ignore_warning = True

    def data(self, indices=None):
        """get header and data in an array form
        Args:
            indices: indices of fitsset, if None, 0 is used

        Returns:
            array of data
        """
        import astropy.io.fits as pyf
        import numpy as np
        if indices == None:
            hdulist = pyf.open(np.array(self.path(string=True))[0])
            hdu = hdulist[0]
            return np.array(hdu.data)

        data_arr = []
        for filen in np.array(self.path(string=True))[np.array(indices)]:
            hdulist = pyf.open(filen)
            hdu = hdulist[0]
            data_arr.append(np.array(hdu.data))
        return data_arr

    def path(self, string=False, check=True):
        """get path array.

        Returns:
            array of paths
        """
        fitsarray = []
        if self.fitsid is None:
            print(self.fitsdir, self.tag, self.extension)
            d = self.fitsdir/(self.tag+self.extension+'.fits')
            if string:
                fitsarray.append(str(d))
            else:
                fitsarray.append(d)
        else:
            for num in self.fitsid:
                d = self.fitsdir/(self.tag+str(num)+self.extension+'.fits')
                if check and not d.exists():
                    print(d, 'does not exist.')
                    if self.not_ignore_warning:
                        print(
                            'If you wanna proceed despite this warning, set self.not_ignore_warning = False in your fitsset.')
                        sys.exit('Error.')
                else:
                    if string:
                        fitsarray.append(str(d))
                    else:
                        fitsarray.append(d)
        return fitsarray

    def clean(self):
        """Clean i.e. remove fits files if exists."""
        import os
        if self.unlock:
            for i in range(len(self.path())):
                if self.path(check=False)[i].exists():
                    os.remove(self.path(check=False, string=True)[i])
                    print('rm old '+self.path(check=False, string=True)[i])

    def at_list(self, listname='tmp'):
        """make at list used in IRAF."""

        listname = listname+'.list'
        atlistname = str(self.fitsdir/listname)
        f = open(atlistname, mode='w')
        for fits in self.path(string=True):
            f.write(fits+'\n')
        f.close()
        atlistname = '@'+atlistname
        return atlistname

#    def imcombine(self,tag,combine="median"):
#        combined_fitsset=FitsSet(tag,self.fitsdir)
#        combined_fitsset.clean()
#        iraf.imcombine(input=self.at_list(listname=tag),output=combined_fitsset.path(check=False)[0],combine=combine) #
#        return combined_fitsset
#
#    def apall(self,apref,wavref=None,extout=".ec",extwout="_1dw"):
#        """apall
#        Args:
#           apref: reference aperture
#        """
#
#        currentdir=os.getcwd()
#        os.chdir(str(self.fitsdir))
#        apref_path=apref.path(string=False,check=True)[0].name #ref_ap
#        iraf.imred()
#        iraf.eche()
#        extoned=[]
#        for i,eachpath in enumerate(tqdm.tqdm(self.path(string=False,check=True))):
#            outfile=(eachpath.name).replace(".fits",extout+".fits")
#            extoned.append(outfile)
#            iraf.apall(input=eachpath.name,output=outfile,find="n",recenter="n",resize="n",edit="n",trace="n",fittrace#="n",extract="y",references=apref_path,review="n",interactive="n")
#
#
#        if wavref != None:
#            #wavref_path=wavref.path(string=False,check=True)[0].name #wav ref
#            #CHECKING ecfile in database
#            wavref_path=directory_util.cp(self.fitsdir,wavref,"ec")
#            for i,eachpath in enumerate(tqdm.tqdm(self.path(string=False,check=True))):
#                outwfile=(eachpath.name).replace(".fits",extwout+".fits")
#                print("WAV",extoned[i],wavref_path)
#                iraf.refs(input=extoned[i],references=wavref_path.name,select="match")
#                iraf.dispcor(input=extoned[i],output=outwfile,flux="no")
#
#
#        os.chdir(currentdir)
#        return

#    def apnormalize(self,tag,ref):
#        currentdir=os.getcwd()
#        os.chdir(anadir)
#
#        apnormalfits=FitsSet(tag,self.fitsdir)
#        iraf.imred.echell.apnormaize(input=self.path(check=False,string=True)[0],output=apnormalfits.path(check=False,string=True)[0],find="n",recenter="n",resize="n",edit="y",trace="n",fittrace="n",references=ref,order=35,niterat=1)
#
#        os.chdir(currentdir)
#        return apnormalfits
