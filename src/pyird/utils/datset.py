"""Set of dat files."""

import sys
__all__ = ['DatSet']


class DatSet(object):

    def __init__(self, dir, prefix='', extension=''):
        self.prefix = prefix
        self.extension = extension
        self.dir = dir
        self.fitsid = None
        self.unlock = True  # unlock for cleanning (i.e. remove)
        self.not_ignore_warning = True

    def path(self, string=False, check=True):
        """get path array.

        Returns:
            array of paths
        """
        datarray = []
        if self.fitsid is not None:
            for num in self.fitsid:
                d = self.dir/(self.prefix+str(num)+self.extension+'.dat')
                if check and not d.exists():
                    print(d, 'does not exist.')
                    if self.not_ignore_warning:
                        print(
                            'If you wanna proceed despite this warning, set self.not_ignore_warning = False in your fitsset.')
                        sys.exit('Error.')
                else:
                    if string:
                        datarray.append(str(d))
                    else:
                        datarray.append(d)
        return datarray

    def flatpath(self, string=False, check=True):
        """get path for the normalized flat.

        Returns:
            path for the normalized flat
        """

        d = self.dir/('nwflat_'+self.band+self.extension+'.dat')
        if check and not d.exists():
            print(d, 'does not exist.')
            if self.not_ignore_warning:
                print(
                    'If you wanna proceed despite this warning, set self.not_ignore_warning = False in your fitsset.')
                sys.exit('Error.')
        else:
            if string:
                d = str(d)
        return d

    def clean(self):
        """Clean i.e. remove files if exists."""
        import os
        if self.unlock:
            for i in range(len(self.path())):
                if self.path(check=False)[i].exists():
                    os.remove(self.path(check=False, string=True)[i])
                    print('rm old '+self.path(check=False, string=True)[i])
