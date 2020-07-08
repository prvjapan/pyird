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

    def check_ready(self):
        check_fitsset(self.comparison)
        check_fitsset(self.targets)
        check_fitsset(self.flat_mmf1)
        check_fitsset(self.flat_mmf2)
            
def check_fitsset(obj):
    print("#:",len(obj["num"]))
    for num in obj["num"]:
        d=obj["dir"]/(obj["tag"]+str(num)+".fits")
        if not d.exists():
            print(d, "does not exist.")
    
