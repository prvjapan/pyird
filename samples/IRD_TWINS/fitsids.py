import pathlib
from pyird import irdstream
from pyird import getinfo

def get_coordinate():
    from astropy.coordinates import EarthLocation
    subaru = EarthLocation.of_site('Subaru')
    ra,dec = getinfo.get_radec("HD 111959")    
    RV=-1.64 #
    info=[ra,dec,RV,'2021-3-17']
        
    return subaru,info
    
def directories(mode="YJ"):
    datadir=pathlib.Path("/media/kawahara/kingyo/IRD_TWINS/data")
    flatdatadir=pathlib.Path("/media/kawahara/kingyo/IRD_G196_3B/data")

    if mode=="YJ":
        print("YJ band")
        anadir=pathlib.Path("/home/kawahara/ird/pyird/samples/IRD_TWINS/ana/YJ")
    elif mode=="H":
        print("H band")
        anadir=pathlib.Path("/home/kawahara/ird/pyird/samples/IRD_TWINS/ana/H")
    else:
        print("no mode")

    return datadir,anadir,flatdatadir

def IRDBD_Stream2D(mode="YJ"):
    datadir,anadir,flatdatadir=directories(mode=mode)
    
    #Define Stream2D from raw data
    targets=irdstream.Stream2D("targets",datadir,anadir)
    targets.fitsid=\
    list(range(42306,42324,2))\
    +list(range(42548,42550,2))
    
    flat_mmf=irdstream.Stream2D("flat_mmf",flatdatadir,anadir)
    flat_mmf.fitsid=\
    list(range(41704,41904,2))

#    thar strong
    thars1=irdstream.Stream2D("thars1",datadir,anadir)
    thars1.fitsid=list(range(42364,42440,2)) #THARSTAR
    
    if mode=="YJ":
        print("YJ band")
    elif mode=="H":
        print("H band")
        targets.fitsid_increment()
        flat_mmf.fitsid_increment()
        thars1.fitsid_increment()
        
    return targets, flat_mmf, thars1


def IRDBD_Stream1D(mode="YJ"):
    datadir,anadir,fdatadir=directories(mode=mode)
    #Define Stream1D from raw data

    SKY=irdstream.Stream1D("SKY",datadir,anadir)
    SKY.fitsid=[]
    
    HD111959A=irdstream.Stream1D("HD111959A",datadir,anadir)
    HD111959A.fitsid=[41308]
    HD111959B=irdstream.Stream1D("HD111959B",datadir,anadir)
    HD111959B.fitsid=[41310]
    TwoMassJ1249A=irdstream.Stream1D("TwoMassJ1249A",datadir,anadir)
    TwoMassJ1249A.fitsid=[41312]
    TwoMassJ1249B=irdstream.Stream1D("TwoMassJ1249B",datadir,anadir)
    TwoMassJ1249B.fitsid=[41314]
    TwoMassJ1323A=irdstream.Stream1D("TwoMassJ1323A",datadir,anadir)
    TwoMassJ1323A.fitsid=[41316]
    TwoMassJ1323B=irdstream.Stream1D("TwoMassJ1323B",datadir,anadir)
    TwoMassJ1323B.fitsid=[41318]
    TwoMassJ1047A=irdstream.Stream1D("TwoMassJ1047A",datadir,anadir)
    TwoMassJ1047A.fitsid=[41320]
    TwoMassJ1047B=irdstream.Stream1D("TwoMassJ1047B",datadir,anadir)
    TwoMassJ1047B.fitsid=[41322]

    if mode=="YJ":
        print("YJ band")
    elif mode=="H":
        print("H band")
        HD111959A.fitsid_increment()
        HD111959B.fitsid_increment()
        TwoMassJ1249A.fitsid_increment()
        TwoMassJ1249B.fitsid_increment()
        TwoMassJ1323A.fitsid_increment()
        TwoMassJ1323B.fitsid_increment()
        TwoMassJ1047A.fitsid_increment()
        TwoMassJ1047B.fitsid_increment()
        
    return HD111959A,HD111959B,TwoMassJ1249A,TwoMassJ1249B,TwoMassJ1323A,TwoMassJ1323B,TwoMassJ1047A,TwoMassJ1047B
