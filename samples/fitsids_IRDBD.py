from pyird import irdstream
import pathlib

def get_coordinate():
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation
    from astropy import units as u

    lat_subaru=19.0 + 49/60.+ 43/3600.0
    lon_subaru=-(155.0 + 28/60.0 + 50/3600.0)
    height_subaru = 4139.0
    subaru = EarthLocation.of_site('Subaru')
    #    subaru = EarthLocation.from_geodetic(lat=lat_subaru*u.deg, lon=lon_subaru*u.deg, height=height_subaru*u.m)
    
    J2224ra=336.1827083
    J2224dec=-1.9814625
    J2224RV=-36.48
    J2224info=[J2224ra,J2224dec,J2224RV,'2020-7-22']

    J1645ra=251.3420
    J1645dec=-13.331189
    J1645RV=27.0
    J1645info=[J1645ra,J1645dec,J1645RV,'2020-7-22']
    
    J1507ra=226.9486250
    J1507dec=-16.4611439
    J1507RV=-39.85    
    J1507info=[J1507ra,J1507dec,J1507RV,'2020-7-13']
    
    return subaru,J2224info,J1507info,J1645info
    
def directories(mode="YJ"):
    datadir=pathlib.Path("/media/kawahara/kingyo/IRDBD/data")
    if mode=="YJ":
        print("YJ band")
        anadir=pathlib.Path("/media/kawahara/kingyo/IRDBD/ana/YJ")
    elif mode=="H":
        print("H band")
        anadir=pathlib.Path("/media/kawahara/kingyo/IRDBD/ana/H")
    else:
        print("no mode")

    return datadir,anadir

def IRDBD_Stream2D(mode="YJ"):
    datadir,anadir=directories(mode=mode)
    #Define Stream2D from raw data
    targets=irdstream.Stream2D("targets",datadir,anadir)
    targets.fitsid=\
    list(range(31062,31070,2))+list(range(31340,31346,2))+list(range(31346,31350,2))
    #J1507, J1645, J2224
    
    flat_mmf=irdstream.Stream2D("flat_mmf",datadir,anadir)
    flat_mmf.fitsid=\
    list(range(31580,31776,2))
#    list(range(31504,31546,2)) #THAR-STAR
#    list(range(31546,31576,2)) #THAR-COMB

#    thar strong
    thars1=irdstream.Stream2D("thars1",datadir,anadir)
    thars1.fitsid=list(range(31504,31540,2))+list(range(31542,31546,2))
    thars2=irdstream.Stream2D("thars2",datadir,anadir)
    thars2.fitsid=list(range(31546,31575,2))

    
    if mode=="YJ":
        print("YJ band")
    elif mode=="H":
        print("H band")
        targets.fitsid_increment()
        flat_mmf.fitsid_increment()
        thars1.fitsid_increment()
        thars2.fitsid_increment()
        
    return targets, flat_mmf, thars1, thars2


def IRDBD_Stream1D(mode="YJ"):
    datadir,anadir=directories(mode=mode)
    #Define Stream1D from raw data
    #J1507, J1645, J2224

    J1507=irdstream.Stream1D("J1507",datadir,anadir)
#    J1507.fitsid=list(range(31062,31070,2))
    J1507.fitsid=[31068]
#    J1507.fitsid=[31066,31068]

    J1645=irdstream.Stream1D("J1645",datadir,anadir)
#    J1645.fitsid=list(range(31340,31346,2))
    J1645.fitsid=[31344]

    J2224=irdstream.Stream1D("J2224",datadir,anadir)
    #    J2224.fitsid=list(range(31346,31350,2))
    J2224.fitsid=[31348]
    if mode=="YJ":
        print("YJ band")
    elif mode=="H":
        print("H band")
        J1507.fitsid_increment()
        J1645.fitsid_increment()
        J2224.fitsid_increment()
        
    return J1507, J1645, J2224

