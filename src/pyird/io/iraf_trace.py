"""Load IRAF-made aperture file

"""
import numpy as np

def read_trace_file(filename):
    """
    Args:
       trace file

    Returns:
        y0, interp_function, xmin, xmax, coeff

    """
    
    with open(filename) as f:
        cont = f.readlines()
    f.close()
    norder=0
    y0=[]
    read_curve=False
    interp_function=[]
    curvepar=[]
    xmin=[]
    xmax=[]
    coeff=[]
    for line in cont:
        arr=line.split()
        if len(arr)>0:
            #identify begin
            if arr[0]=="begin":
                #reset read_curve
                if len(curvepar)>0:
                    interp_function.append(curvepar[0])
                    nlegendreorder=int(curvepar[1])
                    xmin.append(curvepar[2])
                    xmax.append(curvepar[3])
                    coeffval=[float(s) for s in curvepar[4:4+nlegendreorder]]
                    coeff.append(coeffval)
                    read_curve=False
                
                norder=norder+1
                y0.append(float(arr[4]))
            elif arr[0]=="curve":
                read_curve=True
                curvepar=[]
            elif read_curve and arr[0]!="#":
                curvepar.append(float(arr[0]))
                
    interp_function.append(float(curvepar[0]))
    xmin.append(float(curvepar[2]))
    xmax.append(float(curvepar[3]))
    nlegendreorder=int(curvepar[1])
    coeffval=[float(s) for s in curvepar[4:4+nlegendreorder]]
    coeff.append(coeffval)
    interp_function=np.array(interp_function,dtype=int)
    xmin=np.array(xmin,dtype=int)
    xmax=np.array(xmax,dtype=int)
    
    return y0, interp_function, xmin, xmax, coeff


if __name__=="__main__":
    filename="/home/kawahara/pyird/data/iraf/aprefA"
    read_trace_file(filename)

