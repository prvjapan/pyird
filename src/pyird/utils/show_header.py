#!/usr/bin/env python
import argparse
from astropy.io import fits
import numpy as np


def loadtermcol():
    tc = ['\033[94m', '\033[91m', '\033[0m']
    return tc


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', nargs='+', required=True, help='fits file')
    parser.add_argument('-t', nargs='+', help='tag(s)')
    parser.add_argument(
        '-r', help='unique header for -t option', action='store_true')
    parser.add_argument('-s', help='show sequences', action='store_true')
    args = parser.parse_args()

    if args.t:
        print('FITS', end=' ')
        for tag in args.t:
            print(',', tag, end=' ')
        print('')

    tc = loadtermcol()
    unp = '-'
    filep = 'Extract unique header'
    i = 0
    j = 0
    for file in args.f:
        if not args.r:
            print(file, end=' ')

        hduread = fits.open(file)
        if args.t:
            un = ''
            for tag in args.t:
                if args.r:
                    un = un+' '+str(hduread[0].header[tag])
                else:
                    print(',', hduread[0].header[tag], end=' ')

            if un != unp and args.r:
                j = j+1
                if i > 0:
                    print('('+objp+')')
                    if args.s:
                        print('-------------')
                    print(filep, unp)

                print(tc[np.mod(j, len(tc))]+file, un)
                unp = un
                i = 0
            elif file == args.f[-1] and args.r:
                print('('+objp+')')
                print(tc[np.mod(j, len(tc))]+file, un)
            elif args.r:
                i = i+1
                if args.s:
                    if i == 1:
                        print('-------------')
                        print('List: ',  filep.replace('.fits', '').replace(
                            'q', '').replace('v', ''), end=' ')
                    print(file.replace('.fits', '').replace(
                        'q', '').replace('v', ''), end=' ')
                else:
                    print('.', end=' ')

            else:
                print('')
        else:
            ha = hduread[0].header
            print(repr(ha))
        filep = file
        objp = hduread[0].header['OBJECT']
