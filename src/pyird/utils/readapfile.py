def readap(apfile):
    print(apfile)
    f = open(apfile, 'r')
    norder = 0
    ap = []
    for line in f:
        #        print(repr(line))
        if line[0] != '#':
            sline = (line.split('\t'))
            if sline[0] == 'begin':
                norder = norder+1
                tmp = []
                tmp.append(line)
                if norder > 1:
                    ap.append(tmp)
            else:
                tmp.append(line)
    ap.append(tmp)
    f.close()
    print('number of orders=', norder)
    return ap


def writeap(ap, outfile):
    f = open(outfile, 'w')
    for apf in ap:
        for eachline in apf:
            f.write(eachline)
    f.close()


if __name__ == '__main__':
    ap = readap('aptest')
    writeap(ap, 'testapw')
