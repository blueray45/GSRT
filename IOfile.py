# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 19:02:32 2015
The whole code to manage I/O in GSRT

Created by :
Theodosius Marwan Irnaka
marwan.irnaka@gmail.com
GNU Lisence
"""

import commondata as cd

def parsing_input_file(fname):
    """
    Simple subroutine to parse input file, based on the calculation mode
    
    USAGE :
    
    parseddata = parsing_input_file('xxxx.dat',mode)
    
    mode -->
        'linear-elastic'
        'linear-viscoelastic'
        'linear-equivalent'
        'non-linear'
    """
    
    """
    input block numbering :
        1   --> calculation mode
        2   --> input motion
        3   --> angle of incidence
        4   --> location of source
        5   --> number of TF
        6   --> TF pairs
        7   --> number of layer
        8   --> soil layer parameters
    """

    # reading the whole configuration
    data = {}
    with open(fname) as f:
        
        for line in f:
            if not comment_check(line):
                mode = blockreader(1,line)        
                mode = mode.strip()
                data['mode'] = mode
                break
        
        modeID = cd.mode.index(mode)
        
        hl = []; vs = []; dn = []; qs = []
        vp = []; qp = [];
        tfPair = []
        layerID=1; TFID=1
        #if mode.lower()==cd.mode[0] or mode.lower()==cd.mode[1] or \
        #    mode.lower()==cd.mode[2] or mode.lower()==cd.mode[3] or\
        #    mode.lower()==cd.mode[4] or mode.lower()==cd.mode[5] or \
        #    mode.lower()==cd.mode[6] or mode.lower()==cd.mode[7]:
        
        if modeID<8:
            if modeID<4:
                blockSeq = [1,2,3,4,5,6,7,8]
            elif modeID<8:
                blockSeq = [1,2,3,4,5,6,7,9]
            IDSeq = 1
            for line in f:
                if not comment_check(line):
                    if blockSeq[IDSeq]==2:
                        data['inputmotion'] = blockreader(blockSeq[IDSeq],line)
                        IDSeq+=1
                    elif blockSeq[IDSeq]==3:
                        data['iang'] = blockreader(blockSeq[IDSeq],line)
                        IDSeq+=1
                    elif blockSeq[IDSeq]==4:
                        data['sourceloc'] = blockreader(blockSeq[IDSeq],line)
                        IDSeq+=1
                    elif blockSeq[IDSeq]==5:
                        data['ntf'] = blockreader(blockSeq[IDSeq],line)
                        IDSeq+=1
                    elif blockSeq[IDSeq]==6:
                        if TFID<=data['ntf']:
                            stfPair = blockreader(blockSeq[IDSeq],line)
                            tfPair.append(stfPair)
                            TFID+=1
                            if TFID>data['ntf']:
                                data['tfPair']=tfPair
                                IDSeq+=1
                    elif blockSeq[IDSeq]==7:
                        data['nlayer'] = blockreader(blockSeq[IDSeq],line)
                        IDSeq+=1
                    elif blockSeq[IDSeq]==8:
                        if layerID<=data['nlayer']:
                            shl,svs,sdn,sqs = blockreader(blockSeq[IDSeq],line)
                            hl.append(shl)
                            vs.append(svs)
                            dn.append(sdn)
                            qs.append(9999. if modeID==0 or modeID==2 else sqs)
                            layerID+=1
                            if layerID>data['nlayer']:
                                data['hl']=hl
                                data['vs']=vs
                                data['dn']=dn
                                data['qs']=qs
                                IDSeq+=1
                    elif blockSeq[IDSeq]==9:
                        if layerID<=data['nlayer']:
                            shl,svp,svs,sdn,sqp,sqs = blockreader(blockSeq[IDSeq],line)
                            hl.append(shl)
                            vp.append(svp)
                            vs.append(svs)
                            dn.append(sdn)
                            qp.append(9999. if modeID==4 or modeID==6 else sqp)
                            qs.append(9999. if modeID==4 or modeID==6 else sqs)
                            layerID+=1
                            if layerID>data['nlayer']:
                                data['hl']=hl
                                data['vp']=vp
                                data['vs']=vs
                                data['dn']=dn
                                data['qp']=qp
                                data['qs']=qs
                                data['comp']=mode[4].lower()
                                IDSeq+=1
            return data
        #elif :
        #    print('not implemented! Sorry!')
        else:
            raise ValueError('mode is not supported!')
            
def comment_check(stringinput):
    return True if stringinput[0]=='#' else False
    
def blockreader(blockID,line):
    if blockID==1:
        return line.strip()
    elif blockID==2:
        return line.split()
    elif blockID==3:
        return float(line)*0.01745329252
    elif blockID==4:
        return int(line)-1
    elif blockID==5:
        return int(line)
    elif blockID==6:
        s1,s2 = line.split()
        return [int(s1)-1,int(s2)-1]
    elif blockID==7:
        return int(line)
    elif blockID==8:
        shl,svs,srho,sqs = line.split()
        return float(shl),float(svs),float(srho),float(sqs)
    elif blockID==9:
        shl,svp,svs,srho,sqp,sqs = line.split()
        return float(shl),float(svp),float(svs),float(srho),float(sqp),float(sqs)