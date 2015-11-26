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

def parsing_input_file(fname,mode='linear-elastic'):
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
        1   --> number of layer
        2   --> soil layer parameters
    """
    with open(fname) as f:
        hl = []; vs = []; rho = []; qs = []
        tfPair = []
        data = []
        layerID=1; TFID=1
        if mode.lower()==cd.mode[0] or mode.lower()==cd.mode[1]:
            data.append([cd.mode[0]])
            blockID = 1
            for line in f:
                if not comment_check(line):
                    if blockID==1:
                        ntf = blockreader(blockID,line)
                        data.append([ntf])
                        blockID+=1
                    elif blockID==2:
                        if TFID<=ntf:
                            stfPair = blockreader(blockID,line)
                            tfPair.append(stfPair)
                            TFID+=1
                            if TFID>ntf:
                                data.append(tfPair)
                                blockID+=1
                    elif blockID==3:
                        nl = blockreader(blockID,line)
                        data.append([nl])
                        blockID+=1
                    elif blockID==4:
                        if layerID<=nl:
                            shl,svs,srho,sqs = blockreader(blockID,line)
                            hl.append(shl)
                            vs.append(svs)
                            rho.append(srho)
                            qs.append(9999. if mode.lower()==cd.mode[0] else sqs)
                            layerID+=1
                            if layerID>nl:
                                data.append([hl,vs,rho,qs])
                                blockID+=1
            return data
        elif mode.lower()==cd.mode[2]:
            print('not implemented! Sorry!')
        elif mode.lower()==cd.mode[3]:
            print('not implemented! Sorry!')
        else:
            raise ValueError('mode is not supported!')
            
def comment_check(stringinput):
    return True if stringinput[0]=='#' else False
    
def blockreader(blockID,line):
    if blockID==1:
        return int(line)
    elif blockID==2:
        s1,s2 = line.split()
        return [int(s1)-1,int(s2)-1]
    elif blockID==3:
        return int(line)
    elif blockID==4:
        shl,svs,srho,sqs = line.split()
        return float(shl),float(svs),float(srho),float(sqs)