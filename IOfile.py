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
    data['parfile'] = fname
    data['sensitifity'] = False
    with open(fname) as f:
        
        for line in f:
            if not comment_check(line):
                mode = blockreader(1,line)        
                mode = mode.strip()
                data['mode'] = mode
                break
        
        modeID = cd.mode.index(mode)
        data['modeID'] = modeID
        hl = []; vs = []; dn = []; qs = []
        vp = []; qp = [];
        tfPair = []; soiltype = [];
        layerID=1; TFID=1
        #if mode.lower()==cd.mode[0] or mode.lower()==cd.mode[1] or \
        #    mode.lower()==cd.mode[2] or mode.lower()==cd.mode[3] or\
        #    mode.lower()==cd.mode[4] or mode.lower()==cd.mode[5] or \
        #    mode.lower()==cd.mode[6] or mode.lower()==cd.mode[7]:

        if modeID<13:
            if modeID<5:
                blockSeq = [1,2,11,3,4,5,6,7,8]
            elif modeID<11:
                blockSeq = [1,2,11,3,4,5,6,7,9]
            elif modeID==11:
                blockSeq = [1,2,11,10,3,4,5,6,7,12]
            elif modeID==12:
                blockSeq = [1,2,11,10,3,4,5,6,7,13]
            IDSeq = 1
            for line in f:
                if not comment_check(line):
                    if blockSeq[IDSeq]==2:
                        data['inputmotion'] = blockreader(blockSeq[IDSeq],line)
                        IDSeq+=1
                    elif blockSeq[IDSeq]==11:
                        data['inputtype'] = blockreader(2,line)
                        IDSeq+=1
                    elif blockSeq[IDSeq]==10:
                        data['GoverGmaxfile'] = blockreader(2,line)
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
                            qp.append(9999. if modeID==5 or modeID==7 else sqp)
                            qs.append(9999. if modeID==5 or modeID==7 else sqs)
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
                    elif blockSeq[IDSeq]==12:
                        if layerID<=data['nlayer']:
                            shl,svs,sdn,sqs,st = blockreader(blockSeq[IDSeq],line)
                            hl.append(shl)
                            vs.append(svs)
                            dn.append(sdn)
                            qs.append(9999. if modeID==0 or modeID==2 else sqs)
                            soiltype.append(int(st))
                            layerID+=1
                            if layerID>data['nlayer']:
                                data['hl']=hl
                                data['vs']=vs
                                data['dn']=dn
                                data['qs']=qs
                                data['soiltype']=soiltype
                                IDSeq+=1
                    elif blockSeq[IDSeq]==13:
                        if layerID<=data['nlayer']:
                            shl,svp,svs,sdn,sqp,sqs,st = blockreader(blockSeq[IDSeq],line)
                            hl.append(shl)
                            vp.append(svp)
                            vs.append(svs)
                            dn.append(sdn)
                            qp.append(9999. if modeID==5 or modeID==7 else sqp)
                            qs.append(9999. if modeID==5 or modeID==7 else sqs)
                            soiltype.append(int(st))
                            layerID+=1
                            if layerID>data['nlayer']:
                                data['hl']=hl
                                data['vp']=vp
                                data['vs']=vs
                                data['dn']=dn
                                data['qp']=qp
                                data['qs']=qs
                                data['comp']=mode[4].lower()
                                data['soiltype'] = soiltype
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
    elif blockID==12:
        shl,svs,srho,sqs,soiltype = line.split()
        return float(shl),float(svs),float(srho),float(sqs),int(soiltype)
    elif blockID==13:
        shl,svp,svs,srho,sqp,sqs,soiltype = line.split()
        return float(shl),float(svp),float(svs),float(srho),float(sqp),float(sqs),int(soiltype)
        
def parsing_nonlinear_parameter(fname,verbose=False):
    """
    function to read non linear parameter
    file format (ascii):
        "number of soil types --> N"
        "number of strain on soil type I"
        "strain"   "G/Gmax"  "Damping Ratio"
        "number of strain on soil type II"
        "strain"   "G/Gmax"  "Damping Ratio"
        ...
        "number of strain on soil type N"
        "strain"   "G/Gmax"  "Damping Ratio"
    """
    
    with open(fname) as f:
        lineID = 0
        nonlinpar = {}
        for line in f:
            if lineID==0:
                nonlinpar['N soil type']=int(line)
                nonlinpar['n-samples']=[]
                nonlinpar['nonlin strain']=[[] for _ in range(nonlinpar['N soil type'])]
                nonlinpar['nonlin G/Gmax']=[[] for _ in range(nonlinpar['N soil type'])]
                nonlinpar['nonlin damping']=[[] for _ in range(nonlinpar['N soil type'])]
                lineID += 1
                lineSoilType = 0
            elif lineID==1:
                nonlinpar['n-samples'].append(int(line))
                lineID += 1
                lineSoilLine = 1
            elif lineID==2:
                tmp = line.split()
                nonlinpar['nonlin strain'][lineSoilType].append(float(tmp[0]))
                nonlinpar['nonlin G/Gmax'][lineSoilType].append(float(tmp[1]))
                nonlinpar['nonlin damping'][lineSoilType].append(float(tmp[2]))
                if lineSoilLine<nonlinpar['n-samples'][-1]:
                    lineSoilLine += 1
                else:
                    lineSoilLine = 1
                    if lineSoilType+1==nonlinpar['N soil type']:
                        break
                    else:
                        lineSoilType += 1
                        lineID = 1
        if verbose:
            print('')
            print('Number of non linear set of soil profile : %i'%nonlinpar['N soil type'])
            for i in range(nonlinpar['N soil type']):
                print('Soil ID : %2i'%i)
                print('\tMin/Max strain \t: %1.5f and %1.5f'%(min(nonlinpar['nonlin strain'][i]),max(nonlinpar['nonlin strain'][i])))
                print('\tMin/Max G/Gmax \t: %1.5f and %1.5f'%(min(nonlinpar['nonlin G/Gmax'][i]),max(nonlinpar['nonlin G/Gmax'][i]))) 
                print('\tMin/Max Damping ratio \t: %2.4f and %2.4f'%(min(nonlinpar['nonlin damping'][i]),max(nonlinpar['nonlin damping'][i])))
            print('')
        return nonlinpar

def read_ascii_seismogram(fname):
    """
    read seismogram with ascii format
    'time' 'amplitude (m/s^2)'
    """
    time = []
    amp = []
    with open(fname) as f:
        for line in f:
            tmp = line.split()
            time.append(float(tmp[0]))
            amp.append(float(tmp[1]))
    return time,amp