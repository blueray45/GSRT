# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 02:06:13 2016

@author: irnakat
"""

bigmachine = False

# batch test transfer function
import datetime
import sys
import platform as pf
import time
import numpy as np

globalstart = time.clock()
fname = 'log_'+str(datetime.datetime.now().date())+'_'+str(datetime.datetime.now().time())+'.txt'
#fname = 'logtest.txt'
def get_processor_info():
    import platform, subprocess
    if platform.system() == "Windows":
        return platform.processor()
    elif platform.system() == "Darwin":
        return subprocess.check_output(['/usr/sbin/sysctl', "-n", "machdep.cpu.brand_string"]).strip()
    elif platform.system() == "Linux":
        command = "cat /proc/cpuinfo"
        return subprocess.check_output(command, shell=True).strip()
    return ""


if bigmachine:
    # for big machine
    ntest = 20
    inputfileid = [0,1,2,3,4,5,8,9]
    inputfileidpsv = [6,7]
    inputfileidlayers = [0,3]
    inputfileidlayerspsv = [6,7]
    nlayer = np.logspace(0,3,16)
else:
    # for small machine
    ntest = 1
    inputfileid = [0]
    inputfileidpsv = [7]
    inputfileidlayers = [0]
    inputfileidlayerspsv = [6]
    nlayer = np.logspace(0,2,11)

clf_error = 'sys.exit("Calculation is aborted! Check log file on "+fname+" !")'
with open(fname,'w+') as f:
    # 000 - write basic info on log 
    print('000 - basic info')
    f.write('000 - basic info\n')
    f.write('time \t : '+str(datetime.datetime.now())+'\n')
    f.write('num test\t : '+str(ntest)+'\n')
    f.write('platform\t : '+sys.platform+'\n')
    f.write('processor\t : '+str(pf.processor())+'\n')    
    try:
        t = get_processor_info()
        f.write('more proc\t : '+t[t.find('model name	: ')+13:t.find('stepping')-1]+'\n')  
    except:
        print 'more proc failed! move on.'
    f.write('Py build\t : '+str(pf.python_build())+'\n')  
    f.write('PyCompiler\t : '+str(pf.python_compiler())+'\n')
    f.write('Comp name\t : '+str(pf.uname()[1])+'\n')
    f.write('Kernel\t : '+str(pf.uname()[2])+'\n') 
    f.write('OS Version\t : '+str(pf.uname()[3])+'\n\n') 
    
    print(' --> OK')
    
    # 001 - check python modules
    print('001 - modules test')
    modulesList = ['numpy','copy','scipy','commondata','IOfile','GSRTtools',
                   'TFCalculator','TSCalculator','requests',
                   'os']
    try:
        for mid in modulesList:
            exec('import '+mid)
        f.write('001 - modules test\nModules test are complete!\n')
    except ImportError:
        f.write('001 - modules test\nModules "'+mid+'" is not exist!\nCalculation is aborted!\n')
        #raise("Calculation is aborted! Check log file on "+fname+" !")
        eval(clf_error)
    f.write('\n')
    
    print(' --> OK')
    
    #002 - Check Input File
    print('002 - input files')
    f.write('002 - input files\n')
    inputlist = ['sampleinput_linear_elastic_1layer_halfspace.dat',
                 'sampleinput_linear_elastic_6layer_halfspace.dat',
                 'sampleinput_linear_elastic_6layer_halfspace_400vs30.dat',
                 'sampleinput_linear_eq_elastic_1layer_halfspace.dat',
                 'sampleinput_linear_eq_elastic_6layer_halfspace.dat',
                 'sampleinput_linear_eq_elastic_6layer_halfspace_400vs30.dat',
                 'sampleinput_psv_s_linear_elastic_1layer_halfspace.dat',
                 'sampleinput_psv_s_linear_elastic_6layer_halfspace.dat',
                 'ferndale2_parameter.dat',
                 'ferndale2_parameter_lineq.dat']
    for iid in inputlist:
        if not os.path.isfile(iid):
            f.write('Input file --> "'+iid+'" is not exist!')
            eval(clf_error)
    f.write('Input files are complete!\n')
    
    #Import modules
    from TFCalculator import TFCalculator as TFC
    from TSCalculator import TSCalculator as TSC
    f.write('\n')
            
    print(' --> OK')
    
    #003 - Kramer Transfer Function
    print('003 - kramer TF')
    f.write('003 - kramer TF\n')
    f.write('fID\telapsedrd\tstdevrd\telapsed\tstdev\tsizetf\n')
    for iid in inputfileid:
        elapsedRead = []
        elapsedTF = []
        elapsedTFbytes = []
        for i in range(ntest):        
            start = time.clock()
            data = IOfile.parsing_input_file(inputlist[iid])
            elapsedRead.append(time.clock()-start)
            data['inputmotion']='dirac.dat'
            theclass = TFC(data)
            start = time.clock()
            theclass.tf_kramer286_sh()
            elapsedTF.append(time.clock()-start)
            elapsedTFbytes.append(numpy.array(theclass.tf).nbytes)
        f.write('%d\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t\n'%(iid,np.mean(elapsedRead),np.std(elapsedRead),
                                                    np.mean(elapsedTF),np.std(elapsedTF),np.max(elapsedTFbytes)))
    f.write('\n')
        
    print(' --> OK')
    
    #004 - Knopoff SH simple Transfer Function
    print('004 - knopoff SH simple TF')
    f.write('004 - knopoff SH simple TF\n')
    f.write('fID\telapsedrd\tstdevrd\telapsed\tstdev\tsizetf\n')
    for iid in inputfileid:
        elapsedRead = []
        elapsedTF = []
        elapsedTFbytes = []
        for i in range(ntest):        
            start = time.clock()
            data = IOfile.parsing_input_file(inputlist[iid])
            elapsedRead.append(time.clock()-start)
            data['inputmotion']='dirac.dat'
            theclass = TFC(data)
            start = time.clock()
            theclass.tf_knopoff_sh()
            elapsedTF.append(time.clock()-start)
            elapsedTFbytes.append(numpy.array(theclass.tf).nbytes)
        f.write('%d\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t\n'%(iid,np.mean(elapsedRead),np.std(elapsedRead),
                                                    np.mean(elapsedTF),np.std(elapsedTF),np.max(elapsedTFbytes)))
    f.write('\n')
        
    print(' --> OK')
    
    #005 - Knopoff SH advance Transfer Function
    print('005 - knopoff SH advance TF')
    f.write('005 - knopoff SH advance TF\n')
    f.write('fID\telapsedrd\tstdevrd\telapsed\tstdev\tsizetf\n')
    for iid in inputfileid:
        elapsedRead = []
        elapsedTF = []
        elapsedTFbytes = []
        for i in range(ntest):        
            start = time.clock()
            data = IOfile.parsing_input_file(inputlist[iid])
            elapsedRead.append(time.clock()-start)
            data['inputmotion']='dirac.dat'
            theclass = TFC(data)
            start = time.clock()
            theclass.tf_knopoff_sh_adv()
            elapsedTF.append(time.clock()-start)
            elapsedTFbytes.append(numpy.array(theclass.tf).nbytes)
        f.write('%d\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t\n'%(iid,np.mean(elapsedRead),np.std(elapsedRead),
                                                    np.mean(elapsedTF),np.std(elapsedTF),np.max(elapsedTFbytes)))
    f.write('\n')
        
    print(' --> OK')
    
    #006 - Kennet SH Transfer Function
    print('006 - Kennet SH TF')
    f.write('006 - Kennet SH TF\n')
    f.write('fID\telapsedrd\tstdevrd\telapsed\tstdev\tsizetf\n')
    for iid in inputfileid:
        elapsedRead = []
        elapsedTF = []
        elapsedTFbytes = []
        for i in range(ntest):        
            start = time.clock()
            data = IOfile.parsing_input_file(inputlist[iid])
            elapsedRead.append(time.clock()-start)
            data['inputmotion']='dirac.dat'
            theclass = TFC(data)
            start = time.clock()
            theclass.tf_kennet_sh()
            elapsedTF.append(time.clock()-start)
            elapsedTFbytes.append(numpy.array(theclass.tf).nbytes)
        f.write('%d\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t\n'%(iid,np.mean(elapsedRead),np.std(elapsedRead),
                                                    np.mean(elapsedTF),np.std(elapsedTF),np.max(elapsedTFbytes)))
    f.write('\n')
        
    print(' --> OK')
    
    #007 - Knopoff PSV Transfer Function
    print('007 - Knopoff PSV TF')
    f.write('007 - Knopoff PSV TF\n')
    f.write('fID\telapsedrd\tstdevrd\telapsed\tstdev\tsizetf\n')
    for iid in inputfileidpsv:
        elapsedRead = []
        elapsedTF = []
        elapsedTFbytes = []
        for i in range(ntest):        
            start = time.clock()
            data = IOfile.parsing_input_file(inputlist[iid])
            elapsedRead.append(time.clock()-start)
            data['inputmotion']='dirac.dat'
            theclass = TFC(data)
            start = time.clock()
            theclass.tf_knopoff_psv_adv()
            elapsedTF.append(time.clock()-start)
            elapsedTFbytes.append(numpy.array(theclass.tf).nbytes)
        f.write('%d\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t\n'%(iid,np.mean(elapsedRead),np.std(elapsedRead),
                                                    np.mean(elapsedTF),np.std(elapsedTF),np.max(elapsedTFbytes)))
    f.write('\n')
        
    print(' --> OK')
    
    #008 - Kramer SH Time Series - dirac.dat
    print('008 - Kramer SH TS')
    f.write('008 - Kramer SH TS\n')
    f.write('fID\tn-iter\telapsed\tstdev\n')
    method = 'kramer286_sh'
    for iid in inputfileid:
        elapsedTS = []
        for i in range(ntest):
            start = time.clock()
            tsclass = TSC(inputlist[iid],method=method)
            elapsedTS.append(time.clock()-start)
        f.write('%d\t%d\t%.2e\t%.2e\n'%(iid,tsclass.lastiter,np.mean(elapsedTS),np.std(elapsedTS)))        
    f.write('\n')
        
    print(' --> OK')
    
    #009 - Knopoff SH Simple Time Series - dirac.dat
    print('009 - Knopoff SH Simple TS')
    f.write('009 - Knopoff SH Simple TS\n')
    f.write('fID\tn-iter\telapsed\tstdev\n')
    method = 'knopoff_sh'
    for iid in inputfileid:
        elapsedTS = []
        for i in range(ntest):
            start = time.clock()
            tsclass = TSC(inputlist[iid],method=method)
            elapsedTS.append(time.clock()-start)
        f.write('%d\t%d\t%.2e\t%.2e\n'%(iid,tsclass.lastiter,np.mean(elapsedTS),np.std(elapsedTS)))        
    f.write('\n')
        
    print(' --> OK')
    
    #010 - Knopoff SH advance Time Series - dirac.dat
    print('010 - Knopoff SH advance TS')
    f.write('010 - Knopoff SH advance TS\n')
    f.write('fID\tn-iter\telapsed\tstdev\n')
    method = 'knopoff_sh_adv'
    for iid in inputfileid:
        elapsedTS = []
        for i in range(ntest):
            start = time.clock()
            tsclass = TSC(inputlist[iid],method=method)
            elapsedTS.append(time.clock()-start)
        f.write('%d\t%d\t%.2e\t%.2e\n'%(iid,tsclass.lastiter,np.mean(elapsedTS),np.std(elapsedTS)))        
    f.write('\n')
        
    print(' --> OK')
    
    #011 - Kennet SH Time Series - dirac.dat
    print('011 - Kennet SH TS')
    f.write('011 - Kennet SH TS\n')
    f.write('fID\tn-iter\telapsed\tstdev\n')
    method = 'kennet_sh'
    for iid in inputfileid:
        elapsedTS = []
        for i in range(ntest):
            start = time.clock()
            tsclass = TSC(inputlist[iid],method=method)
            elapsedTS.append(time.clock()-start)
        f.write('%d\t%d\t%.2e\t%.2e\n'%(iid,tsclass.lastiter,np.mean(elapsedTS),np.std(elapsedTS)))        
    f.write('\n')
        
    print(' --> OK')
    
    #012 - Knopoff PSV Time Series - dirac.dat
    print('012 - Knopoff PSV TS')
    f.write('012 - Knopoff PSV TS\n')
    f.write('fID\tn-iter\telapsed\tstdev\n')
    method = 'knopoff_psv_adv'
    for iid in inputfileidpsv:
        elapsedTS = []
        for i in range(ntest):
            start = time.clock()
            tsclass = TSC(inputlist[iid],method=method)
            elapsedTS.append(time.clock()-start)
        f.write('%d\t%d\t%.2e\t%.2e\n'%(iid,tsclass.lastiter,np.mean(elapsedTS),np.std(elapsedTS)))        
    f.write('\n')
        
    print(' --> OK')
    
    # extreme case test increase number of layer
    nlayer = [int(i) for i in nlayer]
    idanalysis = 12
    for nl in nlayer:
        #003 - Kramer Transfer Function
        idanalysis+=1
        print('%03d - Kramer SH %d layers'%(idanalysis,nl))
        f.write('%03d - Kramer SH %d layers\n'%(idanalysis,nl))
        f.write('fID\telapsedrd\tstdevrd\telapsed\tstdev\tsizetf\n')
        for iid in inputfileidlayers:
            elapsedRead = []
            elapsedTF = []
            elapsedTFbytes = []
            for i in range(ntest):        
                start = time.clock()
                data = IOfile.parsing_input_file(inputlist[iid])
                elapsedRead.append(time.clock()-start)
                data['inputmotion']='dirac.dat'
                data['sourceloc']=nl
                data['nlayer']=nl+1
                data['tfPair'][0][1]=nl
                newhl = []; newvs = []; newdn = []; newqs = []; newvp = []; newqp = []
                newsoiltype = []
                for j in range(nl):
                    newhl.append(data['hl'][0]/nl)
                    newvs.append(data['vs'][0])
                    newqs.append(data['qs'][0])
                    newdn.append(data['dn'][0])
                newhl.append(data['hl'][-1])
                newvs.append(data['vs'][-1])
                newqs.append(data['qs'][-1])
                newdn.append(data['dn'][-1])
                data['hl']=newhl
                data['vs']=newvs
                data['qs']=newqs
                data['dn']=newdn
                try:
                    for j in range(nl):
                        newvp.append(data['vp'])
                        newqp.append(data['qp'])
                    newvp.append(data['vp'][-1])
                    newqp.append(data['qp'][-1])
                    data['vp']=newvp
                    data['qp']=newqp
                except:
                    pass
                try:
                    for j in range(nl):
                        newsoiltype.append(data['soiltype'][0])
                    newsoiltype.append(data['soiltype'][-1])
                    data['soiltype']=newsoiltype
                except:
                    pass
                theclass = TFC(data)
                start = time.clock()
                theclass.tf_kramer286_sh()
                elapsedTF.append(time.clock()-start)
                elapsedTFbytes.append(numpy.array(theclass.tf).nbytes)
            f.write('%d\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t\n'%(iid,np.mean(elapsedRead),np.std(elapsedRead),
                                                        np.mean(elapsedTF),np.std(elapsedTF),np.max(elapsedTFbytes)))
        f.write('\n')
            
        print(' --> OK')
                      
        #004 - Knopoff SH simple Transfer Function
        idanalysis+=1
        print('%03d - Knopoff SH Simple %d layers'%(idanalysis,nl))
        f.write('%03d - Knopoff SH Simple %d layers\n'%(idanalysis,nl))
        f.write('fID\telapsedrd\tstdevrd\telapsed\tstdev\tsizetf\n')
        for iid in inputfileidlayers:
            elapsedRead = []
            elapsedTF = []
            elapsedTFbytes = []
            for i in range(ntest):        
                start = time.clock()
                data = IOfile.parsing_input_file(inputlist[iid])
                elapsedRead.append(time.clock()-start)
                data['inputmotion']='dirac.dat'
                data['sourceloc']=nl
                data['nlayer']=nl+1
                data['tfPair'][0][1]=nl
                newhl = []; newvs = []; newdn = []; newqs = []; newvp = []; newqp = []
                newsoiltype = []
                for j in range(nl):
                    newhl.append(data['hl'][0]/nl)
                    newvs.append(data['vs'][0])
                    newqs.append(data['qs'][0])
                    newdn.append(data['dn'][0])
                newhl.append(data['hl'][-1])
                newvs.append(data['vs'][-1])
                newqs.append(data['qs'][-1])
                newdn.append(data['dn'][-1])
                data['hl']=newhl
                data['vs']=newvs
                data['qs']=newqs
                data['dn']=newdn
                try:
                    for j in range(nl):
                        newvp.append(data['vp'])
                        newqp.append(data['qp'])
                    newvp.append(data['vp'][-1])
                    newqp.append(data['qp'][-1])
                    data['vp']=newvp
                    data['qp']=newqp
                except:
                    pass
                try:
                    for j in range(nl):
                        newsoiltype.append(data['soiltype'][0])
                    newsoiltype.append(data['soiltype'][-1])
                    data['soiltype']=newsoiltype
                except:
                    pass
                theclass = TFC(data)
                start = time.clock()
                theclass.tf_knopoff_sh()
                elapsedTF.append(time.clock()-start)
                elapsedTFbytes.append(numpy.array(theclass.tf).nbytes)
            f.write('%d\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t\n'%(iid,np.mean(elapsedRead),np.std(elapsedRead),
                                                        np.mean(elapsedTF),np.std(elapsedTF),np.max(elapsedTFbytes)))
        f.write('\n')
            
        print(' --> OK')
                   
        #005 - Knopoff SH advance Transfer Function
        idanalysis+=1
        print('%03d - Knopoff SH advance %d layers'%(idanalysis,nl))
        f.write('%03d - Knopoff SH advance %d layers\n'%(idanalysis,nl))
        f.write('fID\telapsedrd\tstdevrd\telapsed\tstdev\tsizetf\n')
        for iid in inputfileidlayers:
            elapsedRead = []
            elapsedTF = []
            elapsedTFbytes = []
            for i in range(ntest):        
                start = time.clock()
                data = IOfile.parsing_input_file(inputlist[iid])
                elapsedRead.append(time.clock()-start)
                data['inputmotion']='dirac.dat'
                data['sourceloc']=nl
                data['nlayer']=nl+1
                data['tfPair'][0][1]=nl
                newhl = []; newvs = []; newdn = []; newqs = []; newvp = []; newqp = []
                newsoiltype = []
                for j in range(nl):
                    newhl.append(data['hl'][0]/nl)
                    newvs.append(data['vs'][0])
                    newqs.append(data['qs'][0])
                    newdn.append(data['dn'][0])
                newhl.append(data['hl'][-1])
                newvs.append(data['vs'][-1])
                newqs.append(data['qs'][-1])
                newdn.append(data['dn'][-1])
                data['hl']=newhl
                data['vs']=newvs
                data['qs']=newqs
                data['dn']=newdn
                try:
                    for j in range(nl):
                        newvp.append(data['vp'])
                        newqp.append(data['qp'])
                    newvp.append(data['vp'][-1])
                    newqp.append(data['qp'][-1])
                    data['vp']=newvp
                    data['qp']=newqp
                except:
                    pass
                try:
                    for j in range(nl):
                        newsoiltype.append(data['soiltype'][0])
                    newsoiltype.append(data['soiltype'][-1])
                    data['soiltype']=newsoiltype
                except:
                    pass
                theclass = TFC(data)
                start = time.clock()
                theclass.tf_knopoff_sh_adv()
                elapsedTF.append(time.clock()-start)
                elapsedTFbytes.append(numpy.array(theclass.tf).nbytes)
            f.write('%d\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t\n'%(iid,np.mean(elapsedRead),np.std(elapsedRead),
                                                        np.mean(elapsedTF),np.std(elapsedTF),np.max(elapsedTFbytes)))
        f.write('\n')
            
        print(' --> OK')
                   
        #006 - Kennet SH Transfer Function
        idanalysis+=1
        print('%03d - Kennet SH %d layers'%(idanalysis,nl))
        f.write('%03d - Kennet SH %d layers\n'%(idanalysis,nl))
        f.write('fID\telapsedrd\tstdevrd\telapsed\tstdev\tsizetf\n')
        for iid in inputfileidlayers:
            elapsedRead = []
            elapsedTF = []
            elapsedTFbytes = []
            for i in range(ntest):        
                start = time.clock()
                data = IOfile.parsing_input_file(inputlist[iid])
                elapsedRead.append(time.clock()-start)
                data['inputmotion']='dirac.dat'
                data['sourceloc']=nl
                data['nlayer']=nl+1
                data['tfPair'][0][1]=nl
                newhl = []; newvs = []; newdn = []; newqs = []; newvp = []; newqp = []
                newsoiltype = []
                for j in range(nl):
                    newhl.append(data['hl'][0]/nl)
                    newvs.append(data['vs'][0])
                    newqs.append(data['qs'][0])
                    newdn.append(data['dn'][0])
                newhl.append(data['hl'][-1])
                newvs.append(data['vs'][-1])
                newqs.append(data['qs'][-1])
                newdn.append(data['dn'][-1])
                data['hl']=newhl
                data['vs']=newvs
                data['qs']=newqs
                data['dn']=newdn
                try:
                    for j in range(nl):
                        newvp.append(data['vp'])
                        newqp.append(data['qp'])
                    newvp.append(data['vp'][-1])
                    newqp.append(data['qp'][-1])
                    data['vp']=newvp
                    data['qp']=newqp
                except:
                    pass
                try:
                    for j in range(nl):
                        newsoiltype.append(data['soiltype'][0])
                    newsoiltype.append(data['soiltype'][-1])
                    data['soiltype']=newsoiltype
                except:
                    pass
                theclass = TFC(data)
                start = time.clock()
                theclass.tf_kennet_sh()
                elapsedTF.append(time.clock()-start)
                elapsedTFbytes.append(numpy.array(theclass.tf).nbytes)
            f.write('%d\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t\n'%(iid,np.mean(elapsedRead),np.std(elapsedRead),
                                                        np.mean(elapsedTF),np.std(elapsedTF),np.max(elapsedTFbytes)))
        f.write('\n')
            
        print(' --> OK')
                   
        #007 - Knopoff PSV Transfer Function
        idanalysis+=1
        print('%03d - Knopoff PSV %d layers'%(idanalysis,nl))
        f.write('%03d - Knopoff PSV %d layers\n'%(idanalysis,nl))
        f.write('fID\telapsedrd\tstdevrd\telapsed\tstdev\tsizetf\n')
        for iid in inputfileidlayerspsv:
            elapsedRead = []
            elapsedTF = []
            elapsedTFbytes = []
            for i in range(ntest):        
                start = time.clock()
                data = IOfile.parsing_input_file(inputlist[iid])
                elapsedRead.append(time.clock()-start)
                data['inputmotion']='dirac.dat'
                data['sourceloc']=nl
                data['nlayer']=nl+1
                data['tfPair'][0][1]=nl
                newhl = []; newvs = []; newdn = []; newqs = []; newvp = []; newqp = []
                newsoiltype = []
                for j in range(nl):
                    newhl.append(data['hl'][0]/nl)
                    newvs.append(data['vs'][0])
                    newqs.append(data['qs'][0])
                    newdn.append(data['dn'][0])
                newhl.append(data['hl'][-1])
                newvs.append(data['vs'][-1])
                newqs.append(data['qs'][-1])
                newdn.append(data['dn'][-1])
                data['hl']=newhl
                data['vs']=newvs
                data['qs']=newqs
                data['dn']=newdn
                try:
                    for j in range(nl):
                        newvp.append(data['vp'][0])
                        newqp.append(data['qp'][0])
                    newvp.append(data['vp'][-1])
                    newqp.append(data['qp'][-1])
                    data['vp']=newvp
                    data['qp']=newqp
                except:
                    pass
                try:
                    for j in range(nl):
                        newsoiltype.append(data['soiltype'][0])
                    newsoiltype.append(data['soiltype'][-1])
                    data['soiltype']=newsoiltype
                except:
                    pass
                theclass = TFC(data)
                start = time.clock()
                theclass.tf_knopoff_psv_adv()
                elapsedTF.append(time.clock()-start)
                elapsedTFbytes.append(numpy.array(theclass.tf).nbytes)
            f.write('%d\t%.2e\t%.2e\t%.2e\t%.2e\t%d\t\n'%(iid,np.mean(elapsedRead),np.std(elapsedRead),
                                                        np.mean(elapsedTF),np.std(elapsedTF),np.max(elapsedTFbytes)))
        f.write('\n')
            
        print(' --> OK')

    globalelapsed = (time.clock()-globalstart)
    f.write('Total Elapsed Time : %.4f'%globalelapsed)
print('Total elapsed time : %.4f'%globalelapsed)
    
#print ('\nreading log file\n')
#with open(fname,'r') as f:
#    print f.read()