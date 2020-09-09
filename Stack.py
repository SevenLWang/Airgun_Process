"""
Created on Tue Jul 14 09:08:50 2020

Modify History:
July 14, 2020 Li Wang compiled the stacking part, from reattributing to final stacking
July 17, 2020 Li Wang added deconvolution part 
August 23, 2020 Li Wang changed the naming method of stacked records for adapting the kiwi tool

@author: alexwang
"""

import subprocess
import os
from glob import glob
from os.path import join, basename
import sys
import numpy as np
import obspy
import pandas as pd
from tqdm import tqdm

############################################################################### 

###############################################################################

year = sys.argv[1]
stack = sys.argv[2]
day = ['290-292','291-301','293-312']
if int(year) == 2015:
    date = year + day[0]
elif int(year) == 2016:
    date = year + day[1]
elif int(year) == 2017: 
    date = year + day[2]

ddir = '/work/wang_li/Project/Airgun_Process'
dire = join(ddir,'Stack_Process')
yeardir = '%s*' % (year)
rotdir = '%s*_rot' % (year)
targetdir = join(dire,'20stacks_prepared')
low_snr_dir = join(dire,'Low_SNR_Records')
final_stack_dir = join(dire,'Final_Stack_Records')
differ_dir = join(targetdir,date)
differdir = glob(join(differ_dir,'ZDY*'))

dire_deconv = join(ddir,'Deconvolution_Process')
refer_dir = join(dire_deconv,'Refer_Records')
result_dir = join(dire_deconv,'Result_Records')

os.environ['SAC_DISPLAY_COPYRIGHT']='0'

############################################################################### 

############################################################################### 

    #reattribute waveform records based on the channel

cmd = 'mkdir %s %s %s %s %s' % (targetdir,low_snr_dir,final_stack_dir,refer_dir,result_dir)
os.system(cmd)
cmd = 'mkdir %s/%s %s/%s %s/%s %s/%s' % (targetdir,date,final_stack_dir,date,refer_dir,date,result_dir,date)
os.system(cmd)

    
for num in range(1,41):
    num = str(num).zfill(2)
    d1 = targetdir +'/ZDY'+ num + '.r'
    d2 = targetdir +'/ZDY'+ num + '.t'
    d3 = targetdir +'/ZDY'+ num + '.z'
    d4 = differ_dir + '/ZDY'+ num + '.r'
    d5 = differ_dir + '/ZDY'+ num + '.t'
    d6 = differ_dir + '/ZDY'+ num + '.z'
    #print(d)
    cmd = 'mkdir {} {} {} {} {} {}' .format(d1,d2,d3,d4,d5,d6)
    os.system(cmd)

print('Reattribute records:')
for f in tqdm(glob(join(ddir,'Middle_Save',yeardir,rotdir,'*.[rtz]'))):
   filename = basename(f)
   station = filename.split('.')[3]
   channel = filename.split('.')[5]
   for d in glob(join(targetdir,'ZDY*')):
       ss = basename(d).split('.')[0]
       cc = basename(d).split('.')[1]
       if (station == ss and channel == cc):
           cmd = 'cp {} {}'.format(f,d)
           os.system(cmd)
           #print(f)
           break

for d in glob(join(targetdir,'ZDY*')):
    sta = basename(d).split('.')[0]
    cmd = "ls -l %s|grep '%s'|wc -l" % (d,sta)
    filenumber = os.popen(cmd).read()
    if int(filenumber) == 0:
        cmd = 'rm -r {}'.format(d)
        os.system(cmd)
        continue
    else:
        continue

############################################################################## 

############################################################################## 
    
    #stack every 20 channels record

p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
print('Stack every 20-channel records:')
for d in tqdm(glob(join(targetdir,'ZDY*'))):
    #    print(d)
    count = 0
    for root,dirs,files in os.walk(d):
        for f in files:
            if os.path.splitext(f)[1] == '.r' or os.path.splitext(f)[1] == '.t' or os.path.splitext(f)[1] == '.z':
                count += 1
    if count < int(stack):
        max_rounds = 1
    else:
        max_rounds = count//int(stack)
    ss = basename(d).split('.')[0]
    cc = basename(d).split('.')[1]
    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
    stacknumber = 0
    rounds = 0
    fdir = glob(join(d,'*.[rtz]'))
    fdir.sort()
    for root,dirs,files in os.walk(d):
        for f in files:
            if os.path.splitext(f)[1] == '.r' or os.path.splitext(f)[1] == '.t' or os.path.splitext(f)[1] == '.z':
                count += 1
    if count >= int(stack):
        max_rounds = count//int(stack)
    #    print(max_rounds)
        for f in (fdir):
            # open the sac stack program 
            if rounds < max_rounds:
                if stacknumber == 0:
                    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                    s = ''
                    s += 'sss\n'
                    s += 'zerostack\n'
                    cmd = 'saclst evlo evla stlo stla f %s ' % (f)
                    junk,evlo,evla,stlo,stla = os.popen(cmd).read().split()
                    s += 'addstack %s\n' % (f)
                    stacknumber += 1
                    continue
                else:
                    # continue the stack
                    if stacknumber < int(stack):
                        s += 'addstack %s\n' % (f)
                        stacknumber += 1
                    # output the stacked record
                    else:
                        rounds += 1
                        name = ['linestack',date,str(stacknumber),str(rounds).zfill(2),ss,cc]
                        name = '.'.join(name)
                        s += 'timewindow -20 180\n'
                        s += 'sumstack\n'
                        s += 'writestack %s\n' % (name)
                        s += 'quitsub\n'
                        s += 'r %s\n' % (name)
                        s += 'ch evlo %.6f evla %.6f stlo %.6f stla %.6f\n' % (float(evlo),float(evla),float(stlo),float(stla))
                        if int(year) == 2015:
                            s += 'ch nzyear 2015 nzjday 290 nzhour 00 nzmin 00 nzsec 00 nzmsec 000\n'
                        elif int(year) == 2016:
                            s += 'ch nzyear 2016 nzjday 291 nzhour 00 nzmin 00 nzsec 00 nzmsec 000\n'
                        elif int(year) == 2017:
                            s += 'ch nzyear 2017 nzjday 293 nzhour 00 nzmin 00 nzsec 00 nzmsec 000\n'
                        s += 'ch o 20\n'
                        s += 'wh\n'
                        s += 'q\n'
                        p.communicate(s.encode())
                        stacknumber = 0
                    continue
            else:
                s += 'quit\n'
                p.communicate(s.encode())

    cmd = 'rm %s/*.[rtz]' % (d)
    os.system(cmd)
############################################################################### 

############################################################################### 

    #select the data according to the SNR

stations = pd.read_table('/work/wang_li/Project/Airgun_Process/arrival_time.dat',index_col=False, names=['Station','P_arrival','Delta'], encoding='gb2312', sep=' ')

files = glob(join(dire,'linestack.*'))
files.sort()
#print(files)  
print('Calculate the Signal-Noise ratio of every stacked record:')  
for f in tqdm(files):

### read the data   
#    print(f)
    filename = basename(f)
    station = filename.split('.')[4]
    channel = filename.split('.')[5]
    stack = filename.split('.')[2]
    st = obspy.read(f)
    tr = st[0]
    tr_filt = tr.copy()
    tr_filt.filter('bandpass',freqmin = 3, freqmax = 8 ,corners = 4, zerophase = True)

### match the data    
    for i in range(len(stations)):
        staname = stations['Station'][i]
        Ptime = stations['P_arrival'][i]; Ptime = float(Ptime)
        dt = stations['Delta'][i]; dt = float(dt)
        sta = staname.split('.')[2]
        if station == sta:
            signal_winlen = 30; noise_winlen = 20
            s_index1 = int(Ptime/dt); s_index2 = int((Ptime + signal_winlen)/dt)
            n_index1 = int((Ptime-noise_winlen)/dt); n_index2 = int(Ptime/dt)
            signal = tr_filt[s_index1:s_index2]
            noise = tr_filt[n_index1:n_index2]

### cauculate the SNR
            max_signal = max(signal)
            std_noise = np.std(noise)
            ratio = max_signal/std_noise
            #ratio = 10 * math.log(ratio,10)
            #snr.append(ratio)
            p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            s = ''
            s += 'r %s\n' % (f)
            s += 'ch user0 %.4f\n' % (float(ratio))
            s += 'wh\n'
            s += 'q\n'
            p.communicate(s.encode())
            break
        else:
            continue

### select the high SNR record
print('Select the high Signal-Noise ratio records:') 
for f in tqdm(glob(join(dire,'linestack.*'))):
    cmd = "saclst user0 f %s | awk '{print $2}'" % (f)
    snr = os.popen(cmd).read().strip()
    snr = float(snr)
    if snr >= 5.0:         
        sta = basename(f).split('.')[4]
        cha = basename(f).split('.')[5]
        for d in glob(join(differ_dir,'ZDY*')):
            ss = basename(d).split('.')[0]
            cc = basename(d).split('.')[1]
            if (sta == ss and cha == cc):
                cmd = 'mv {} {}'.format(f,d)
                os.system(cmd)
    else:
        cmd = 'mv {} {}'.format(f,low_snr_dir)
        os.system(cmd)

### calculate the available stack number for every channel of each staion

cmd = 'rm %s_stacknumber.dat\n' % (date)
os.system(cmd)
print('Write the total stack number of every station:')
differdir.sort()        
for d in tqdm(differdir):
    number = 0
    for filename in glob(join(d,'*.[rtz]')):
        name = basename(filename).split('.')[2]
        number += int(name)
    staname = basename(d)
    cmd = 'echo %s %d >> %s_stacknumber.dat' % (staname,number,date)
    os.system(cmd)
    
###############################################################################        

###############################################################################  

    # stack the selected records

print('Stack the selected high Signal-Noise ratio records')
for d in tqdm(differdir):
    ss = basename(d).split('.')[0]
    cc = basename(d).split('.')[1]
    cmd = "ls -l %s|grep '%s'|wc -l" % (d,ss)
    filenumber = os.popen(cmd).read()
    stacknumber = 0
    if int(filenumber) != 0:
        p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        s = ''
        s += 'sss\n'
        s += 'zerostack\n'
        for f in glob(join(d,'linestack.*')):
            number = int(basename(f).split('.')[2])
            print(number,f)
            cmd = 'saclst evlo evla stlo stla f %s ' % (f)
            junk,evlo,evla,stlo,stla = os.popen(cmd).read().split()
            s += 'addstack %s\n' % (f)
            stacknumber = stacknumber + number
#        name = ['linestack',date,'final',str(stacknumber),ss,cc]
        if cc == 'r':
            cc = 'BHR'
        elif cc == 't':
            cc = 'BHT'
        elif cc == 'z':
            cc = 'BHZ'
        name = ['DISPL',ss,cc]
        name = '.'.join(name)
        s += 'timewindow 0 200\n'
        s += 'sumstack\n'
        s += 'writestack %s\n' % (name)
        s += 'quitsub\n'
        s += 'r %s\n' % (name)
#        s += 'bp c 2 10 n 4 p 2'
        s += 'ch evlo %.6f evla %.6f stlo %.6f stla %.6f\n' % (float(evlo),float(evla),float(stlo),float(stla))
        if int(year) == 2015:
            s += 'ch nzyear 2015 nzjday 290 nzhour 00 nzmin 00 nzsec 00 nzmsec 000\n'
        elif int(year) == 2016:
            s += 'ch nzyear 2016 nzjday 291 nzhour 00 nzmin 00 nzsec 00 nzmsec 000\n'
        elif int(year) == 2017:
            s += 'ch nzyear 2017 nzjday 293 nzhour 00 nzmin 00 nzsec 00 nzmsec 000\n'
        s += 'ch o 20\n'
        s += 'wh\n'
        s += 'q\n'
        p.communicate(s.encode())

cmd = 'mv DISPL.* %s' % (final_stack_dir)
os.system(cmd)

for f in glob(join(final_stack_dir,'DISP*')):
    st = obspy.read(f)
    tr = st[0]
    delta = tr.stats.delta
    time = 0.0
    filename = basename(f)
    with open(filename,'w') as table:
        for i in range(20001):
            line = ''
            time = float(delta) * i
            line += str(time) + ' ' + str(tr.data[i]) + ' ' +'\n'
            table.write(line)

for f in glob('DISP*'):
    filename = basename(f)
    name1 = filename.split('.')[1]
    name2 = filename.split('.')[2]
    name = ['DISPL',name1,name2]
    name = '.'.join(name)
    os.rename(filename,name)

############################################################################## 

##############################################################################

    deconvolute the final stacked records

### cut the signal window of records as the reference waveform in deconvolution from station ZDY22 
files = join(final_stack_dir,'*.ZDY22.*')
for f in glob(files):
    filename = basename(f)
    stacktype = filename.split('.')[0]
    number = filename.split('.')[3]
    station = filename.split('.')[4]
    channel = filename.split('.')[5]
    p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    s = ''
    s += 'cut 19 23\n'
    s += 'r %s\n' % (f)
    s += 'w %s.%s.%s.%s.refer.%s\n' % (stacktype,date,number,station,channel)
    s += 'q\n'
    p.communicate(s.encode())
cmd = 'mv *.[rtz] %s/%s' % (refer_dir,date)
os.system(cmd)

### deconvolution
print('Deconvolution:')
files = glob(join(final_stack_dir,'*.[rtz]'))
refer_files = glob(join(refer_dir,date,'*.[rtz]'))
for f in tqdm(files):
    filename = basename(f)
    stacktype = filename.split('.')[0]
    date = filename.split('.')[1]
    number = filename.split('.')[3]
    station = filename.split('.')[4]
    channel = filename.split('.')[5]
    name = [stacktype,date,station,'deconv',channel]
    name = '.'.join(name)
    for refer_file in (refer_files):
        STACKTYPE = basename(refer_file).split('.')[0]
        DATE = basename(refer_file).split('.')[1]
        NUMBER = basename(refer_file).split('.')[3]
        STATION = basename(refer_file).split('.')[4]
        CHANNEL = basename(refer_file).split('.')[5]
        if date == DATE and channel == CHANNEL:
            cmd = '/work/wang_li/Project/Airgun_Process/wtdeconv -s %s -d %s -w 0.05 -o %s' % (f,refer_file,name)
            os.system(cmd)
            p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            s = ''
            s += 'r %s\n' % (name)
#            s += 'bp c 3 8 n 4 p 2\n'
            s += 'ch KSTNM {} KCMPNM {} KNETWK G3\n'.format(station,channel)
            s += 'w over\n'
            s += 'q\n'
            p.communicate(s.encode())
            break
        else:
            continue

cmd = 'mv *.[rtz] %s/%s' % (result_dir,date)
os.system(cmd)

### supply the header information of records
stainfor_dir = join(ddir,'station_infor.dat')
stations = pd.read_table(stainfor_dir,names=['Stations','CMPAZ','CMPINC'],encoding='gb2312',sep=' ')

print('Supply some header information:')
for f in tqdm(glob(join(result_dir,date,'*.[rtz]'))):
    sacname = basename(f)
    temp = sacname.split('.')
    for i in range(len(stations)):
        station_name = stations['Stations'][i]
        cmpaz = stations['CMPAZ'][i]
        cmpinc = stations['CMPINC'][i]
        temp1 = station_name.split('.')
        if (temp[3] == temp1[3] and temp[5] == temp1[5]):
            p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            s = ""
            s += "r %s\n" % (f)
            s += "ch cmpaz %.6f\n" % float(cmpaz)
            s += "ch cmpinc %.6f\n" % float(cmpinc)
            s += "wh\n"
            s += "q\n"
            p.communicate(s.encode())
            break
        else:
            continue
                    
###############################################################################  
        
        
        


