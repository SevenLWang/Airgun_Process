"""
Created on Tue Jul 18 14:35:00 2020

Modify History:
July 20, 2020 Li Wang compiled the Pre_Process part, from deconstructing to Rotatio.
August 23, 2020 Li Wang changed the type of instrument response removal (vel to none) and the unit (cm to m) for the usage of kiwi tool.
September 8, 2020 Li Wang removed the component rotation part for maintain the NEZ coordinate.
September 9, 2020 Li Wang added the part of file format conversion part (convert format SAC to Miniseed), and modularization.


@author: alexwang
"""

import subprocess
import os
from glob import glob
from os.path import join, basename
import sys
from obspy import UTCDateTime
import numpy as np
import pandas as pd
from tqdm import tqdm
import shutil
from pyrocko import io

###############################################################################

###############################################################################

year = sys.argv[1]
min_day = sys.argv[2]
max_day = sys.argv[3]

dire = '/work/wang_li/Project/Airgun_Process'
middle_save = join(dire,'Middle_Save')
mseed_dir = join(dire,'mseed_file')
mseed_file_dir = join(mseed_dir,year+'*')
sac_dir = join(dire,'SAC_data')
# sac_year_dir = join(sac_dir,year)
yeardir = join(sac_dir,year+'*')
sac_file_dir = join(yeardir,'*.SAC')

os.environ['SAC_DISPLAY_COPYRIGHT']='0'


###############################################################################

    # deconstructe the original records 

def ReattributeMiniSEED():
    ### Reattribute the mseed records
    cmd = ''
    for day in range(1,366):
        if int(day) >= int(min_day) and int(day) <= int(max_day):
            day = str(day).zfill(3)
            date = year + day
            d = join(sac_dir,date)
            #print(d)
            cmd += 'mkdir %s\n' % (d)
        #    print(cmd)
    os.system(cmd)


    print('Reattribute the original mseed files according to the date')
    for f in tqdm(glob(join(mseed_file_dir,'*.mseed'))):
        filename = basename(f)
        year = filename.split('.')[5]
        day = filename.split('.')[6]
        for d in glob(yeardir):
            yy = basename(d)[0:4]
            dd = basename(d)[4:7]
            if (year == yy and day == dd):
                cmd = 'cp {} {}'.format(f,d)
                os.system(cmd)

    ### Convert the mseed records into SAC files
    print('Convert mseed files into SAC files:')
    cmd = ''
    for d in tqdm(glob(yeardir)):
        cmd1 = ''
        for f in glob(join(d,'*.mseed')):
            cmd1 += 'mseed2sac %s\n' % (f)
        os.system(cmd1)
        cmd += 'rm %s/*.mseed\n ' % (d)
    os.system(cmd)

    for f in glob(join(dire,'*.SAC')):
        filename = basename(f)
        year = filename.split('.')[5]
        day = filename.split('.')[6]
        for d in glob(yeardir):
            yy = basename(d)[0:4]
            dd = basename(d)[4:7]
            if (year == yy and day == dd):
                cmd = 'mv {} {}'.format(f,d)
                os.system(cmd)

###############################################################################

    # Pre_Process: rmean,rtrend,taper and transfer the instrument response

def PreProcess():
    print('Pre_Process of sac files:')
    for f in tqdm(glob(sac_file_dir)):
        p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        s = ''
        s += 'r *.SAC\n'
        s += 'rmean\n'
        s += 'rtrend\n'
        s += 'taper\n'
        s += 'trans from evalresp fname /work/wang_li/Project/Airgun_Process/RESP.ALL to none freq 0.004 0.005 30 32\n'
        s += 'mul 1e-9\n'
        s += 'w over\n'
        s += 'q\n' 
        p.communicate(s.encode())

###############################################################################

    # write the header information into SAC files

def HeaderInfoWrite():
    stations = pd.read_table('/work/wang_li/Project/Airgun_Process/meta.G3.V20191013',names=['Net','Sta','Loc','Chan','Lat','Lon','Elev','Depth','Az','Inc','Inst','Scale','ScaleFreq','ScaleUnits','SampleRate','Start','End'],encoding='gb2312',sep=',')
    print('write the header information into SAC files')

    for f in tqdm(glob(join(sac_file_dir))):
    #    print(f)
        sacname = os.path.basename(f)
        temp = sacname.split('.')
        for i in range(len(stations)):
            net = stations['Net'][i]
            station = stations['Sta'][i]
            loc = stations['Loc'][i]
            cha = stations['Chan'][i]
            lon = stations['Lon'][i]
            lat = stations['Lat'][i]
    #        print(net,station,loc,cha,lon,lat)
    #        print(temp[0],temp[1],temp[2])
            if (net == temp[0] and station == temp[1] and int(loc) == int(temp[2]) and cha == temp[3]):
    #            print(lat,lon)
                p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                s = ""
                s += "r %s\n" % f
                s += "ch stla %.6f\n" % float(lat)
                s += "ch stlo %.6f\n" % float(lon)
                s += "ch evlo 100.098 evla 38.7485\n"
                s += "wh\n"
                s += "q\n"
                p.communicate(s.encode())
            else:
                continue

##############################################################################

    # cut the selected event(min_day:max_day) according to the excitation timetable 

def CutEvents():
    events = pd.read_table(join(dire, 'otime_2015_2018'), index_col=False, names=['Year','Julday','Hour','Minute','Second','MSecond','','Number','Coeffcient'], encoding='gb2312', sep=' ')
    print('cut the selected event(min_day:max_day) according to the excitation timetable')

    for d in tqdm(glob(yeardir)):
        date = basename(d)
        cmd = 'mkdir %s/%s\n' % (middle_save,date)
        os.system(cmd)
        for j in range(1,41):
            input_data = glob(join(d, '*G3.ZDY%02d.*' % (j)))
            for f in (input_data):
                number = 1
                filename = basename(f)
                #print(filename)
                y = filename.split('.')[5];jd = filename.split('.')[6];net = filename.split('.')[0];sta = filename.split('.')[1];loc = filename.split('.')[2];cha = filename.split('.')[3]
                #print(y,jd,net,sta,loc,cha)
                for i in range(len(events)):
                    year = events['Year'][i]
                    jday = events['Julday'][i]
                    hour = events['Hour'][i]
                    minute = events['Minute'][i]
                    second = events['Second'][i]
                    msecond = events['MSecond'][i]
                    #print(year,y,jday,jd)
                    if int(year) == int(y) and int(jday) == int(jd) :
                        cmd = "saclst nzyear nzjday nzhour nzmin nzsec nzmsec f %s" % (f)
                        junk, kzyy, kzjd, kzhou, kzmin, kzsec, kzmsec = os.popen(cmd).read().split()
                        Time = UTCDateTime(year=int(year),julday=int(jday),hour=int(hour),minute=int(minute),second=int(second),microsecond=(int(msecond)*1000))
                        kztime = UTCDateTime(year=int(kzyy),julday=int(kzjd),hour=int(kzhou),minute=int(kzmin),second=int(kzsec),microsecond=(int(kzmsec)*1000))
                        otime = Time - kztime 
                        btime = otime - 100
                        etime = otime + 300
                        #print(otime,btime,etime)
                        p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                        s = " "
                        s += "cut %.3f %.3f\n" % (float(btime),float(etime))
                        s += "r %s\n" % (f)
                        s += "ch b %.3f e %.3f\n" % (float(btime),float(etime))
                        s += "ch o %.3f\n" % (float(otime))
                        s += "ch allt (0 - &1,o&) iztype IO\n"
                        s += "w %d%d.%s.%s.%s.%s.%s.SAC\n" % (int(y),int(jd),str(number).zfill(3),net,sta,loc,cha)
                        s += "q\n"
                        p.communicate(s.encode())
                        number += 1
                    else:
                        continue

    cmd = ''
    for f in glob(join(dire,'*.SAC')):
        filename = basename(f)
        date = filename.split('.')[0]
        for d in glob(join(middle_save,'*')):
            dname = basename(d)
            if (date == dname):
                cmd = 'mv %s %s\n' % (f,d)
                os.system(cmd)
            else:
                continue

##############################################################################

    # Convert the file format from SAC into MiniSEED, and Reattribution & Rename

def FileFormatConversion():
    print('Convert the file format from SAC into MiniSEED')
    for f in tqdm(glob(join(middle_save,year+'*','*.SAC'))):
        traces = io.load(f, format='sac')
        out_filename = f[:-4] + '.mseed'
        io.save(traces, out_filename)

    print('Reattribute the cut waveforms according to the excitation time:')
    for d in glob(join(middle_save,year+'*')):
        sequences = []
        for f in glob(join(d,'*.SAC')):
            filename = basename(f)
            sequence = filename.split('.')[1]
            sequences.append(sequence)
        maxexcitation = max(sequences)
        cmd = ''
        for i in range(int(maxexcitation)):
            newdir = 'No.' + str(i+1).zfill(3)
            newdir = join(d,newdir)
            cmd += 'mkdir %s\n' % (newdir)
        os.system(cmd)

    cmd = ''
    for d in glob(join(middle_save,year+'*')):
        date = basename(d)
        for dd in glob(join(d,'No.*'))
            dirname = basename(dd)
            seq = dirname.split('.')[1]
            for f in tqdm(glob(join(d,'*.mseed'))):
                filename = basename(f)
                datetime = filename.split('.')[0]
                sequence = filename.split('.')[1]
                if int(datetime) == int(date) and int(sequence) == int(seq):
                    cmd = 'cp %s %s' % (f,dd)
                    os.system(cmd)
                else:
                    continue

    print('Rename the MiniSEED files for compatible Kiwi Tool use')
    for d in tqdm(glob(join(middle_save,year+'*','No.*'))):
        for f in glob(join(d,'*.mseed')):
            filename = basename(f)
            staname = filename.split('.')[3]
            chaname = filename.split('.')[5]
            chaname = 'B' + chaname[1:3]
            name = ['DISPL',staname,chaname]
            name = '.'.join(name)
            path = os.path.realpath(d)
            os.rename(f,join(path,name))

    for d in glob(join(middle_save,year+'*'))
        cmd = 'rm %s/*.SAC' % (d)
        cmd += 'rm %s/*.mseedd' % (d)
        os.system(cmd)

###############################################################################

    # components rotation

def ComponentRotation():
    print('Rotate the components from ENZ to RTZ:')
    ###check the number of components in every station and delete the stations which loss at least one component
    for d in glob(join(middle_save,year+'*')):
        filenames = glob(join(d,'*.SAC'))
        for num in range(1,41):
            n = 0
            for files in (filenames):
                sta = basename(files).split(".")[5]; check = 'ZDY%02d' % (num)
                if sta == check :
                    n += 1
            if n > 0 and n < 3 :
                cmd = "rm %s/*ZDY%02d*" % (d,num) 
                os.system(cmd)
            
    ### Rotation
    for input_dire in tqdm(glob(join(middle_save,year+'*'))):
        number = 0
        for filename in os.listdir(input_dire):
            if os.path.splitext(filename)[1]=='.SAC':
                number += 1
        time = number//3
        SAC = subprocess.Popen(['sac'],stdin=subprocess.PIPE) 
        s = " "
        output_dire = os.system("mkdir %s/%s_rot/" % (input_dire,basename(input_dire)))
        input_dire1=glob(join(input_dire,"*.SH?.SAC"))
        input_dire1.sort()
        for freq in range(time):
            file_1 = input_dire1[freq*3]
            file_2 = input_dire1[freq*3+1]
            file_3 = input_dire1[freq*3+2]
            filename = os.path.basename(file_1)
            net = filename.split(".")[2]
            sta = filename.split(".")[3]
            loc = filename.split(".")[4]
            date = filename.split(".")[0]
            number = filename.split(".")[1]
            year = date[0:4]
            jday = date[4:7]

            #path1 = file_1.split('/')[-1]
            #path2 = file_2.split('/')[-1]
            #path3 = file_3.split('/')[-1]

    ### Cmpaz Correction
            s += "r %s\nch cmpaz %8.2f\nwh\n" % (file_1, 90)
            s += "r %s\nch cmpaz %8.2f\nwh\n" % (file_2, 0)
            s += "r %s\nch cmpaz %8.2f\nwh\n" % (file_3, 0)

    ### Cmpinc Correction
            s += "r %s %s\n ch cmpinc 90\n wh\n" % (file_1,file_2)
            s += "r %s\n ch cmpinc 0\n wh\n" % (file_3)

    ### Header Write    
            s += "r %s %s\n" % (file_1,file_2)
            s += "rot to gcp\n"
            s += "w %s %s\n" % (join(year+jday+'.'+number+'.'+net+'.'+sta+'.'+loc+'.r'), join(year+jday+'.'+number+'.'+net+'.'+sta+'.'+loc+'.t'))
            s += "cut off\n"
            shutil.copy(file_3, join(year+jday+'.'+number+'.'+net+'.'+sta+'.'+loc+'.z'))
        s += "q\n"
        SAC.communicate(s.encode())
        os.system("mv *.[rtz] %s/%s_rot/" % (input_dire,basename(input_dire)))

###############################################################################

    # main program

ReattributeMiniSEED()
PreProcess()
HeaderInfoWrite()
CutEvents()
FileFormatConversion()
# ComponentRotation()
