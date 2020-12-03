"""
Created on Tue Jul 18 14:35:00 2020

Modify History:
July 20, 2020 Li Wang compiled the Pre_Process part, from deconstructing to Rotatio.
August 23, 2020 Li Wang changed the type of instrument response removal (vel to none) and the unit (cm to m) for the usage of kiwi tool.
September 8, 2020 Li Wang removed the component rotation part for maintain the NEZ coordinate.
September 9, 2020 Li Wang added the part of file format conversion part (convert format SAC to Miniseed), and modularization.
October 13, 2020 Li Wang added the part of stacking function for improve the SNR level of signals.
November 15, 2020 Li Wang adjusted the program for offering the waveforms adapted to Grond program.
November 18, 2020 Li Wang modified the wrong process of the pre-process part.
November 21, 2020 Li Wang modified the Pre-Processed part and separated the FileConversion part into two parts for Stack and NoStack


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
middle_for_stack = join(middle_save,'middledir_for_stack')
mseed_dir = join(dire,'mseed_file')
mseed_file_dir = join(mseed_dir,year+'*')
sac_dir = join(dire,'SAC_data')
# sac_year_dir = join(sac_dir,year)
yeardir = join(sac_dir,year+'???')
sac_file_dir = join(yeardir,'*.SAC')

os.environ['SAC_DISPLAY_COPYRIGHT']='0'


###############################################################################

    # deconstructe the original records 

def ReattributeMiniSEED():
    ## Reattribute the mseed records
    cmd = ''
    for day in range(1,366):
        if int(day) >= int(min_day) and int(day) <= int(max_day):
            for d in glob(join(mseed_dir,year+'*')):
                dd = basename(d)[4:7]
                if int(day) == int(dd):      
                    day = str(day).zfill(3)
                    date = year + day
                    d = join(sac_dir,date)
                    #print(d)
                    cmd += 'mkdir %s\n' % (d)
            # print(cmd)
            os.system(cmd)


    print('Reattribute the original mseed files according to the date')
    for d in glob(yeardir):
        yy = basename(d)[0:4]
        dd = basename(d)[4:7]
        if int(yy) == int(year) and int(dd) >= int(min_day) and int(dd) <= int(max_day):
            for f in tqdm(glob(join(mseed_file_dir,'*.mseed'))):
                filename = basename(f)
                y = filename.split('.')[5]
                day = filename.split('.')[6]
                if (y == yy and day == dd):
                    cmd = 'cp {} {}'.format(f,d)
                    os.system(cmd)

    ### Convert the mseed records into SAC files
    print('Convert mseed files into SAC files:')
    cmd = ''
    for d in glob(join(sac_dir,year+'???')):
        y = basename(d)[0:4]
        jd = basename(d)[4:7]
        if int(y) == int(year) and int(jd) >= int(min_day) and int(jd) <= int(max_day):
            cmd1 = ''
            for f in glob(join(d,'*.mseed')):
                cmd1 += 'mseed2sac %s\n' % (f)
            os.system(cmd1)
            cmd += 'rm %s/*.mseed\n ' % (d)
        os.system(cmd)

    for f in glob(join(dire,'*.SAC')):
        filename = basename(f)
        y = filename.split('.')[5]
        day = filename.split('.')[6]
        for d in glob(yeardir):
            yy = basename(d)[0:4]
            dd = basename(d)[4:7]
            if (y == yy and day == dd):
                cmd = 'mv {} {}'.format(f,d)
                os.system(cmd)

###############################################################################

    # write the header information into SAC files

def HeaderInfoWrite():
    stations = pd.read_table('/work/wang_li/Project/Airgun_Process/meta.G3.V20191013',names=['Net','Sta','Loc','Chan','Lat','Lon','Elev','Depth','Az','Inc','Inst','Scale','ScaleFreq','ScaleUnits','SampleRate','Start','End'],encoding='gb2312',sep=',')
    print('write the header information into SAC files')

    for d in glob(join(sac_dir,year+'???')):
        y = basename(d)[0:4]
        jd = basename(d)[4:7]
        if int(y) == int(year) and int(jd) >= int(min_day) and int(jd) <= int(max_day):
            for f in tqdm(glob(join(d,'*.SAC'))):
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
        y = basename(d)[0:4]
        jd = basename(d)[4:7]
        if int(y) == int(year) and int(jd) >= int(min_day) and int(jd) <= int(max_day):
            date = basename(d)
            cmd = 'mkdir %s/%s\n' % (middle_save,date)
            os.system(cmd)
            for j in range(1,41):
                input_data = glob(join(d, '*G3.ZDY%02d.*' % (j)))
                for f in (input_data):
                    number = 1
                    filename = basename(f)
                    # print(filename)
                    y = filename.split('.')[5];jd = filename.split('.')[6];net = filename.split('.')[0];sta = filename.split('.')[1];loc = filename.split('.')[2];cha = filename.split('.')[3]
                    #print(y,jd,net,sta,loc,cha)
                    for i in range(len(events)):
                        yy = events['Year'][i]
                        jday = events['Julday'][i]
                        hour = events['Hour'][i]
                        minute = events['Minute'][i]
                        second = events['Second'][i]
                        msecond = events['MSecond'][i]
                        # print(year,y,jday,jd)
                        if int(yy) == int(y) and int(jday) == int(jd) :
                            cmd = "saclst nzyear nzjday nzhour nzmin nzsec nzmsec f %s" % (f)
                            junk, kzyy, kzjd, kzhou, kzmin, kzsec, kzmsec = os.popen(cmd).read().split()
                            Time = UTCDateTime(year=int(y),julday=int(jday),hour=int(hour),minute=int(minute),second=int(second),microsecond=(int(msecond)*1000))
                            kztime = UTCDateTime(year=int(kzyy),julday=int(kzjd),hour=int(kzhou),minute=int(kzmin),second=int(kzsec),microsecond=(int(kzmsec)*1000))
                            otime = Time - kztime 
                            btime = otime - 20
                            etime = otime + 180
                            #print(otime,btime,etime)
                            p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                            s = " "
                            s += "cut %.3f %.3f\n" % (float(btime),float(etime))
                            s += "r %s\n" % (f)
                            s += "ch b %.3f e %.3f\n" % (float(btime),float(etime))
                            s += "ch o %.3f\n" % (float(otime))
                            s += "ch allt (0 - &1,b&) iztype IB\n"
                            # s += "interpolate delta 0.05\n"
                            s += "w %d%03d.%s.%s.%s.%s.%s.SAC\n" % (int(yy),int(jd),str(number).zfill(3),net,sta,loc,cha)
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

###############################################################################

    # Pre_Process: rmean,rtrend,taper and transfer the instrument response

def PreProcess_and_Reattribution():
    print('Pre_Process of SAC files:')
    for d in glob(join(middle_save,year+'???')):
        y = basename(d)[0:4]
        jd = basename(d)[4:7]
        if int(y) == int(year) and int(jd) >= int(min_day) and int(jd) <= int(max_day):
            for f in tqdm(glob(join(d,'*.SAC'))):
                path,files = os.path.split(f)
                newname,extension = os.path.splitext(f)
                newfile = join(path, newname+'.SAC_RESP')
                p = subprocess.Popen(['sac'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                s = ''
                s += 'r %s\n' % (f)
                s += 'rmean\n'
                s += 'rtrend\n'
                s += 'taper\n'
                s += 'trans from evalresp fname /work/wang_li/Project/Airgun_Process/RESP.ALL to none freq 0.004 0.005 40 42\n'
                # s += 'interpolate delta 0.05\n'
                s += 'w %s\n' % (newfile)
                s += 'q\n' 
                p.communicate(s.encode())

    print('Reattribute the processed waveforms according to the excitation time:')
    for d in glob(join(middle_save,year+'???')):
        y = basename(d)[0:4]
        jd = basename(d)[4:7]
        if int(y) == int(year) and int(jd) >= int(min_day) and int(jd) <= int(max_day):
            sequences = []
            for f in glob(join(d,'*.SAC_RESP')):
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
            for dd in glob(join(d,'No.*')):
                dirname = basename(dd)
                seq = dirname.split('.')[1]
                for f in glob(join(d,'*.SAC_RESP')):
                    filename = basename(f)
                    # datetime = filename.split('_')[5].split('-')[0] + filename.split('_')[10].split('.')[0]
                    # sequence = filename.split('_')[9]
                    sequence = filename.split('.')[1]
                    if int(sequence) == int(seq):
                        cmd = 'cp %s %s' % (f,dd)
                        os.system(cmd)
                    else:
                        continue

##############################################################################

    # Stack the single shot file, including file format conversion and rename

def Stack():
    cmd = 'mkdir {}'.format(middle_for_stack)
    os.system(cmd)

    for i in range(int(min_day),int(max_day)+1,1):
        d1 = join(middle_for_stack,year+str(i))
        d2 = join(middle_for_stack,year+str(i))
        d3 = join(middle_for_stack,year+str(i))
        cmd = 'mkdir {} {} {}' .format(d1,d2,d3)
        os.system(cmd)
        for num in range(1,41):
            num = str(num).zfill(2)
            dd1 = join(middle_for_stack,year+str(i),'ZDY'+num+'.SHE')
            dd2 = join(middle_for_stack,year+str(i),'ZDY'+num+'.SHN')
            dd3 = join(middle_for_stack,year+str(i),'ZDY'+num+'.SHZ')
            #print(d)
            cmd = 'mkdir {} {} {}' .format(dd1,dd2,dd3)
            os.system(cmd)

    print('Reattribute all single shots with the same channel into one directory :')
    for f in tqdm(glob(join(dire,'*.SAC'))):
        filename = basename(f)
        datetime = filename.split('.')[0]
        station = filename.split('.')[3]
        channel = filename.split('.')[5]
        for d in glob(join(middle_for_stack,'201*')):
            dt = basename(d)
            if (datetime == dt):
                for dd in glob(join(d,'ZDY*')):
                    ss = basename(dd).split('.')[0]
                    cc = basename(dd).split('.')[1]
                    if (station == ss and channel == cc):
                        cmd = 'mv {} {}'.format(f,dd)
                        os.system(cmd)
                        #print(f)
                        break
                    else:
                        continue
                break
            else:
                continue
    

    print('Stack records:')
    for d in tqdm(glob(join(middle_for_stack,'201*'))):
        for dd in glob(join(d,'ZDY*')):
            staname = basename(dd).split('.')[0]
            chaname = basename(dd).split('.')[1]
            # for root,dirs,files in os.walk(d):
            #     for f in files:
            #         if os.path.splitext(f)[1] == '.SAC':
            #             count += 1
            stacknumber = 0
            p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
            s = ''
            s += 'sss\n'
            s += 'zerostack\n'
            stlo = 0
            stla = 0
            for f in glob(join(dd,'*.SAC')):
                cmd = 'saclst stlo stla f %s ' % (f)
                junk,stlo,stla = os.popen(cmd).read().split()
                s += 'addstack %s\n' % (f)
                stacknumber += 1
            name = ['DISPL',staname,chaname,'SAC']
            name = '.'.join(name)
            s += 'timewindow -20 180\n'
            s += 'sumstack\n'
            s += 'writestack %s\n' % (name)
            s += 'quitsub\n'
            s += 'r %s\n' % (name)
            if int(year) == 2015:
                s += 'ch nzyear 2015 nzjday 290 nzhour 00 nzmin 00 nzsec 00 nzmsec 000\n'
            elif int(year) == 2016:
                s += 'ch nzyear 2016 nzjday 291 nzhour 00 nzmin 00 nzsec 00 nzmsec 000\n'
            elif int(year) == 2017:
                s += 'ch nzyear 2017 nzjday 293 nzhour 00 nzmin 00 nzsec 00 nzmsec 000\n'
            s += 'ch evlo 100.098 evla 38.7485 stlo %.6f stla %.6f\n' % (float(stlo),float(stla))
            # s += 'ch o 20\n'
            s += 'wh\n'
            s += 'q\n'
            p.communicate(s.encode())
            cmd = 'rm -r %s' % (dd)
            os.system(cmd)

        print('Convert the file format from SAC into MiniSEED:')
        for f in tqdm(glob('DISPL.*')):
            traces = io.load(f, format='sac')
            out_filename = f[:-4] + '.mseed'
            io.save(traces, out_filename)
        
        print('Rename the MiniSEED files and Reattribution:')
        cmd = ''
        for f in tqdm(glob('*.mseed')):
            amplitude_type = basename(f).split('.')[0]
            staname = basename(f).split('.')[1]
            chaname = basename(f).split('.')[2]
            name = [amplitude_type,staname,chaname]
            name = '.'.join(name)
            os.rename(f,join(dire,name))
            
        cmd = 'mv %s %s\n' % ('*.SH?',d)
        cmd += 'rm DISPL.*\n'
        os.system(cmd)

##############################################################################

    # Convert the file format from SAC into MiniSEED, Rename & Reattribution

def FileFormatConversion_NoStack():

    print('Convert the file format from SAC into MiniSEED:')
    for f in tqdm(glob(join(middle_save,year+'*','No.*','*.SAC_RESP'))):
        cmd = "saclst nzyear nzjday nzhour nzmin nzsec nzmsec b e f %s" % (f)
        junk, kzyy, kzjd, kzhou, kzmin, kzsec, kzmsec, b, e = os.popen(cmd).read().split()
        kztime = UTCDateTime(year=int(kzyy),julday=int(kzjd),hour=int(kzhou),minute=int(kzmin),second=int(kzsec),microsecond=(int(kzmsec)*1000))
        btime = kztime + float(b)
        etime = kztime + float(e)
        byear = btime.year; bmon = btime.month; bday = btime.day; bhour = btime.hour; bmin = btime.minute; bsec = btime.second
        eyear = etime.year; emon = etime.month; eday = etime.day; ehour = etime.hour; emin = etime.minute; esec = etime.second
        traces = io.load(f, format='sac')
        out_filename = f + '.mseed'
        io.save(traces, out_filename)
        number = basename(f).split('.')[1]
        net = basename(f).split('.')[2]
        sta = basename(f).split('.')[3]
        cha = basename(f).split('.')[5]
        y = basename(f).split('.')[0][0:4]
        jday = basename(f).split('.')[0][4:7]
        name = 'waveform_%s_%s__%s_%d-%02d-%02d_%02d-%02d-%02d_%d-%02d-%02d_%02d-%02d-%02d_%s_%s.mseed' % (net,sta,cha,byear,bmon,bday,bhour,bmin,bsec,eyear,emon,eday,ehour,emin,esec,number,jday)
        cmd = 'mv %s %s' % (out_filename,name)
        os.system(cmd)

    print('Reattribute the MiniSEED files according to the excitation time:')
    for f in glob(join(dire,'*.mseed')):
        filename = basename(f)
        datetime = filename.split('_')[5].split('-')[0] + filename.split('_')[10].split('.')[0]
        sequence = filename.split('_')[9]
        exit_flag = False
        for d in glob(join(middle_save,year+'*')):
            date = basename(d)
            for dd in glob(join(d,'No.*')):
                seq = basename(dd).split('.')[1]
                if int(datetime) == int(date) and int(sequence) == int(seq):
                    cmd = 'mv %s %s' % (f,dd)
                    os.system(cmd)
                    exit_flag = True
                    break
            if exit_flag:
                break 


    # print('Rename the MiniSEED files for compatible Kiwi Tool use')
    # for d in tqdm(glob(join(middle_save,year+'*','No.*'))):
    #     for f in glob(join(d,'*.mseed')):
    #         filename = basename(f)
    #         staname = filename.split('.')[3]
    #         chaname = filename.split('.')[5]
    #         chaname = 'B' + chaname[1:3]
    #         name = ['DISPL',staname,chaname]
    #         name = '.'.join(name)
    #         path = os.path.realpath(d)
    #         os.rename(f,join(path,name))

    # for d in glob(join(middle_save,year+'*')):
    #     cmd = 'rm %s/*.SAC\n' % (d)
    #     cmd += 'rm %s/*.mseed\n' % (d)
    #     os.system(cmd)

##############################################################################

    # Convert the file format from SAC into MiniSEED, Reattribution & Rename

def FileFormatConversion_Stack():
    print('Stack the single day\' all shots into one record:')
    for d in glob(join(middle_save,year+'???')):
        y = basename(d)[0:4]
        jd = basename(d)[4:7]
        if int(y) == int(year) and int(jd) >= int(min_day) and int(jd) <= int(max_day):
            julday = basename(d)[4:7]
            stack_dir = join(d,'stack_dir') 
            cmd = 'mkdir %s' % (stack_dir)
            os.system(cmd)
            for i in range(1,41):
                sta = 'ZDY%02d' % (i)
                d1 = join(os.path.realpath(d),'%s.SHE') % (sta)
                d2 = join(os.path.realpath(d),'%s.SHN') % (sta)
                d3 = join(os.path.realpath(d),'%s.SHZ') % (sta)
                cmd = 'mkdir %s %s %s' % (d1,d2,d3)
                os.system(cmd)
        
            for dd in glob(join(d,'ZDY*')):
                sta = basename(dd).split('.')[0]
                cha = basename(dd).split('.')[1]
                for f in glob(join(d,'No.*','*.SAC_RESP')):
                    station = basename(f).split('.')[3]
                    channel = basename(f).split('.')[5]
                    stlo = 0
                    stla = 0
                    if station == sta and channel == cha:
                        cmd = 'cp %s %s' % (f,dd)
                        os.system(cmd)

            for dd in glob(join(d,'ZDY*')):
                sta = basename(dd).split('.')[0]
                cha = basename(dd).split('.')[1]
                stacknumber = 0
                p = subprocess.Popen(['sac'], stdin=subprocess.PIPE)
                s = ''
                s += 'sss\n'
                s += 'zerostack\n'
                for f in glob(join(dd,'*.SAC_RESP')):
                    cmd = 'saclst stlo stla f %s ' % (f)
                    junk,stlo,stla = os.popen(cmd).read().split()
                    s += 'addstack %s\n' % (f)
                    stacknumber += 1
                name = ['linestack',basename(d),sta,cha,str(stacknumber),'SAC_RESP']
                name = '.'.join(name)
                s += 'timewindow 0 100\n'
                s += 'sumstack\n'
                s += 'writestack %s\n' % (name)
                s += 'quitsub\n'
                s += 'r %s\n' % (name)
                s += 'ch nzyear %d nzjday %d nzhour 00 nzmin 00 nzsec 00 nzmsec 000\n' % (int(year),int(julday))
                s += 'ch evlo 100.098 evla 38.7485 stlo %.6f stla %.6f\n' % (float(stlo),float(stla))
                s += 'ch o 20\n'
                # s += 'bp c 1 9 n 4 p 2\n'
                s += 'wh\n'
                s += 'q\n'
                p.communicate(s.encode())
            cmd = 'mv linestack.* %s\n' % (stack_dir)
            os.system(cmd)

        # print('Convert the file format from SAC into MiniSEED:')
        # for f in tqdm(glob(join(stack_dir,'*.SAC_RESP'))):
        #     cmd = "saclst nzyear nzjday nzhour nzmin nzsec nzmsec b e f %s" % (f)
        #     junk, kzyy, kzjd, kzhou, kzmin, kzsec, kzmsec, b, e = os.popen(cmd).read().split()
        #     kztime = UTCDateTime(year=int(kzyy),julday=int(kzjd),hour=int(kzhou),minute=int(kzmin),second=int(kzsec),microsecond=(int(kzmsec)*1000))
        #     btime = kztime + float(b)
        #     etime = kztime + float(e)
        #     byear = btime.year; bmon = btime.month; bday = btime.day; bhour = btime.hour; bmin = btime.minute; bsec = btime.second
        #     eyear = etime.year; emon = etime.month; eday = etime.day; ehour = etime.hour; emin = etime.minute; esec = etime.second
        #     traces = io.load(f, format='sac')
        #     out_filename = f + '.mseed'
        #     io.save(traces, out_filename)
        #     total_stack = basename(f).split('.')[3]
        #     net = 'G3'
        #     sta = basename(f).split('.')[1]
        #     cha = basename(f).split('.')[2]
        #     y = basename(f).split('.')[0][0:4]
        #     jday = basename(f).split('.')[0][4:7]
        #     name = 'waveform_%s_%s__%s_%d-%02d-%02d_%02d-%02d-%02d_%d-%02d-%02d_%02d-%02d-%02d_%s.mseed' % (net,sta,cha,byear,bmon,bday,bhour,bmin,bsec,eyear,emon,eday,ehour,emin,esec,total_number)
        #     cmd = 'mv %s %s\n' % (out_filename,name)
        #     cmd += 'mv %s %s\n' % (name,stack_dir)
        #     os.system(cmd)

###############################################################################

    # Components Rotation

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
        input_dire1=glob(join(input_dire,"*.SAC"))
        input_dire1.sort()
        for freq in range(time):
            file_1 = input_dire1[freq*3]
            file_2 = input_dire1[freq*3+1]
            file_3 = input_dire1[freq*3+2]
            filename = os.path.basename(file_1)
            date = filename.split(".")[0]
            number = filename.split(".")[1]
            net = filename.split(".")[2]
            sta = filename.split(".")[3]
            loc = filename.split(".")[4]
            # number = filename.split(".")[1]
            # y = date[0:4]
            # jday = date[4:7]

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
            s += "w %s %s\n" % (join(date+'.'+number+'.'+net+'.'+sta+'.'+loc+'.R'), join(date+'.'+number+'.'+net+'.'+sta+'.'+loc+'.T'))
            s += "cut off\n"
            shutil.copy(file_3, join(date+'.'+number+'.'+net+'.'+sta+'.'+loc+'.Z'))
        s += "q\n"
        SAC.communicate(s.encode())
        os.system("mv *.[rtz] %s/%s_rot/" % (input_dire,basename(input_dire)))

###############################################################################
    
    ###########################

    # MAIN PROGRAM
    
    ###########################

ReattributeMiniSEED()
HeaderInfoWrite()
CutEvents()
PreProcess_and_Reattribution()
# Stack()
# ComponentRotation()
# FileFormatConversion_NoStack()
FileFormatConversion_Stack()

