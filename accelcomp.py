#!/usr/bin/env python

# ###################################
# Code: accelcomp.py
# Adam Ringler
#
# This code compares broadband data to accelerometers by deconvolving all of the data to displacement
# ###################################
# ##################################
# Here are the different methods
# getorientation()
# getdip()
# rotatehorizontal()
# choptocommon()
# getlatlon()
# getstalist()
# readcmt()
# getdata()
# getPAZ2()
# getcolor()
# writestats()
# ##################################

import sys
import os
import glob
import numpy
import matplotlib
import math
import warnings
import argparse
from matplotlib.pyplot import (figure,axes,plot,xlabel,ylabel,title,subplot,legend,savefig,show,xlim)
from obspy import read, Stream, read_events
from obspy.core import UTCDateTime
from obspy.io.xseed import Parser
from time import gmtime, strftime
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.cross_correlation import xcorr
from obspy.taup import TauPyModel
from obspy.clients.fdsn import Client
from shutil import copyfile

# Here are the various parameters a user might change
# We might want to put these in a config file to avoid them sitting in the code

datalessloc = '/APPS/metadata/SEED/'
# Here is the data location use True for xs0 otherwise use false






def getargs():

    parser = argparse.ArgumentParser(description = "Program to compare Accelerometer and Broadband data")

    parser.add_argument('-n', type = str, action = "store", \
        dest = "network", required = True, help = "Network name Example: IU")

    parser.add_argument('-cmt', type = str, action = "store", \
        dest = "cmt", required = True, help = "CMT Solution")

    parser.add_argument('-resDir',type = str, action = "store", \
        dest = "resDir", required = True, help = "Result directory name Example: blah")

    parser.add_argument('-debug', action = "store_true", dest = "debug", \
        default = False, help = "Run in debug mode")

    parser.add_argument('-sta', type = str, action = "store", \
        dest = "sta", required = False, help = "Stations to use Example with a comma (,) separator : TUC,ANMO")

    parser.add_argument('-tslen', type = int, action ="store", \
        dest = "lents", required = False, help = "Length of time series in seconds Example:  2000, default is 2000 s")

    parser.add_argument('-dataloc', action = "store_true", dest = "dataloc", \
        default = False, help = "Use /msd data location, otherwise use /tr1 also")

    parser.add_argument('-filter', action = "store", nargs = 3, dest = "filter", required = False, \
        help = "Filter parameters using minimum period maximum period and number of corners Example: 100 200 4, " + \
            "default is 4 25 4")

    parser.add_argument('-PS', action = "store", nargs = 2, dest = "ps", required = False, \
        help = "Time before P and after S, default is 120 before P and 600 after S")

    parserval = parser.parse_args()
    return parserval








def choptocommon(stream):
#A function to chop the data to a common time window
    debugchoptocommon = False
    stimes = []
    etimes = []

#Lets get the start and end time for each trace in the stream
    for trace in stream:
        stimes.append(trace.stats.starttime)
        etimes.append(trace.stats.endtime)
    newstime = stimes[0]
    newetime = etimes[0]

#Lets now find the latest start time
    for curstime in stimes:
        if debug:
            print(curstime)
        if curstime >= newstime:
            newstime = curstime

#Lets find the earliest end time    
    for curetime in etimes:
        if debug:        
            print(curetime)
        if curetime <= newetime:
            newetime = curetime

    if debugchoptocommon:
        print(newstime)
        print(newetime)
        print(stream)

#Now we trim each trace by our latest start time and earliest end time
    for trace in stream:    
        trace.trim(starttime=newstime,endtime=newetime)
    if debug:
        print(stream)
    return stream





def readcmt(cmt, debug=False):
    """ read the CMT information contained in event object """
    hdur = cat[0].focal_mechanisms[0].moment_tensor.source_time_function.duration 
    cmtlat = cat[0].origins[0].latitude
    cmtlon = cat[0].origins[0].longitude
    eventtime = cat[0].origins[0].time
    dep = cat[0].origins[0].depth
    # The time shift is included in the event time for the cmt
    # This should eventually be removed
    tshift = cat[0].origins[0].time - cat[0].origins[1].time
    return cmtlat, cmtlon, eventtime, tshift, hdur, dep 

def getdata(net, sta, eventtime, lents, dataloc, debug=False):
#This function goes to one of the archives and gets the data

#Here are some hard coded values probably only necessary when messing with this function
    preeventday = eventtime - 24*60*60
    posteventday = eventtime + 24*60*60
    prepostwin= 3000
#If II get off of /tr1 else get the data from /xs0 or /xs1
    if net == 'II':
        dataprefix = '/tr1/telemetry_days/'
    else:
#Not using /tr1 data so lets get it from the archive

        dataprefix = '/msd/'
    if dataloc:
        dataprefix = '/tr1/telemetry_days/'   
    if debug:
        print('Here is the dataprefix:' + dataprefix)
    datastime = eventtime-prepostwin
    dataetime = eventtime+lents+prepostwin
    fstring = dataprefix + net + '_' + sta + '/'
    path = fstring + str(datastime.year) + '/*' + str(datastime.julday).zfill(3) +'*/*'
    st = read(path + 'LH*.seed')
    st += read(path + 'LN*.seed')
    if (dataetime.year > datastime.year) or (dataetime.julday > datastime.julday):
        path = fstring + str(dataetime.year) + '/*' + str(dataetime.julday).zfill(3) +'*/*'
        st +=  read(path + 'LH*.seed')
        st += read(path + 'LN*.seed')
    st.trim(starttime=datastime, endtime=dataetime)

    st.merge(fill_value='latest')
    if debug:
        print('We have data')
    return st



def getcolor(chan,loc):
#This function sets the color of the trace by channel and location

#If we have an accelerometer we use black
    if chan in set(['LXN','LXE','LXZ']):
        color = 'k'

#Primary sensors will be green
    elif (loc == '00' or loc ==''):
        color = 'g'

#Secondary sensors will be red
    elif loc == '10':
        color = 'r'

#Third sensors this could be an odd broadband will be cyan
    elif loc == '20':
        color = 'c'
    else:

#Everything else is black
        color = 'b'
    return color



def writestats(statfile,streamin,comp):
#This function does the final stat computations for the accelerometer plots just produced
#This was taken from the synthetics code so the accel plays the role of the synthetic
    debugwstats = False
    try:
        syncomp = "LN" + comp    
        datacomp = "LH" + comp
        
        if debugwstats:
            print(streamin)
            print('Here is the comp:' + syncomp)
        syn = streamin.select(channel = syncomp)
        if debugwstats:
            print(syn)
        for tr in streamin.select(channel = datacomp):    
            if debugwstats:
                print('Here is the trace value:' + str(numpy.sum(tr.data*syn[0].data)))
                print('Here is the accel value:' + str(numpy.sum(numpy.square(syn[0].data))))
#Here we compute a residual scale factor
            resi = "{0:.2f}".format(numpy.sum(tr.data*syn[0].data)/numpy.sum(numpy.square(syn[0].data)))
            if debugwstats:
                print('Here is the resi:' + str(resi))
#Lets compute a cross-correlation and lag
            lag, corr = xcorr(tr,syn[0],50)
            corr = "{0:.2f}".format(corr)
            if debugwstats:
                print('Here is the corr:' + str(corr))
                print('Here are the results:' + tr.stats.network + "," + tr.stats.station)
                print(' Here are more:' + "," + tr.stats.location + "," + tr.stats.channel + "," +  str(resi))
                print('And more:' + "," + str(lag) + "," + str(corr) + "\n")

#Now we want to write to a file
            statfile.write(tr.stats.network + "," + tr.stats.station)
            statfile.write("," + tr.stats.location + "," + tr.stats.channel + "," +  str(resi))
            statfile.write("," + str(lag) + "," + str(corr) + "\n")
    
    except:    
        if debug:
            print('No residual for' + cursta + ' ' + 'LH' + comp) 
    return



#Start of the main part of the program

parserval = getargs()

# Eventually allow this to change in the arguments
client = Client("IRIS")

if parserval.dataloc:
    dataloc = True
else:
    dataloc = False


debug = parserval.debug


if parserval.lents:
    lents = parserval.lents
else:
    lents = 2000

if parserval.filter:
    userminfre = 1.0/float(parserval.filter[1])
    usermaxfre = 1.0/float(parserval.filter[0])
    filtercornerpoles = int(parserval.filter[2])
else:
    userminfre = .05
    usermaxfre = .25
    #Use half the value you think you want e.g. 2 gives you a total of 4 poles
    filtercornerpoles = 4

if parserval.ps:
    bfarrival = int(parserval.ps[0])
    afarrival = int(parserval.ps[1])
else:
    #Here is the number of seconds before the P-wave arrival
    bfarrival = 120
    #Here is the number of seconds after the S-wave arrival
    afarrival = 600



# Read in the CMT solution from the synthetic directory
if debug:
    print("We are using local synthetics")
if not os.path.isfile(parserval.cmt):
    print("No CMT found")
    exit(0)
try:
    cat = read_events(parserval.cmt)
except:
    try:

        # Here we have a use case where we are missing a space
        copyfile(parserval.cmt, 'CMTTEMP')
        f=open('CMTTEMP','r')
        CMT = f.read()
        f.close()
        f = open('CMTTEMP', 'w')
        f.write(' ' + CMT)
        f.close()
        cat = read_events('CMTTEMP')
        os.remove('CMTTEMP')
    except:
        print("No CMT found")
        exit(0)
if debug:
    print(cat)        

cmtlat, cmtlon, eventtime, tshift, hdur, dep = readcmt(cat)
if eventtime.year <= 22:
    eventtime.year += 2000
elif eventtime.year <= 70:
    eventtime.year += 1900
#Lets make a local results directory
resultdir = parserval.resDir
if resultdir[-1] == '/':
    resultdir = resultdir[:-1]
if not os.path.exists(os.getcwd() + '/' + resultdir):
    os.mkdir(os.getcwd() + '/' + resultdir)

#Lets get the current network
curnet = parserval.network

#Lets open the results file to write
statfile = open(os.getcwd() + '/' + resultdir + '/Results' + curnet + '.csv' ,'w')
statfile.write('net,sta,loc,chan,scalefac,lag,corr\n')

model=TauPyModel(model="ak135")


#If we arent doing manual station lists we need to get one for the network
if parserval.sta:
    manstalist = True
    if debug: 
        print("We are using a manual station list")
    stalist = parserval.sta.split(",")
    stations = []
    for sta in stalist:
        stations.append(sta)
    if debug:
        print(stations) 
else:
    manstalist = False
    stations = client.get_stations(network=parserval.network, starttime=eventtime, endtime=eventtime)
    stations = [ sta.code for sta in stations[0]]

if debug:
    print("Here are the stations we found") 
    for sta in stations:
        print("Here is a station:" + sta)


#Lets start by using a station list and then move to a different approach
for sta in stations:
    
    net = parserval.network
    cursta = sta

#Now we get the data for the event
    try:
        st = getdata(net,cursta,eventtime,lents,dataloc)
    except:
        print('No data for ' + net + ' ' + cursta)
        continue
        
        
    # We should get some metadata
    try:
        inv = client.get_stations(starttime=eventtime, endtime=eventtime+5., network=net,
            sta = sta, channel="L*", level="response")
    except:
        continue
    #Lets go through each trace in the stream and deconvolve and filter
    for tr in st:
        #Here we get the response and remove it
        
        try:

            tr.taper(max_percentage=0.05, type='cosine')
            #If we have an accelerometer we want an extra zero to go to displacement
            tr.remove_response(inventory = inv, output='DISP')
            #Here we filter, integrate, taper, trim, detrend, and filter
            tr.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=filtercornerpoles)
            tr.taper(max_percentage=0.05, type='cosine')
            tr.trim(starttime=eventtime + tshift/2,endtime=(eventtime+lents + tshift/2))
            tr.detrend()
            tr.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=filtercornerpoles)
        except:
            print('Can not find the response')
            st.remove(tr)


    try:
        st.rotate('->ZNE', inventory=inv)
        st=choptocommon(st)
    except:
        continue

    #We want to get the distance of the event and of the station
    #We also want the back-azimuth
    lat = inv[0][0].latitude
    lon = inv[0][0].longitude

    dist= gps2dist_azimuth(float(cmtlat),float(cmtlon),lat,lon)
    bazi ="{0:.1f}".format(dist[2])
    distDeg = dist[0]
    dist ="{0:.1f}".format( 0.0089932 * dist[0] / 1000)
    if debug:
        print('Here is the distance:' + str(dist))
        print('Here is the depth:' + str(dep))

#Here is the travel time so we can do the final trim
#Should this be in a function to avoid it being in the main loop?
    tt = model.get_travel_times(distance_in_degree =float(distDeg), source_depth_in_km = dep/1000.) 
    firstarrival = tt[0].time
    for ttphase in tt:
        phasename = ttphase.name
        phasename = phasename[:1]
        if phasename == 'S':
            secondarrival = ttphase.time
            break
    newstime = st[0].stats.starttime + firstarrival - bfarrival
    newetime = st[0].stats.starttime + secondarrival + afarrival
    st.trim(starttime=newstime,endtime=newetime)


    if debug:
        print('Here is the first arrival time: ' + str(firstarrival))
        print('Here is the second arrival time: ' + str(secondarrival))

        
    # Here is the mess of plotting info  
    # Set the time series
    tz=numpy.arange(0,st[0].stats.npts / st[0].stats.sampling_rate, st[0].stats.delta)

    # Get a legend and plot the vertical
    synplot = figure(1)

    # Here we setup subplot 1 and do a title
    subplot(311)
    titlelegend = st[0].stats.network + ' ' + st[0].stats.station + ' '

    # Here is the start time of the plot
    stime = str(st[0].stats.starttime.year) + ' ' + str(st[0].stats.starttime.julday) + ' ' + \
    str(st[0].stats.starttime.hour) + ':' + str(st[0].stats.starttime.minute) + \
    ':' + str("{0:.2f}".format(st[0].stats.starttime.second))

    titlelegend = titlelegend + stime + ' ' 
    
    if debug:
        print("Latitude:" + str(lat))
        print("Longitude:" + str(lon)) 
        print("CMT Latitude:" + str(cmtlat))
        print("CMT Longitude:" + str(cmtlon))
    
    #Here we add the distance and the back-azimuth to the legend    
    titlelegend = titlelegend + 'Dist:' + str(dist) 
    titlelegend = titlelegend + ' BAzi:' + str(bazi) 

    #Here we add the frequencies
    minper = "{0:.0f}".format(1/usermaxfre)
    maxper = "{0:.0f}".format(1/userminfre)
    titlelegend = titlelegend + ' ' + str(minper) + '-' + str(maxper) + ' s per.'
    title(titlelegend,fontsize=12)
    st.sort(['location'])
    for comps in st.select(component="Z"):
        curcolor = getcolor(comps.stats.channel,comps.stats.location)
        plot(tz,(comps.data*(10**3)), curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
    legend(prop={'size':6})
    xlim(0,len(tz))
    st = choptocommon(st)
    st.sort(['location', 'channel'])
    if debug:
        print("Here is the final stream:")
        print(st)

    # We now plot the N/S component of the data
    subplot(312)
    tne=numpy.arange(0,st[0].stats.npts / st[0].stats.sampling_rate, st[0].stats.delta)
    for comps in st.select(component="N"):
        curcolor = getcolor(comps.stats.channel,comps.stats.location)
        plot(tne, (comps.data*(10**3)), curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
    legend(prop={'size':6})
    ylabel('Displacement (mm)')    
    xlim(0, len(tne))
    
    # Now we plot the E/W compoent of the data
    subplot(313)
    for comps in st.select(component = "E"):
        curcolor = getcolor(comps.stats.channel, comps.stats.location)
        plot(tne, (comps.data*(10**3)),curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
    legend(prop={'size':6})
    xlabel('Time (s)')
    xlim(0, len(tne))

    # Finally we need to save the figure
    savefig(os.getcwd() + '/' + resultdir + '/' + st[0].stats.network + cursta + \
    str(st[0].stats.starttime.year) + str(st[0].stats.starttime.julday) + \
    str(st[0].stats.starttime.hour) + str(st[0].stats.starttime.minute) + '.png', format = 'png', dpi=400)

    # Lets clear the plot so we have no residual
    synplot.clear()

    # Time to write some info into the statfile
    # Write the network and the station
    writestats(statfile, st, 'Z')
    writestats(statfile, st, 'N')
    writestats(statfile, st, 'E')
    
    
# Lets get an RMS from the synthetic and the data

statfile.close()
