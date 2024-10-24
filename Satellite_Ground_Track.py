#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 03:16:10 2023

@author: bouboujr
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import time as t
import calendar as calendar
from datetime import datetime

def Ground_trace(M0,T,e,w,i,Lw,Wt,V,interval,satellite_name,proj,t0,anim_globe,blue_marble):
    """
    Calculate and plot the ground trace of a satellite depending of the choosing projection.

    Parameters
    ----------
    M0 : float
        Initial mean anomaly of the satellite (degrees).
    T : float
        Orbital period (seconds).
    e : float
        Orbital eccentricity (how elongated the orbit is).
    w : float
        Argument of periapsis (degrees, orientation of the closest approach).
    i : float
        Inclination (degrees, the tilt of the orbit relative to Earth's equator).
    Lw : float
        Longitude of the ascending node (degrees).
    Wt : float
        Earth's rotational velocity (degrees/second).
    V : int
        True anomaly (degrees, position of the satellite in its orbit).
    interval : int
        Interval of true anomalies to compute.
    satellite_name : str
        The name of the satellite being tracked.
    projection : str
        Map projection type ('mercator' or 'globe').
    t0 : list
        Initial time in [year, month, day, hour, minute, second] format.
    anim_globe : bool
        If True, creates an animated globe projection of the satellite's motion.
    blue_marble : bool
        If True, uses the famous Blue Marble's picture on the projection
    
    Returns
    -------
    Table : np.array
        Table containing satellite position data (anomalies, latitude, longitude).
    """
    # Initialize an empty table to store satellite data
    Table=np.zeros((6,V//interval+1))
    
    j=1 

    for v in range(0,V+1,interval): # Loop over true anomaly values
        
        # Compute the latitude based on inclination and argument of periapsis
        latitude=np.degrees(np.arcsin(np.sin((w+v)*np.pi/180)*np.sin(i*np.pi/180)))
            
        
        # If the sine of the true anomaly is positive or zero
        if np.sin(v*np.pi/180)>=0:
            
            if v>360: # Adjust if the true anomaly exceeds 360 degrees
                
                E=np.arccos((e+np.cos(v*np.pi/180))/(1+e*np.cos(v*np.pi/180)))*180/np.pi + 360
                M=E*(np.pi/180)-e*np.sin(E*np.pi/180)
                M=M*180/np.pi
                dt=(M-M0)*T/360

            else: # Normal computation
                
                E=np.arccos((e+np.cos(v*np.pi/180))/(1+e*np.cos(v*np.pi/180)))*180/np.pi 
                M=E*(np.pi/180)-e*np.sin(E*np.pi/180)
                M=M*180/np.pi
                dt=(M-M0)*T/360
        
        # If the sine of the true anomaly is negative
        elif np.sin(v*np.pi/180)<0:
            
            E=np.arccos(-(e+np.cos((v)*np.pi/180))/(1+e*np.cos((v)*np.pi/180)))*180/np.pi +180
            M=E*(np.pi/180)-e*np.sin(E*np.pi/180)
            M=M*180/np.pi
            dt=(M-M0)*T/360
            
            j+=1
        
        # Compute the longitude based on the Earth's rotation and the current anomaly
        if np.cos((w+v)*np.pi/180)>0:
        
            longitude=np.degrees(np.arcsin(np.tan(latitude*np.pi/180)/np.tan(i*np.pi/180)))+Lw-Wt*dt
        
        elif np.cos((w+v)*np.pi/180)<0:
            
            longitude=180-np.degrees(np.arcsin(np.tan(latitude*np.pi/180)/np.tan(i*np.pi/180)))+Lw-Wt*dt
            
        # Fill in the table with the calculated data    
        Table[0, v//interval]=v
        Table[1, v//interval]=E
        Table[2, v//interval]=M
        Table[3, v//interval]=dt
        Table[4, v//interval]=latitude
        Table[5, v//interval]=longitude
    
    # Find the closest point to the initial anomaly (M0)
    j=0  
    inter=10000
    
    while j!=(np.shape(Table)[1]-1):
        
        if abs(Table[3,j])<inter:
            inter2=Table[3,j]
            inter=abs(Table[3,j])
            c=j
            
        j+=1
    
    # If projection is Mercator, plot the data on a Mercator map
    if proj=='mercator':
        
        if blue_marble:
            # Set up the map
            map = Basemap(projection='merc',llcrnrlat=-85,urcrnrlat=85,
                        llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution = None)
            map.bluemarble()
        
        else:
            # Set up the map
            map = Basemap(projection='merc',llcrnrlat=-85,urcrnrlat=85,
                        llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution = 'c')
            
            map.drawcoastlines()
            map.fillcontinents(color='green',lake_color='blue')
            
            map.drawparallels(np.arange(-90.,90.,10.),labels=[1,1,0,0])
            map.drawmeridians(np.arange(-180.,181.,30),labels=[0,0,0,1])
            map.drawmapboundary(fill_color='blue')
        
        
        # Plot the satellite's trajectory
        map.scatter(Table[5,:],Table[4,:],latlon=True,marker='.',color='r',linewidth=1)
        map.scatter(Table[5,c],Table[4,c],latlon=True,marker='.',color='w',linewidth=2)
        map.scatter(Table[5,0],Table[4,0],latlon=True,marker='.',color='w',linewidth=2)
        
        # Plot annotations for key points
        bbox_props = dict(fc='white', ec='green', lw=2, boxstyle='round,pad=0.1')
        font = {'fontname': 'Arial'}
        date = datetime.fromtimestamp(t.mktime(t.gmtime(calendar.timegm(t.struct_time(t0)))))
    
        CS=map.nightshade(date,alpha=0.4,delta=0.25)
                
        
        plt.annotate('{}\n(UTC) noté (t0 + {} s)'.format(t.asctime(t.gmtime(calendar.timegm(t.struct_time(t0))+inter2)), inter2), 
                    map(Table[5,c],Table[4,c]), xytext=(-0.57,-0.125), textcoords='axes fraction',
                    fontsize=20, weight=1000, **font,
                    arrowprops=dict(arrowstyle='-', connectionstyle='angle3'),
                    bbox = bbox_props)
        
        plt.annotate('{}\n(UTC)'.format(t.asctime(t.gmtime(calendar.timegm(t.struct_time(t0))+int(Table[3,1])))), 
                    map(Table[5,0],Table[4,0]), xytext=(-0.57,0.85), textcoords='axes fraction', 
                    fontsize=20, weight=1000, **font,
                    arrowprops=dict(arrowstyle='-', connectionstyle='angle3'),
                    bbox = bbox_props)
        
        plt.title("Ground trace of satellite {} with {} projection".format(satellite_name,proj))
        plt.show()
    
    
    elif proj=='globe':
        for i in range(np.size(Table[1,:])):

            if anim_globe:
                if blue_marble:
                    # Set up the map
                    map = Basemap(projection='ortho',lat_0=Table[4,i],lon_0=Table[5,i],resolution=None)
                    map.bluemarble()
                
                else:
                    # Set up the map
                    map = Basemap(projection='ortho',lat_0=Table[4,i],lon_0=Table[5,i],resolution='c')
                    map.drawcoastlines()
                    map.fillcontinents(color='green',lake_color='blue')
                    
                    map.drawparallels(np.arange(-90.,90.,10.),labels=[1,1,0,0])
                    map.drawmeridians(np.arange(-180.,181.,30),labels=[0,0,0,1])
                    map.drawmapboundary(fill_color='blue')
                    
            else:
                
                
                if blue_marble:
                    # Set up the map
                    map = Basemap(projection='ortho',lat_0=Table[4,c],lon_0=Table[5,c],resolution=None)
                    map.bluemarble()
                
                else:
                    # Set up the map
                    map = Basemap(projection='ortho',lat_0=Table[4,i],lon_0=Table[5,i],resolution='c')
                    
                    map.drawcoastlines(linewidth=0.25)
                    map.drawcountries(linewidth=0.25)
                    map.fillcontinents(color='green',lake_color='blue')
                    
                    map.drawmapboundary(fill_color='blue')
                    
                    map.drawmeridians(np.arange(0,360,30))
                    map.drawparallels(np.arange(-90,90,30))
            
            date = datetime.fromtimestamp(t.mktime(t.gmtime(calendar.timegm(t.struct_time(t0)))))
            
            CS=map.nightshade(date,alpha=0.4,delta=0.25)
            
            map.scatter(Table[5,:],Table[4,:],latlon=True,marker='.',color='r',linewidth=1)
            map.scatter(Table[5,c],Table[4,c],latlon=True,marker='.',color='w',linewidth=1)
            map.scatter(Table[5,-1],Table[4,-1],latlon=True,marker='.',color='w',linewidth=1)
            map.scatter(Table[5,0],Table[4,0],latlon=True,marker='.',color='w',linewidth=1)
            
            map.scatter(Table[5,i],Table[4,i],latlon=True,marker='s',color='c',linewidth=1)
            
            plt.title("Ground trace of satellite {} with {} projection".format(satellite_name,proj))
            plt.pause(0.1)
            plt.cla()
            
    Table=np.hstack(([['v'],['E'],['M'],['dt'],['latitude'],['longitude']],Table))
    
    return Table

def Ground_trace_anim_merca(M0, T, e, w, i, Lw, Wt, V, interval, satellite_name,t0):
    """
    Calculate and plot the ground trace of a satellite over time using mercator's projection.

    Parameters
    ----------
    M0 : float
        Initial mean anomaly of the satellite (degrees).
    T : float
        Orbital period (seconds).
    e : float
        Orbital eccentricity (how elongated the orbit is).
    w : float
        Argument of periapsis (degrees, orientation of the closest approach).
    i : float
        Inclination (degrees, the tilt of the orbit relative to Earth's equator).
    Lw : float
        Longitude of the ascending node (degrees).
    Wt : float
        Earth's rotational velocity (degrees/second).
    V : int
        True anomaly (degrees, position of the satellite in its orbit).
    interval : int
        Interval of true anomalies to compute.
    satellite_name : str
        The name of the satellite being tracked.
    projection : str
        Map projection type ('mercator' or 'globe').
    t0 : list
        Initial time in [year, month, day, hour, minute, second] format.

    Returns
    -------
    Table : np.array
        Table containing satellite position data (anomalies, latitude, longitude).

    """
    Table=np.zeros((6,V//interval+1))
    
    j=1

    for v in range(0,V+1,interval):
        
        latitude=np.degrees(np.arcsin(np.sin((w+v)*np.pi/180)*np.sin(i*np.pi/180)))
            
        if np.sin(v*np.pi/180)>=0:
            
            if v>360:
                
                E=np.arccos((e+np.cos(v*np.pi/180))/(1+e*np.cos(v*np.pi/180)))*180/np.pi + 360
                M=E-e*np.sin(E)
                dt=(M-M0)*T/360

            else:
                
                E=np.arccos((e+np.cos(v*np.pi/180))/(1+e*np.cos(v*np.pi/180)))*180/np.pi 
                M=E-e*np.sin(E)
                dt=(M-M0)*T/360
        
        elif np.sin(v*np.pi/180)<0:
            
            E=np.arccos((e+np.cos((v)*np.pi/180))/(1+e*np.cos((v)*np.pi/180)))*180/np.pi +(2*interval)*j-interval
            M=E-e*np.sin(E)
            dt=(M-M0)*T/360
            
            j+=1
        
        if np.cos((w+v)*np.pi/180)>0:
        
            longitude=np.degrees(np.arcsin(np.tan(latitude*np.pi/180)/np.tan(i*np.pi/180)))+Lw-Wt*dt-180
        
        elif np.cos((w+v)*np.pi/180)<0:
            
            longitude=np.degrees(np.pi-np.arcsin(np.tan(latitude*np.pi/180)/np.tan(i*np.pi/180)))+Lw-Wt*dt-180
            
            
        Table[0, v//interval]=v
        Table[1, v//interval]=E
        Table[2, v//interval]=M
        Table[3, v//interval]=dt
        Table[4, v//interval]=latitude
        Table[5, v//interval]=longitude
    
    j=0  
    inter=50
    c=0
    
    while j!=(np.shape(Table)[1]-1):
        
        if abs(Table[3,j])<inter:
            inter=abs(Table[3,j])
            c=j
            
        j+=1
        
    map = Basemap(projection='merc',llcrnrlat=-85,urcrnrlat=85,\
                llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
        
    map.drawcoastlines()
    map.fillcontinents(color='green',lake_color='blue')
    # draw parallels and meridians.
    map.drawparallels(np.arange(-90.,90.,10.),labels=[1,1,0,0])
    map.drawmeridians(np.arange(-180.,181.,30),labels=[0,0,0,1])
    map.drawmapboundary(fill_color='blue')
    
    date = datetime.fromtimestamp(t.mktime(t.gmtime(calendar.timegm(t.struct_time(t0)))))
    
    CS=map.nightshade(date,alpha=0.4,delta=0.25)
    
    bbox_props = dict(fc='white', ec='green', lw=2, boxstyle='round,pad=0.1')
    
    font = {'fontname': 'Arial'}
    
    for k in range(np.size(Table[1,:])):
        
        
        if k==0 or k==c or k==(np.size(Table[1,:])-1):
            
            map.scatter(Table[5,k],Table[4,k],latlon=True,marker='.',color='w',linewidth=2)
            
            if k==0:
                plt.annotate('{}\n(UTC)'.format(t.asctime(t.gmtime(calendar.timegm(t.struct_time(t0))+int(Table[3,1])))), 
                            map(Table[5,k],Table[4,k]), xytext=(-0.57,0.85), textcoords='axes fraction', 
                            fontsize=20, weight=1000, **font,
                            arrowprops=dict(arrowstyle='-', connectionstyle='angle3'),
                            bbox = bbox_props)
            elif k==c:
                plt.annotate('{}\n(UTC) noté (t0)'.format(t.asctime(t.gmtime(calendar.timegm(t.struct_time(t0))))), 
                            map(Table[5,k],Table[4,k]), xytext=(-0.57,-0.005), textcoords='axes fraction',
                            fontsize=20, weight=1000, **font,
                            arrowprops=dict(arrowstyle='-', connectionstyle='angle3'),
                            bbox = bbox_props)
            
            elif k==np.size(Table[1,:])-1:
                plt.annotate('{}\n(UTC)'.format(t.asctime(t.gmtime(calendar.timegm(t.struct_time(t0))+int(Table[3,-1])))), 
                            map(Table[5,-1],Table[4,-1]), xytext=(1.07,0.85), textcoords='axes fraction', 
                            fontsize=20, weight=1000, **font,
                            arrowprops=dict(arrowstyle='-', connectionstyle='angle3'),
                            bbox = bbox_props)
        
        else:
            map.scatter(Table[5,k],Table[4,k],latlon=True,marker='.',color='y',linewidth=1)
        
        plt.title("Ground trace of satellite {} with {} projection".format(satellite_name,proj))
        
        plt.pause(0.0001)
        
    Table = np.hstack(([['v'],['E'],['M'],['dt'],['latitude'],['longitude']],Table))
    
    return Table
    
#%% PNEO 3 Mercator Projection

M0 = 284.5280                # deg   Anomalie moyenne
T = 97.06 * 60               # s     Période orbitale
e = 0.0001305                #       Excentricité
w=75.6067                    # deg   Argument du périgée
i=97.8934                    # deg   Inclinaison
Lw=85.5377                   # deg   Longitude du noeud ascendant à t0
Wt=360/86164                 # deg/s Vitesse angulaire
V=480                        # deg   Anomalie vrai
interval=2                   #       Intervalle des anomalies vrais
satellite_name='PNEO3'       # str   Nom du satellite étudié
proj='mercator'              # str   Type de projection ("mercator" ou "globe")
t0=[2023,1,8,18,10,8,6,8,-1] # list  [Année,mois,jour,heure,minute,secondes,n-ième jour de la semaine, n-ième jour de l'année, laisser sur -1]
anim_globe = True            # bool  Suivi du satellite en projection "globe"
blue_marble = True           # bool  

#Table = Ground_trace(M0, T, e, w, i, Lw, Wt, V, interval, satellite_name,proj,t0,anim_globe,blue_marble)    

Table = Ground_trace_anim_merca(M0, T, e, w, i, Lw, Wt, V, interval, satellite_name, t0)
