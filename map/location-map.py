# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 13:20:01 2018

@author: samantha
"""

import csv
import matplotlib.pyplot as plt
from time import gmtime, strftime

def plot_map(lon, lat):
    m = plt.imread("ukphsoil.png")
    plt.imshow(m, extent=[-9,2,50,59])
    plt.axis([-9,2,50,59])
    plt.plot(lon, lat, marker='x', markersize='10', color='white')
    plt.axvline(x=lon, linewidth=2, color="white")
    plt.axhline(y=lat, linewidth=2, color="white")
    plt.savefig('location_soil-ph_%s.png' % (time_stamp))

def look_in_csv(city, csv_path):
    lon_lat = []
    with open(csv_path) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        try:
            for row in reader:
                if row[0] == city:
                    lon_lat = [float(row[2]), float(row[1])]
        except IndexError:
            print("Location not in database...")
    return lon_lat
    
##############################

if __name__ == "__main__":   
    time_stamp = strftime("%Y-%m-%d_%H-%M-%S", gmtime())
    print("Please note: due to the fact that the Earth is spherical and maps are 2-dimensional, there will be some distortion when plotting locations.")
    
    path = input("Insert location of csv filename: ") # /home/samantha/Documents/practice/latlon.csv
    city = input("Insert city (case-sensitive): ")
    lon_lat = look_in_csv(city, path)
    #plot_map(-4.064598, 52.41616)
    try:
        plot_map(lon_lat[0], lon_lat[1])
    except IndexError:
        print("Location not in database or spelt incorrectly (remember input is case sensitive)")