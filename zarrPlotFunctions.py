import numpy as np
import xarray as xr

metadata = 'PATH_TO_H5_FILE'
x = xr.open_dataset(metadata)
lat = x.latitude.data
lon = x.longitude.data

# Input Lat/Lon you are interested in to find chunk

def findChunk(latPoint,lonPoint):
    absLat = np.abs(lat - latPoint)
    absLon = np.abs(lon - lonPoint)
    c = np.maximum(absLon,absLat)
    index = np.where(c == c.min())

    x = int(index[0])
    y = int(index[1])
    
    # Divide by chunk size - ignoring remainder
    chunkX = (x//150)
    chunkY = (y//150)
    
    return chunkX, chunkY


def getLatLons(chunkX,chunkY):
    x0 = chunkX*150
    y0 = chunkY*150
    xRange = x0+150
    yRange = y0+150
    
    chunkLats = lat[x0:xRange,y0:yRange]
    chunkLons = lon[x0:xRange,y0:yRange]
    
    return chunkLats, chunkLons

def nearestGridPoint(ptLat,ptLon,latGrid,lonGrid):
    abslat = np.abs(ptLat - latGrid)
    abslon = np.abs(ptLon - lonGrid)
    c = np.maximum(abslon,abslat)
    ind = np.where(c == c.min())

    i = int(ind[0])
    j = int(ind[1])
    
    return i, j


def perdelta(start, end, delta):
    current = start
    while current < end:
        yield current
        current +=delta
