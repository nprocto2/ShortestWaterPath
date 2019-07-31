import math
import global_land_mask as globe
import networkx as nx
import shapely.geometry as sgeom
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from geopy.distance import geodesic
import pickle


# function to return the path after being passed in the graph, starting point, ending point and type of algorithm
# theta = 0 is Lazy Theta
# theta = 1 is A*
# theta = 2 is A* with smoothing
def thetastarattempt4(G, start, end, theta, weight='weight'):    
    def h(p1, p2):
        d = geodesic(p1, p2).meters
        return d
        
    if theta == '1' or theta == '2':
        path = nx.astar_path(G, start, end, h, weight = 'weight')
        return path
    
    g = {}
    parent = {}
    opendict = {}
    next_to = {}
    temp = {}
    closeset = set()
    holder = set()
    
    g[start] = 0
    parent[start] = start
    opendict[start] = g[start] + h(start, end)
    while opendict is not None:
        s = min(opendict, key=opendict.get)
        del opendict[s]       
        
        for n, d in G[s].items():
            next_to.setdefault(s, set())
            next_to[s].add(n)
            
        if LOS(parent[s], s) == False:
            holder = next_to[s] & closeset
            hold = list(holder)
            for x in range(0, len(hold)):
                temp = dict.fromkeys(hold, g[hold[x]]+h(hold[x], s))
            t2 = min(temp, key=temp.get)
            parent[s] = t2
            g[s] = temp[t2]
            
        if s == end:
            total_path = [s]
            while s != parent[s]:
                total_path.append(parent[s])
                s = parent[s]
            total_path.reverse()
            return total_path
        
        closeset.add(s)
        
        for neighbor, w in G[s].items(): 
            if neighbor not in closeset:
                if neighbor not in opendict:
                    g[neighbor] = math.inf
                    parent[neighbor] = None
                
                if g[parent[s]] + h(parent[s], neighbor) < g[neighbor]:
                    g[neighbor] = g[parent[s]] + h(parent[s], neighbor)  
                    parent[neighbor] = parent[s]
                    
                    if neighbor in opendict:
                        del opendict[neighbor]
                    opendict[neighbor] = g[neighbor] + h(neighbor, end)
    return 0    

# function to return geodetic distance between two lat/lon points
def Geodesic(lat1, lon1, lat2, lon2):
   a = lat1, lon1
   b = lat2, lon2
   d = geodesic(a, b).meters
   return d

# function to determine if there is Line of Sight between two coordinants
# returns false if any coordinate between the two points is land and true if only water 
def LOS(a, b):
    lat1, lon1 = a
    lat2, lon2 = b
    fraction = 1/20
    temp = [IntermediatePointF(lat1, lon1, lat2, lon2, x*fraction) for x in range(1, int(1/fraction))]
    
    latplot, lonplot = list(zip(*temp))
    
    value = [globe.is_land(latplot[x], lonplot[x]) for x in range(len(temp))]
    
    if any(value) == True:
        return False
    
    else:
        return True

# function to interpolate between two lat/lon coordinants to allow more accurate line of sight checking
# returns a list of coordinants based on great circles data 
def IntermediatePointF(lat1, lon1, lat2, lon2, f):
    r = 6371008.0
    d = Geodesic(lat1, lon1, lat2, lon2)

    if d == 0:
        return [lat1, lon1]
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)

    a = math.sin((1-f)*d/r)/(math.sin(d/r))
    b = math.sin(f*d/r)/math.sin(d/r)
    
    x = a*math.cos(lat1)*math.cos(lon1) + b*math.cos(lat2)*math.cos(lon2)
    y = a*math.cos(lat1)*math.sin(lon1) + b*math.cos(lat2)*math.sin(lon2)
    
    z = a*math.sin(lat1) + b*math.sin(lat2)
    
    newlat = math.degrees(math.atan2(z, math.sqrt(x*x + y*y)))
    newlon = math.degrees(math.atan2(y, x))
    
    return [newlat, newlon]

# function to allow any coordinate to be chosen as a start or ending point as long as it is within 5 degrees of water
def fixcoords(start, watercoords):
    lat, lon = start
    a = lat - 5
    b = lat + 5
    c = lon - 5
    d = lon + 5
    l = []
    holder = [(x,y) for (x,y) in watercoords if a < x < b and c < y < d]
    l = [geodesic(start, x).meters for x in holder]
    return holder[l.index(min(l))]

# function to smooth A* by returning a list of unnessary points
def DeletePath(path, a):
    deletepath = []
    for z in range(a+2, len(path)):
        if LOS(path[a], path[z]):
           deletepath.append(path[z-1])
    return deletepath 

# plotting    
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_extent([-180, 180, -90, 90], ccrs.Geodetic())

# for theta: if '0' then lazy theta algorithm, if '1' then A*, and if '2' then A* with smoothing
start = (38.134752, 2.499829)
end = (27.958174, -75.800273)
theta = '0'

# loading files in
infile1 = open('Water Half Coords', 'rb')
watercoords = pickle.load(infile1)
infile1.close()
G = nx.read_gpickle("Half Degree Graph V4")

if start not in watercoords or end not in watercoords:
    start = fixcoords(start, watercoords)
    end = fixcoords(end, watercoords)

path = thetastarattempt4(G, (start), (end), theta, weight='weight')

# smoothing A* path
if theta == '2':
    for z in range(0, len(path)):
        Del = DeletePath(path, z)
        path = [d for d in path if d not in Del]


latplot, lonplot = list(zip(*path))
track = sgeom.LineString(zip(lonplot, latplot))
ax.add_geometries([track], ccrs.Geodetic(), facecolor='none', edgecolor='blue')

plt.scatter(start[1], start[0], color='lime', transform=ccrs.Geodetic())
plt.scatter(end[1], end[0], color='red',transform=ccrs.Geodetic())

plt.show()