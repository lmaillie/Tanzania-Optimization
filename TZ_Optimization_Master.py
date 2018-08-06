############################################################
##                        PART 1
############################################################
## Write all road edges to graph H with weight = distance
############################################################

## Write all road edges to graph H with weight as distance

import networkx as nx
import matplotlib.pyplot as plt

# Start with using nx.read_shp to generate G
# Set simplify = False so that you get every road node, not just road segment ends
G=nx.read_shp('/Users/Luke/Documents/Tanzania_GIS/2018-07-02_ReUTM/roads_utm.shp', simplify=False) #CHANGE FILE PATH

from shapely.geometry import LineString #need for line.length

H = nx.Graph() #H is graph where new edges will be written to

# Convert G into edge list with weighted edges
for u,v in G.edges():
    n1 = u
    n2 = v
    
    #Need to have nodes as list
    line_seg = [] 
    line_seg.append(n1)
    line_seg.append(n2)
    
    #Length calculations
    line = LineString(line_seg)
    L = line.length
    
    #Create new weighted graph
    H.add_edge(n1,n2,distance=L)

print('First 3 Edges',list(H.edges(data=True))[:3]) #See data format

print('Distance of first road seg',list(H.edges(data=True))[0][2]['distance']) #Accessing distance

print('Number of road segments',len(list(H.edges()))) #Should be 700881 so 55 short...

## ADD ATTRIBUTE TO NODES TO IDENTIFY AS ROAD
## Note identifier key: 0 = road, 1 = ward center, 2 = health facility

nx.set_node_attributes(H, 0, 'identifier') #Want to know what kind of node it is (road)

print('Identifer for first node',list(H.nodes(data=True))[0][1]['identifier']) #Accessing identifier


## PAIR EACH WARD CENTER WITH ITS POPULATION

from osgeo import ogr
import json

ALL_coord_pair_wc = []
ALL_pop_wc = []

l = range(3644) #Total of 3644 ward centers (wc's)
for a in l:
    file = ogr.Open("/Users/Luke/Documents/Tanzania_GIS/2018-07-02_ReUTM/TZ_CentroidsWith2012Pop_UTM.shp") #CHANGE FILE PATH
    shape = file.GetLayer(0)
    feature = shape.GetFeature(a)
    
    if feature == None:
        coord_pair=[0,0]
        ward_pop = 0

    else:
        coord_pair = json.loads(feature.ExportToJson())['geometry']['coordinates']
        ward_pop = json.loads(feature.ExportToJson())['properties']['total_both']
    
    a = tuple(coord_pair) #Converts coordinate pair from list to tuple
    b = ward_pop
    
    ALL_coord_pair_wc.append(a)
    ALL_pop_wc.append(b)

wc_pop_dict = list(zip(ALL_coord_pair_wc,ALL_pop_wc)) # Completes pairing ward center with pop

print('1st WC Population,'(wc_pop_dict)[0][1]) #Gives you population for 1st ward center

############################################################
##                        PART 2
############################################################
##               Add all WC's to graph H
############################################################

## Store list of linestrings of road segments to access faster later

from osgeo import ogr
import json

l = range(34144) #Total of 34144 road segments
LIST_ROAD_SEG = []

for a in l:
    file = ogr.Open("/Users/Luke/Documents/Tanzania_GIS/2018-07-02_ReUTM/roads_utm.shp") #CHANGE FILE PATH
    shape = file.GetLayer(0)
    feature = shape.GetFeature(a)
    first = json.loads(feature.ExportToJson())['geometry']['coordinates']
    LIST_ROAD_SEG.append(first)

print('LIST_ROAD_SEG Complete')

#Length before adding WC's = 682262

print('Length Before Adding WCs',len(list(H.nodes())))

## ADD EDGE TO H BETWEEN WARD CENTER AND NEAREST ROAD NODE IN H

import numpy as np
from shapely.geometry import Point, LineString
from scipy.spatial import distance

def closest_node(node, nodes):
    closest_index = distance.cdist([node], nodes).argmin()
    return nodes[closest_index]

l = range(34144) #Total of 34144 features

for node in ALL_coord_pair_wc:
    point = Point(node)
    cutting_list = []
    longest_distance = 10000000000000000000000000 # ~ = infinity
    c = 0 #c will be the index of road segment selected
    
    #Parse list of road segments and update until closest one selected (c = index)
    for a in l:
        road_seg_a = LIST_ROAD_SEG[a]
        line_a = LineString(road_seg_a)
        dist_to_road = point.distance(line_a) #how far is the road from this wc?
        if dist_to_road < longest_distance:
            longest_distance = dist_to_road
            c = a #this gives us index number of the road to access for adding point
    
    #Get linestring for the closest road segment
    road_seg_c = LIST_ROAD_SEG[c]
    line_t = LineString(road_seg_c)
    
    #Get nearest point on line to point
    break_node_temp = line_t.interpolate(line_t.project(point))
    coords = list(break_node_temp.coords)
    
    #Turn point into list to match graph data format
    break_node = []
    break_node.append(coords[0][0])
    break_node.append(coords[0][1])
        
    ## Sometimes it picks break nodes that aren't in the linestring (but extremely close)
    ## This corrects for the slight error
    
    if break_node in road_seg_c:
        i = road_seg_c.index(break_node) #i will be the index of node that is break point
        node_on_road = tuple(break_node) #so you can add edge
    else:
        node_on_road = tuple(closest_node(break_node, road_seg_c))
    
    #Get distance between point and nearest point on road
    line_seg_to_node = []
    line_seg_to_node.append(node_on_road)
    line_seg_to_node.append(node)
    line_to_node = LineString(line_seg_to_node)
    L_to_node = line_to_node.length
    
    #Add edge between points to H with weight = distance
    H.add_edge(node_on_road,node,distance=L_to_node)
    
    #Get index of new node in H
    updated_H = list(H.nodes())
    num_new = updated_H.index(node)
    
    #Change the id on this node to reflect the type of node (ward center)
    H.node[node]['identifier'] = 1

print('Ward Centers added to H')  

#Length after adding WC's = 685906 (3644 added)

print('Length After Adding WCs',len(list(H.nodes())))

############################################################
##                        PART 3
############################################################
##               Add all HF's to graph H
############################################################

## Start by getting list of all health facilities coordinates
ALL_coord_pair_hf = []

l = range(6182) #Total of 6182 features
for a in l:
    file = ogr.Open("/Users/Luke/Documents/Tanzania_GIS/2018-07-02_ReUTM/health_facilities_aUREmdl_UTM.shp")  #CHANGE FILE PATH
    shape = file.GetLayer(0)
    feature = shape.GetFeature(a)
    
    if feature == None:
        coord_pair=[0,0]

    else:
        coord_pair = json.loads(feature.ExportToJson())['geometry']['coordinates']
    
    a = tuple(coord_pair) #Converts coordinate pair from list to tuple
    
    ALL_coord_pair_hf.append(a)

#Length before adding HF's = 685906

print('Length Before Adding HFs',len(list(H.nodes())))

## ADD EDGE TO H BETWEEN HEALTH FACILITY AND NEAREST ROAD NODE IN H

import numpy as np
from shapely.geometry import Point, LineString
from scipy.spatial import distance

def closest_node(node, nodes):
    closest_index = distance.cdist([node], nodes).argmin()
    return nodes[closest_index]

l = range(34144) #Total of 34144 features

for node in ALL_coord_pair_hf:
    point = Point(node)
    cutting_list = []
    longest_distance = 10000000000000000000000000 # ~ = infinity
    c = 0 #c will be the index of road segment selected
    
    #Parse list of road segments and update until closest one selected (c = index)
    for a in l:
        road_seg_a = LIST_ROAD_SEG[a]
        line_a = LineString(road_seg_a)
        dist_to_road = point.distance(line_a) #how far is the road from this wc?
        if dist_to_road < longest_distance:
            longest_distance = dist_to_road
            c = a #this gives us index number of the road to access for adding point
    
    #Get linestring for the closest road segment
    road_seg_c = LIST_ROAD_SEG[c]
    line_t = LineString(road_seg_c)
    
    #Get nearest point on line to point
    #hf_to_new_dist = point.distance(line_t) THIS MIGHT BE NEEDED??
    break_node_temp = line_t.interpolate(line_t.project(point))
    coords = list(break_node_temp.coords)
    
    #Turn point into list to match graph data format
    break_node = []
    break_node.append(coords[0][0])
    break_node.append(coords[0][1])
        
    ## Sometimes it picks break nodes that aren't in the linestring (but extremely close)
    ## This corrects for the slight error
    
    if break_node in road_seg_c:
        i = road_seg_c.index(break_node) #i will be the index of node that is break point
        node_on_road = tuple(break_node) #so you can add edge
    else:
        node_on_road = tuple(closest_node(break_node, road_seg_c))
    
    #Get distance between point and nearest point on road
    line_seg_to_node = []
    line_seg_to_node.append(node_on_road)
    line_seg_to_node.append(node)
    line_to_node = LineString(line_seg_to_node)
    L_to_node = line_to_node.length
    
    #Add edge between points to H with weight = distance
    H.add_edge(node_on_road,node,distance=L_to_node)
    
    #Get index of new node in H
    updated_H = list(H.nodes())
    num_new = updated_H.index(node)
    
    #Change the id on this node to reflect the type of node (ward center)
    H.node[node]['identifier'] = 2

print('Health Facilities added to H')

#Length after adding HF's = 692063 (6157 added... some missing)

print('Length After Adding HFs',len(list(H.nodes())))

############################################################
##                        PART 4
############################################################
##            Build Graph F with just WC and HF
############################################################

# Create list of ward center nodes
ward_centers = []

for (p, d) in H.nodes(data=True):
    if d['identifier'] == 1:
        ward_centers.append(p)

print('Total ward centers',len(ward_centers))

# Create list of health facility nodes
health_facilities = []

for (p, d) in H.nodes(data=True):
    if d['identifier'] == 2:
        health_facilities.append(p)

print('Total health facilities',len(health_facilities))

# Check data format

print('First 3 WC',ward_centers[:3]) #Just to see data format

## BUILD GRAPH F BY TAKING SHORTEST PATH BETWEEN EACH WC AND HF

# Example of obtaining one shortest path between HF and WC

a = ward_centers[1110]
b = health_facilities[50]
p_length = nx.shortest_path_length(H, source=a, target=b, weight='distance')
print('Shortest Path Length between WC 1110 and HF 50',p_length)

# Find nearest neighboring wc to each hf, store in F

F = nx.Graph() # F will be graph of just wc and hf

for a in health_facilities:
    for b in ward_centers:
        try:
            p_length = nx.shortest_path_length(H, source=a, target=b, weight='distance')
            F.add_edge(a,b,distance=p_length)
        except:
            pass
        
print('Finished Building F')
