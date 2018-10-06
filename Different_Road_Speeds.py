############################################################
##                        PART 1
############################################################
## Write all road edges to graph H with weight = distance
############################################################

## Write all road edges to graph H with weight as distance

import networkx as nx
import os
import pickle
import fiona
from osgeo import ogr
import json

##Variables for ease
try:
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
except:
    wd = "E:\\GIT_Checkouts\\Python\\Tanzania-Optimization\\"
    os.chdir(wd)
Hexport1 = "./Temp/H_part1.nx"
Hexport2 = "./Temp/H_part2.nx"
Hexport3 = "./Temp/H_part3.nx"
roadsExport = "./Temp/LIST_ROAD_SEG.p"
matrix_file = "/Users/Luke/Documents/Tanzania_GIS/matrix_file.p" # COME BACK AND FIX THIS!!

# Start with using nx.read_shp to generate G
# Set simplify = False so that you get every road node, not just road segment ends

# SPEED LIMITS (meters/hour)- can change later
primary_speed = 80000
secondary_speed = 40000
tertiary_speed = 30000
walking_speed = 4000

from shapely.geometry import LineString #need for line.length

# Convert G into edge list with weighted edges
if os.path.isfile(Hexport3): #checks if the final graph is created
    print("Files processed Skipping to Part 4")
elif os.path.isfile(Hexport1):
    print("Loading saved part 1 graph")
    H = nx.read_gpickle(Hexport1)
else:
    print("Calculating Part 1 Graph")
     
    ## ## ## ## ## ## ## ## ## 
    ## Need to create 3 shapefiles that will become primary, secondary, tertiary roads
    ## ## ## ## ## ## ## ## ## 
    
    #### First create three of the same shapefiles
    with fiona.open('/Users/Luke/Documents/Tanzania_GIS/shp_for_dist_matrix/roads_utm.shp') as source:
        source_schema = source.schema
        source_driver = source.driver
        source_crs = source.crs
        print(source_schema) # attribute fields & geometry def as dict
        print(source_driver) # "ESRI Shapefile"
        print(source_crs) # coordinate system

        with fiona.open('/Users/Luke/Documents/Tanzania_GIS/shp_for_dist_matrix/roads_PRIMARY.shp', 'w',
                        driver=source_driver,
                        crs=source_crs,
                        schema=source_schema) as shpout:
            for feature in source:
                # if feature should be written:
                shpout.write(feature)

        with fiona.open('/Users/Luke/Documents/Tanzania_GIS/shp_for_dist_matrix/roads_SECONDARY.shp', 'w',
                        driver=source_driver,
                        crs=source_crs,
                        schema=source_schema) as shpout:
            for feature in source:
                # if feature should be written:
                shpout.write(feature)

        with fiona.open('/Users/Luke/Documents/Tanzania_GIS/shp_for_dist_matrix/roads_TERTIARY.shp', 'w',
                        driver=source_driver,
                        crs=source_crs,
                        schema=source_schema) as shpout:
            for feature in source:
                # if feature should be written:
                shpout.write(feature)
                
    #### Next make them match their name(sort by road type)
    ## Primary
    file_roads = ogr.Open('/Users/Luke/Documents/Tanzania_GIS/shp_for_dist_matrix/roads_PRIMARY.shp',True) #CHANGE FILE PATH
    shape_roads = file_roads.GetLayer(0)
    r = range(len(shape_roads))

    print("Total roads: {}".format(shape_roads.GetFeatureCount()))

    for a in r:  
        feature_roads = shape_roads.GetFeature(a)
        first = json.loads(feature_roads.ExportToJson())['properties']['type']
        if (first == 'primary') or (first == 'primary_link') \
            or (first == 'trunk') or (first == 'trunk_link'):
            pass
        else:
            shape_roads.DeleteFeature(a)

    file_roads.ExecuteSQL('REPACK ' + shape_roads.GetName())
    file_roads.ExecuteSQL('RECOMPUTE EXTENT ON ' + shape_roads.GetName())

    print("Primary roads: {}".format(shape_roads.GetFeatureCount()))

    del file_roads
    
    ## Secondary
    file_roads = ogr.Open('/Users/Luke/Documents/Tanzania_GIS/shp_for_dist_matrix/roads_SECONDARY.shp',True) #CHANGE FILE PATH
    shape_roads = file_roads.GetLayer(0)
    r = range(len(shape_roads))

    print("Total roads: {}".format(shape_roads.GetFeatureCount()))

    for a in r:  
        feature_roads = shape_roads.GetFeature(a)
        first = json.loads(feature_roads.ExportToJson())['properties']['type']
        if (first == 'secondary') or (first == 'secondary_link'):
            pass
        else:
            shape_roads.DeleteFeature(a)

    file_roads.ExecuteSQL('REPACK ' + shape_roads.GetName())
    file_roads.ExecuteSQL('RECOMPUTE EXTENT ON ' + shape_roads.GetName())

    print("Secondary roads: {}".format(shape_roads.GetFeatureCount()))

    del file_roads
    
    ## Tertiary
    file_roads = ogr.Open('/Users/Luke/Documents/Tanzania_GIS/shp_for_dist_matrix/roads_TERTIARY.shp',True) #CHANGE FILE PATH
    shape_roads = file_roads.GetLayer(0)
    r = range(len(shape_roads))

    print("Total roads: {}".format(shape_roads.GetFeatureCount()))

    for a in r:  
        feature_roads = shape_roads.GetFeature(a)
        first = json.loads(feature_roads.ExportToJson())['properties']['type']
        if (first == 'road') or (first == 'residential') \
            or (first == 'tertiary_link') or (first == 'tertiary') \
            or (first == 'living_street') or (first == 'dirtroad'):
            pass
        else:
            shape_roads.DeleteFeature(a)

    file_roads.ExecuteSQL('REPACK ' + shape_roads.GetName())
    file_roads.ExecuteSQL('RECOMPUTE EXTENT ON ' + shape_roads.GetName())

    print("Tertiary roads: {}".format(shape_roads.GetFeatureCount()))

    del file_roads
    
    #### Finally, add the roads to H step by step, dividing distance by travel time
    ## First, add primary roads to H
    G_primary=nx.read_shp('/Users/Luke/Documents/Tanzania_GIS/shp_for_dist_matrix/roads_PRIMARY.shp', simplify=False) #CHANGE FILE PATH
    H = nx.Graph() #H is graph where new edges will be written to
    
    for u,v in G_primary.edges():
        n1 = u
        n2 = v
        
        #Need to have nodes as list
        line_seg = [] 
        line_seg.append(n1)
        line_seg.append(n2)
        
        #Length calculations
        line = LineString(line_seg)
        L = (line.length)/primary_speed #note addition of primary_speed
        #if L < 5:
        #    print("Shorty")
        
        #Create new weighted graph
        H.add_edge(n1,n2,distance=L)
        
    ## Next, add secondary roads to H
    G_secondary=nx.read_shp('/Users/Luke/Documents/Tanzania_GIS/shp_for_dist_matrix/roads_SECONDARY.shp', simplify=False) #CHANGE FILE PATH
    
    for u,v in G_secondary.edges():
        n1 = u
        n2 = v
        
        #Need to have nodes as list
        line_seg = [] 
        line_seg.append(n1)
        line_seg.append(n2)
        
        #Length calculations
        line = LineString(line_seg)
        L = (line.length)/secondary_speed #note addition of secondary_speed
        #if L < 5:
        #    print("Shorty")
        
        #Create new weighted graph
        H.add_edge(n1,n2,distance=L)
        
    ## Finally, add tertiary roads to H
    G_tertiary=nx.read_shp('/Users/Luke/Documents/Tanzania_GIS/shp_for_dist_matrix/roads_TERTIARY.shp', simplify=False) #CHANGE FILE PATH
    
    for u,v in G_tertiary.edges():
        n1 = u
        n2 = v
        
        #Need to have nodes as list
        line_seg = [] 
        line_seg.append(n1)
        line_seg.append(n2)
        
        #Length calculations
        line = LineString(line_seg)
        L = (line.length)/tertiary_speed #note addition of secondary_speed
        #if L < 5:
        #    print("Shorty")
        
        #Create new weighted graph
        H.add_edge(n1,n2,distance=L)
        
    ## Note identifier key: 0 = road, 1 = ward center, 2 = health facility
    nx.set_node_attributes(H, 0, 'identifier') #Want to know what kind of node it is (road)
    nx.write_gpickle(H,Hexport1)
    print('First 3 Edges',list(H.edges(data=True))[:3]) #See data format
    print('Distance of first road seg',list(H.edges(data=True))[0][2]['distance']) #Accessing distance
    print('Number of road segments',len(list(H.edges()))) #Should be 700881 so 55 short...
    print('Identifer for first node',list(H.nodes(data=True))[0][1]['identifier']) #Accessing identifier


## PAIR EACH WARD CENTER WITH ITS POPULATION

    from osgeo import ogr
    import json
    
    ALL_coord_pair_wc = []
    ALL_pop_wc = []
    
    file_cent = ogr.Open("/Users/Luke/Documents/Tanzania_GIS/2018-07-02_ReUTM/TZ_CentroidsWith2012Pop_UTM.shp") #CHANGE FILE PATH
    shape_cent = file_cent.GetLayer(0)
    l = range(len(shape_cent)) #Total of 3644 ward centers (wc's)

    for a in l:
        feature_cent = shape_cent.GetFeature(a)
        
        if feature_cent == None:
            coord_pair=[0,0]
            ward_pop = 0
    
        else:
            coord_pair = json.loads(feature_cent.ExportToJson())['geometry']['coordinates']
            ward_pop = json.loads(feature_cent.ExportToJson())['properties']['total_both']
        
        a = tuple(coord_pair) #Converts coordinate pair from list to tuple
        b = ward_pop
        
        ALL_coord_pair_wc.append(a)
        ALL_pop_wc.append(b)
    
    wc_pop_dict = list(zip(ALL_coord_pair_wc,ALL_pop_wc)) # Completes pairing ward center with pop
    
    #print('1st WC Population,'(wc_pop_dict)[0][1]) #Gives you population for 1st ward center

###########################################################
#                        PART 2
###########################################################
#               Add all WC's to graph H
###########################################################

# Create shapefile of "roads" that does not contain ones we consider walking

file_roads = ogr.Open("/Users/Luke/Documents/Tanzania_GIS/road_to_delete/roads_utm.shp",True) #CHANGE FILE PATH
shape_roads = file_roads.GetLayer(0)
r = range(len(shape_roads))

print("Features before: {}".format(shape_roads.GetFeatureCount()))

for a in r:  
    feature_roads = shape_roads.GetFeature(a)
    first = json.loads(feature_roads.ExportToJson())['properties']['type']
    if (first == 'bridleway') or (first == 'construction') \
        or (first == 'cycleway') or (first == 'footway') \
        or (first == 'path') or (first == 'pedestrian') \
        or (first == 'proposed')  or (first == 'service') \
        or (first == 'steps') or (first == 'track') \
        or (first == 'unclassified') or (first == 'unknown') \
        or (first == 'unsurfaced'):
        shape_roads.DeleteFeature(a)
    else:
        pass

file_roads.ExecuteSQL('REPACK ' + shape_roads.GetName())
file_roads.ExecuteSQL('RECOMPUTE EXTENT ON ' + shape_roads.GetName())

print("Features after: {}".format(shape_roads.GetFeatureCount()))

del file_roads

# Store list of linestrings of road segments to access faster later
print("Starting Part 2")


if os.path.isfile(roadsExport):
    print("Opening preexisting layer of linestrings and road segments")
    with open (roadsExport, 'rb') as fp:
        LIST_ROAD_SEG = pickle.load(fp)
else:
    print("creating list of linestrings and road segments")
    from osgeo import ogr
    import json
    LIST_ROAD_SEG = []
    file_roads = ogr.Open("/Users/Luke/Documents/Tanzania_GIS/road_to_delete/roads_utm.shp") #CHANGE FILE PATH
    shape_roads = file_roads.GetLayer(0)
    count = 0
    l = range(len(shape_roads)) #Did this programatically
    for a in l:
        count +=1
        feature_roads = shape_roads.GetFeature(a)
        first = json.loads(feature_roads.ExportToJson())['geometry']['coordinates']
        LIST_ROAD_SEG.append(first)
        if (count % 5000) == 0 : #and (used > 0):
            print(count)
    with open(roadsExport, 'wb') as fp:
        pickle.dump(LIST_ROAD_SEG, fp)
    print('LIST_ROAD_SEG Complete')

#Length before adding WC's = 682262

    print('Length Before Adding WCs',len(list(H.nodes())))

## ADD EDGE TO H BETWEEN WARD CENTER AND NEAREST ROAD NODE IN H

if os.path.isfile(Hexport3):
    pass
elif os.path.isfile(Hexport2):
    print("Loading saved part 2 graph")
    H = nx.read_gpickle(Hexport2)
else:
    from shapely.geometry import Point, LineString
    from scipy.spatial import distance
    
    def closest_node(node, nodes):
        closest_index = distance.cdist([node], nodes).argmin()
        return nodes[closest_index]
    print("Adding edges to H between ward center and nearest road node")
    l = range(len(shape_roads)) #Did this programatically
    count = 0
    for node in ALL_coord_pair_wc:
        count+=1
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
        L_to_node = (line_to_node.length)/walking_speed #Note addition of walking_speed
        
        #Add edge between points to H with weight = distance
        H.add_edge(node_on_road,node,distance=L_to_node)
        
        #Get index of new node in H
        updated_H = list(H.nodes())
        num_new = updated_H.index(node)
        
        #Change the id on this node to reflect the type of node (ward center)
        H.node[node]['identifier'] = 1
        if (count % 1000) == 0 : #and (used > 0):
            print(count)
    nx.write_gpickle(H,Hexport2)
#
    print('Ward Centers added to H')  

#Length after adding WC's = 685906 (3644 added)

    print('Length After Adding WCs',len(list(H.nodes())))

############################################################
##                        PART 3
############################################################
##               Add all HF's to graph H
############################################################

## Start by getting list of all health facilities coordinates
print("Part 3")

#Length before adding HF's = 685906

#print('Length Before Adding HFs',len(list(H.nodes())))

## ADD EDGE TO H BETWEEN HEALTH FACILITY AND NEAREST ROAD NODE IN H


if os.path.isfile(Hexport3):
    print("Loading saved part 3 graph")
    H = nx.read_gpickle(Hexport3)
else:
    print("Getting list of all heath facility coordinates")
    import numpy as np
    from shapely.geometry import Point, LineString
    from scipy.spatial import distance

    def closest_node(node, nodes):
        closest_index = distance.cdist([node], nodes).argmin()
        return nodes[closest_index]
    ALL_coord_pair_hf = []
    file_facilities = ogr.Open("/Users/Luke/Documents/Tanzania_GIS/2018-07-02_ReUTM/health_facilities_aUREmdl_UTM.shp")  #CHANGE FILE PATH
    shape_facilities = file_facilities.GetLayer(0)
    l = range(len(shape_facilities)) #Total of 6182 features
    for a in l:
        feature_facilities = shape_facilities.GetFeature(a)        
        if feature_facilities == None:
            coord_pair=[0,0]    
        else:
            coord_pair = json.loads(feature_facilities.ExportToJson())['geometry']['coordinates']        
        a = tuple(coord_pair) #Converts coordinate pair from list to tuple       
        ALL_coord_pair_hf.append(a)  
    print("Adding edges to H between ward center and nearest road node")
    file_roads = ogr.Open("/Users/Luke/Documents/Tanzania_GIS/road_to_delete/roads_utm.shp") #CHANGE FILE PATH
    shape_roads = file_roads.GetLayer(0)
    l = range(len(shape_roads)) #Total of 34144 features
    count = 0
    for node in ALL_coord_pair_hf:
        count +=1
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
        L_to_node = (line_to_node.length)/walking_speed #Note addition of walking_speed
        
        #Add edge between points to H with weight = distance
        H.add_edge(node_on_road,node,distance=L_to_node)
        
        #Get index of new node in H
        updated_H = list(H.nodes())
        num_new = updated_H.index(node)
        
        #Change the id on this node to reflect the type of node (ward center)
        H.node[node]['identifier'] = 2
        if (count % 1000) == 0 : #and (used > 0):
            print(count)
    nx.write_gpickle(H,Hexport3)

print('Health Facilities added to H')

#Length after adding HF's = 692063 (6157 added... some missing)

print('Length After Adding HFs',len(list(H.nodes())))

############################################################
##                        PART 4
############################################################
##            Build Graph F with just WC and HF
############################################################
Fexport = "F_temp.nx"
list_done_file = "list_done.p"
# Create list of ward center nodes
print("Part 4")
print("Creating list of ward center nodes")
ward_centers = []

for (p, d) in H.nodes(data=True):
    if d['identifier'] == 1:
        ward_centers.append(p)

print('Total ward centers',len(ward_centers))

# Create list of health facility nodes
print("Creating list of health facility nodes")

health_facilities = []

for (p, d) in H.nodes(data=True):
    if d['identifier'] == 2:
        health_facilities.append(p)

print('Total health facilities',len(health_facilities))

# Check data format

print('First 3 WC',ward_centers[:3]) #Just to see data format


############################################################
##                        PART 5
############################################################
##            Build Graph with just WC and HF
############################################################

# A) Convert H to igraph

import igraph as ig

nx.write_gml(H,'graph22.gml', stringizer= nx.readwrite.gml.literal_stringizer) # Export NX graph to file

F_temp = ig.read('graph22.gml',format="gml") # Create new IG graph from file

### Accessing Data in igraph/ Checking Conversion from NetworkX
#print('Checking number of edges',F_temp.ecount()) #Checking
#print('Checking node identifiers',F_temp.vs[692060]['identifier']) #Should be 2 for this node
#print('Checking if length is stored for edges',F_temp.es[0]['distance']) #Should be ~145.06
#print('Checking if nodes preserved as coordinates',F_temp.vs[692060]['id']) #Accessing node coordinates

# B) Lists to iterate over

ward_centers = []
health_facilities = []

for v in F_temp.vs:
    #print(v['identifier'])
    if v['identifier'] == 1.0:
        ward_centers.append(v)
    if v['identifier'] == 2.0:
        health_facilities.append(v)

### Checking Lists      
#print('Number ward centers:',len(ward_centers))
#print('First 3 ward centers',ward_centers[:3])
#print('Number health facilities:',len(health_facilities))
#print('First 3 health facilities',health_facilities[:3])

# C) Shortest path/ Build F

F = ig.Graph() #Create empty igraph

c = 1 #for adding health facilities
d = 2 #for adding ward centers

z = 0 #counter
import time

if os.path.isfile(matrix_file):
    with open (matrix_file, 'rb') as fp:
        path_matrix = pickle.load(fp)
    
else:
    start= time.time()
    path_matrix = F_temp.shortest_paths_dijkstra(source=health_facilities,target=ward_centers,weights='distance', mode = "all")
    end = time.time()
    print(end - start)
    with open(matrix_file, 'wb') as fp:
        pickle.dump(path_matrix, fp)

# Note: Rows are labeled with HF, columns with WC
