# Imports
import os
import numpy as np
import math
import osmnx as ox
import networkx as nx
import pandas as pd
import random
import matplotlib.pyplot as plt
import community as com
from networkx.algorithms import community
import wntr
from heapq import nlargest
import geopy.distance
from shapely import geometry, ops
from shapely.geometry import Point, Polygon
from pyproj import Proj
import utm
from matplotlib.pyplot import cm
import wntr.network.controls as controls

# Class definition
class RandomWaterDistributionNetwork:

    def __init__(self, gmap_key, world_cities, min_population, distance, demand_values, roughness_values, reservoir_heads, pipes_diameters):

        '''
        This class is for constructing Random Virtual Water Distribution Networks from OpenStreetMap streets layout using osmnx module, and using Wntr module for diameter sizing   


        gmap_key(str): key for google maps API for downloading the elevation values for the nodes

        world_cities(pd.DataFrame): pandas DataFrame from the free csv file in: https://simplemaps.com/data/world-cities , it is used for
        extracting the coordinates of the cities.

        min_population(int): the minimal population for the cities to be used for the generation of the Virtual Water Network

        distance(int): is the distance of the nodes from the center point that are going to be kept for constructing the graph

        demand_values(dict): a dictionnary where the keys are the type of soil occupation, and values are peak hour demands in m3/h/hectar

        roughness_values(list): is a list of roughness values that will be assigned randomly to pipes

        reservoir_heads(list): is a list of heads that will be assigned randomly to the reservoirs 

        pipe_diameters(list): list of possible diameters of pipes in the networks (in m)
        '''

        self._gmap_key = gmap_key
        self._world_cities = world_cities
        self._min_population = min_population
        self._distance = distance
        self._demand_values = demand_values
        self._roughness_values = roughness_values
        self._reservoir_heads = reservoir_heads
        self._pipes_diameters = pipes_diameters
        # A list to strore problematic pipes (pipes where the start node is also the end node)
        self._problematic_pipes = []
        self._number_of_reservoirs = 2
        self._main_valves={}
        self._main_pipes=[]
        self._random_point = (32.933699, -5.666542)
        self._center = None
        self._out_points = None
        self._velocity = None
        self._pressure = None
        self._valves_per_sect = {}

    ###################################################################################################################
    # PUBLIC METHODS AND PROPERTIES
    ###################################################################################################################

    def generate_random_graph(self):
        '''
        Returns a random networkx graph using the osmnx module
        '''
        # Filter the world_cities dataframe to keep just the cities with population over _min_population
        big_cities = self._world_cities[self._world_cities['population'] > self._min_population]
        # Choose a random city coordinates
        random_index = random.randint(0, len(big_cities))
        random_city_coord = (big_cities.loc[random_index, 'lat'], big_cities.loc[random_index, 'lng'])
        # Add randomness to the coordinates
        r1 = random.randint(-100, 100)/5000
        r2 = random.randint(-100, 100)/5000
        self._random_point = (random_city_coord[0]+r1, random_city_coord[1]+r2)
        # Generate the street network for the bonding box that have center = _random_point and side length = 2*_distance
        return ox.graph_from_point(self._random_point, distance=self._distance, network_type='drive')


    def add_elevation(self, G):
        '''
        Add elevation attribute to the nodes of the input graph
        '''
        # Add elevation using osmnx module method
        G = ox.add_node_elevations(G, api_key=self._gmap_key)
        

    def clean_graph(self, G, tol): 
        '''
        Merges the clustered nodes that are less than tol apart 
        and deletes self loops and parallel edges from the input graph
        partial(bool): True if used for cleaning the main network subgraph,
        in this case the graph also is cleaned to be compatible with the subgraph

        tol(int): the distance within which the clustered nodes would be merged
        '''      
        # Merge clustered node
        counter = 0            
        all_cluster = []
        nodes_list = G.nodes() # we need to store the list of nodes, because the list will be changing inside the loop
        k = 1
        iteration = 0
        while k>0 and iteration < 10:
            k=0
            iteration += 1
            for n in nodes_list:
                previous_nodes = [y for x in all_cluster for y in x]
                if n not in previous_nodes:
                    cluster = self.tooclose(G, n, tol)
                    if len(cluster)>1:
                        k+=1
                        all_cluster.append(cluster)
                        for x in cluster[:len(cluster)-1]:
                            if x not in previous_nodes:
                                if G.has_node(x):
                                    G = nx.contracted_nodes(G, cluster[len(cluster)-1], x, self_loops=False)
                                    counter +=1
                                    previous_nodes.append(x)        
        nodes_to_remove = []
        # Remove parallel edges and self loops:
        edges_set = list(G.edges())
        for (u, v) in edges_set:
            if G.edges[u, v]['length'] > 1000:
                G.remove_edge(u, v)
            elif G.edges[u, v]['length'] > 150:
                if len(list(G.neighbors(u))) < 2:
                    nodes_to_remove.append(u)
                if len(list(G.neighbors(v))) < 2:
                    nodes_to_remove.append(v)                    
            if u == v:
                while G.number_of_edges(u, v)>0:
                    G.remove_edge(u, v)
            while G.number_of_edges(u, v)>1:
                G.remove_edge(u, v)
        G = self.main_connected(G)
        return G


    def add_node_demands(self, G, area=0.2):
        '''
        Add demand attribute to the nodes of the input graph. The graph is divided to communities and each community is assigned a demand value randomly from demand_values dictionnary

        area(float): The average area in hectare covered by each of the nodes
        '''
        G_copy = nx.Graph(G)
        n_nodes = len(G)
        n_communities = round(n_nodes/4)
        # Using fluid communties algorithm for creating 4 nodes clusters
        community_generator = community.asyn_fluidc(G_copy, n_communities)
        demands = {}
        # Converting demand values to peak demands in m3/s
        demand_values = {key: 4*area*0.00379*v/(3600*24) for key,v in self._demand_values.items()}           
        for i in range(n_communities):
            # Construct a dictionary of demands for each community
            # (community nodes are the keys and demand values are the values)
            demands_list = list(demand_values.values())
            demands.update(dict.fromkeys(next(community_generator), {'demand' : demands_list[i%len(demand_values)]}))
        # Add the demand as attributes to the nodes in the graph
        nx.set_node_attributes(G, demands)


    def main_network(self, G):
        '''
        Create the main distribution network
        '''
        zones = self.zonage(G)
        edges_list = self.connect_zones(G, zones)
        subG = nx.Graph(G).edge_subgraph(edges_list).copy()
        return subG


    def clean_cycles(self, G):
        '''
        Cleans the input graph G from double lines, in these steps:
        1- small cycles where the nodes are aligned
        2- parallel edges and self loops
        '''
        # 1- clean small cycles where the nodes are aligned
        for c in range(2):
            for node in G.nodes():
                nodedata = G.nodes[node]
                nodex, nodey = nodedata['x'], nodedata['y']
                neighbors = G.neighbors(node)
                for neigh in list(neighbors):                
                    neighdata = G.nodes[neigh]                
                    neix, neiy = neighdata['x'], neighdata['y']
                    d = self.proj_distance(neix, neiy, nodex, nodey)
                    cluster = self.tooclose(G, node, d)
                    nodes_in_way = self.in_way(G, node, neigh, cluster)
                    if len(list(nodes_in_way.values())) > 0:
                        G.remove_edge(node, neigh)
                        for k in range(len(nodes_in_way)):
                            current = int(list(nodes_in_way)[k])
                            if k == 0: 
                                G.add_edge(node, current)
                                G.edges[node, current]['length'] = list(nodes_in_way.values())[k]
                            if k == len(nodes_in_way)-1:
                                G.add_edge(current, neigh)
                                G[current][neigh]['length'] = list(nodes_in_way.values())[k]
                            if 0 < k < len(nodes_in_way)-1:
                                nnext = int(list(nodes_in_way)[k+1])
                                G.add_edge(current, nnext)
                                G[current][nnext]['length'] = list(nodes_in_way.values())[k+1]-list(nodes_in_way.values())[k]        
        # 2- remove parallel edges and self loops
        nodes_to_remove = []
        for u in G.nodes():
            for v in list(G.neighbors(u)):
                for n in list(G.neighbors(u)):
                    if (G.has_edge(v,n) and len(list(G.neighbors(n))) < 3):
                        G.remove_edge(v,n)
                        G.remove_edge(u,n)
                        nodes_to_remove.append(n)
        for n in nodes_to_remove:
            G.remove_node(n)
        G = self.main_connected(G)     
        return G


    def add_reservoirs(self, G):
        '''
        Adds number_of_reservoirs reservoirs to the graph G
        '''        
        n_edges = G.number_of_edges()
        self._number_of_reservoirs = max(round(n_edges/100), round(self._distance/500), 2)        
        n_communities = self._number_of_reservoirs
        community_generator = community.asyn_fluidc(self.main_connected(G), n_communities)
        for i in range(n_communities):
            # Make a subgraph from the community nodes
            subG = G.subgraph(next(community_generator))
            # Get the maximum elevation node from the subgraph
            max_elevation_node = max(dict(subG.nodes).items(), key=lambda x: x[1]['elevation'])[0]
            Reservoir_attr = {'x': G.nodes[max_elevation_node]['x']+ random.randint(10, 30), 'y': G.nodes[max_elevation_node]['y'] + random.randint(10, 30), 'elevation': random.choice(self._reservoir_heads) + G.nodes[max_elevation_node]['elevation'], 'demand': 0}
            G.add_node('Reservoir{}' .format(i), **Reservoir_attr)
            G.add_edge('Reservoir{}' .format(i), max_elevation_node)
            G['Reservoir{}'.format(i)][max_elevation_node]['length'] = 50
            G['Reservoir{}'.format(i)][max_elevation_node]['roughness'] = random.choice(self._roughness_values)
        return G


    def add_edge_roughness(self, G):
        '''
        Adds the roughness as attribute to the edges in the graph
        '''
        nx.set_edge_attributes(G, random.choice(self._roughness_values), 'roughness')  


    def create_wn(self, G, subG):
        '''
        Creates wntr.network.WaterNetworkModel for the graph G
        '''
        wn = wntr.network.WaterNetworkModel()
        G_junctions = G.copy()
        # main pipes diameters
        main_pipes = nlargest(3, self._pipes_diameters)
        for i in range(self._number_of_reservoirs):
            G_junctions.remove_node('Reservoir{}' .format(i))
        n_nodes = len(G_junctions)
        n_edges = G.number_of_edges()
        nodes_list = list(G_junctions.nodes(data=True))
        edges_list = list(G.edges(data=True))
        # Adding demand pattern
        wn.add_pattern('pat1', [0.56, 0.37, 0.28, 0.33, 0.45, 0.67, 1.13, 2.02, 1.78, 1.6, 1.56, 1.45, 1.18, 1.19, 1.07, 0.97, 0.67, 0.9, 0.82, 1.22, 1.57, 1.3, 1.05, 0.92])
        # Adding the nodes: the demand junctions and reservoirs
        for i in range(n_nodes):            
            node = nodes_list[i]
            wn.add_junction('{}'.format(node[0]), base_demand=node[1]['demand'], elevation=node[1]['elevation'], coordinates=(node[1]['x'], node[1]['y']), demand_pattern='pat1')
        for i in range(self._number_of_reservoirs):
            reservoir = G.nodes['Reservoir{}' .format(i)]
            wn.add_reservoir('Reservoir{}' .format(i), base_head=reservoir['elevation'], coordinates=(reservoir['x'], reservoir['y']))   
        # Adding edges
        connected_nodes = []
        for i in range(n_edges):
            edge = edges_list[i]
            if edge[0] == edge[1] or edge[2]['length'] > 2000: self._problematic_pipes.append(i)
            if subG.has_edge(edge[0], edge[1]):
                diam = random.choice(main_pipes)
                wn.add_pipe('{}' .format(i), edge[0], edge[1], length=edge[2]['length'], diameter=0.6, roughness=edge[2]['roughness'], minor_loss=0.0, status='OPEN')
                connected_nodes.append('{}' .format(edge[0]))
                connected_nodes.append('{}' .format(edge[1])) 
                self._main_pipes.append(i)                
            else:
                wn.add_pipe('{}' .format(i), edge[0], edge[1], length=edge[2]['length'], diameter=0.2, roughness=edge[2]['roughness'],
                minor_loss=0.0, status='OPEN')
                connected_nodes.append('{}' .format(edge[0]))
                connected_nodes.append('{}' .format(edge[1]))
        disconnected_nodes = set(wn.node_name_list)-set(connected_nodes)
        for i in disconnected_nodes:
            wn.remove_node('{}'.format(i))           
        return wn


    def pipe_sizing(self, wn):
        '''
        Adjusting the diameters of the pipes to get realistic values for the pressures and velocities in the network
        '''
        # Simulate hydraulics
        sim = wntr.sim.EpanetSimulator(wn)
        results = sim.run_sim()
        velocity = results.link['velocity'].values[0]
        pressure = results.node['pressure'].values[0][:-self._number_of_reservoirs].mean()
        min_pressure = results.node['pressure'].values[0][:-self._number_of_reservoirs].min()
        max_velocity = max(velocity)        
        pipes_list = set(range(len(velocity)))-set(self._problematic_pipes)
        # Uncomment the next line if you want to avoid changing the diameters in the main network
        pipes_list = pipes_list-set(self._main_pipes)
        def bigger_diameter(diameter):
            i = 0
            while diameter >= self._pipes_diameters[i] and i < len(self._pipes_diameters)-1 : i+=1
            return self._pipes_diameters[i]

        def smaller_diameter(diameter):
            i = len(self._pipes_diameters)-1
            while diameter <= self._pipes_diameters[i] and i > 0: i-=1
            return self._pipes_diameters[i]

        counter = 0
        starting_diam = self._pipes_diameters[0]
        while counter < 4*len(pipes_list) and max_velocity > 6:
            sim = wntr.sim.EpanetSimulator(wn)
            results = sim.run_sim()
            velocity = results.link['velocity'].values[0]
            max_velocity = max(velocity)
            for link in pipes_list:
                pipe = wn.get_link('{}'.format(link))
                node1 = pipe.start_node_name
                node2 = pipe.end_node_name
                pressure1 = results.node['pressure'][node1][0]
                pressure2 = results.node['pressure'][node2][0]
                pressure = (pressure1 + pressure2)/2
                if pressure < min_pressure: min_pressure = pressure
                if velocity[link] < 0.5 and pressure > 70:
                    pipe.diameter = smaller_diameter(pipe.diameter)
                    counter += 1
                elif velocity[link] > 1.5 and pressure < 40:
                    pipe.diameter = bigger_diameter(pipe.diameter)
                    counter += 1
                elif velocity[link] > 2.5:
                    pipe.diameter = bigger_diameter(pipe.diameter)
                    counter += 1
                elif velocity[link] < 0.01:
                    pipe.diameter = smaller_diameter(pipe.diameter)
                    counter += 1
                elif pressure < 10:
                    pipe.diameter = bigger_diameter(pipe.diameter)
                    counter += 1
                elif pressure > 100:
                    pipe.diameter = smaller_diameter(pipe.diameter)
                    counter += 1 
        self._velocity = results.link['velocity'].values[0]
        self._pressure = results.node['pressure'].values[0][:-self._number_of_reservoirs]
        return wn


    def add_valves(self, wn, subG, n_sectors):
        '''
        Divide the network to sectors using Louvain-Algorithm, and add valves between sectors.
        '''
        G = wn.get_graph()
        G = nx.Graph(G)
        comm = community.asyn_fluidc(G, n_sectors)
        connected_sect = {}
        edges_bw_sect = {}
        sectors = {}
        k=0
        while True:
            try:
                sectors['{}'.format(k)] = list(next(comm))
                k+=1
            except StopIteration:
                break
        
        self.sectors = sectors    
        for sect1 in sectors:
            self._valves_per_sect[sect1] = []
            self._main_valves[sect1] = []
            edges_bw_sect[sect1] = []

        valves = []
        # Find connection between sectors (one for each pair)
        for (u,v) in list(G.edges()):
            for sect1 in sectors:
                if (u in sectors[sect1] and v not in sectors[sect1]) or (v in sectors[sect1] and u not in sectors[sect1]):
                    edges_bw_sect[sect1].append((u,v))
                    for sect2 in sectors:
                        if sect2!=sect1:
                            if v in sectors[sect2]:
                                valves.append((u,v))
        # Add valves
        counter = 0
        pipes = []
        removed_pipes = []
        pipes_list = set(range(wn.num_links))-set(self._problematic_pipes)
        for link in pipes_list:
            pipe = wn.get_link('{}'.format(link))
            node1 = pipe.start_node_name
            if node1.startswith('R'): no1 = node1
            else: no1 = int(node1)
            node2 = pipe.end_node_name
            if node2.startswith('R'): no2 = node2
            else: no2 = int(node2)
            diameter = pipe.diameter         
            if (node1, node2) in valves:
                n1, n2 = None, None
                for sect in sectors:
                    if (node1, node2) in edges_bw_sect[sect] and (n1, n2) != (node1, node2):
                        n1, n2 = node1, node2
                        removed_pipes.append('{}'.format(link))
                        wn.add_valve('S{}_{}'.format(sect, counter), node1, node2, diameter, 'TCV', 0)
                        if subG.has_edge(no1, no2) or subG.has_edge(no2, no1):
                            self._main_valves[sect].append('S{}_{}'.format(sect, counter))
                        else:
                            self._valves_per_sect[sect].append('S{}_{}'.format(sect, counter))
            counter += 1
        # Connect the valves to the closest node in the main distribution network
        min_dist = self._distance
        edges_list = []
        for sect in sectors:
            if len(self._main_valves[sect]) < 1:
                edge_sect = random.choice(edges_bw_sect[sect])
                node_sect = edge_sect[1]
                xr = G.nodes[node_sect]['pos'][0]
                yr = G.nodes[node_sect]['pos'][1]
                closest = list(subG.nodes())[0]
                for n in subG.nodes():
                    xn = subG.nodes[n]['x']
                    yn = subG.nodes[n]['y']
                    dist = self.proj_distance(yn, xn, yr, xr)
                    if dist < min_dist:
                        min_dist = dist
                        closest = n
                path = nx.algorithms.shortest_path(G, node_sect, '{}'.format(closest),'length')
                edges = [edge_sect]
                edges += [(path[j], path[j+1]) for j in range(len(path)-1)]
                edges_list += edges
        # Make the diameters of pipes of the connection to the main d=600mm
        for pipe_name, pipe in wn.pipes():
            node1 = pipe.start_node_name
            node2 = pipe.end_node_name
            if (node1, node2) in edges_list or (node2, node1) in edges_list:
                pipe.diameter = 0.6
        # save the lists of main pipes(control pipes) and pipes between sectors (isolation pipes)
        for sect in self._valves_per_sect:
            for valve_name in self._valves_per_sect[sect]:
                valve = wn.get_link(valve_name)
                node1 = valve.start_node_name
                node2 = valve.end_node_name
                if (node1, node2) in edges_list or (node2, node1) in edges_list:
                    self._valves_per_sect[sect].remove(valve_name)
                    self._main_valves[sect].append(valve_name)  
        # Remove the pipes parallel to the valves to avoid redundancy
        for pipe in removed_pipes:
            wn.remove_link(pipe)
        return wn

            
    def save_wn(self, wn, filename):
        ''' 
        Save the network as inp file
        '''
        # Change the base_demand back to average demand instead of peak demand, and introduce demand paterns
        for node_name in wn.node_name_list[:-self._number_of_reservoirs]:
            node = wn.get_node(node_name)
            # to convert to m3/h
            node.demand_timeseries_list[0].base_value = 0.25* node.base_demand *3600
        for sect in self._valves_per_sect:
            for valve_name in self._valves_per_sect[sect]:
                valve = wn.get_link(valve_name)
                valve.status = 0
                valve.initial_status = wntr.network.base.LinkStatus(valve.status)
        # Set the time duration
        wn.options.time.duration = 86400
        wn.write_inpfile(filename, units='LPS')
        #wn.write_inpfile(filename, units='CMH')
        return wn


    def stats(self, wn):
        '''
        Calculates some stats of the network
        '''
        G = wn.get_graph()
        G = nx.Graph(G)
        degrees = nx.degree_histogram(G)
        # Number of nodes
        n_nodes = len(G)
        # Number of edges
        n_edges = G.size()
        # Total length of the network
        lengths = [pipe.length for pipe_name, pipe in wn.pipes()]
        # Pipe diameters
        pipes_diameters = [pipe.diameter for pipe_name, pipe in wn.pipes()]        
        return {'Number of nodes': n_nodes, 'Number of edges': n_edges, 'Number of reservoirs':self._number_of_reservoirs, 'Velocity' : self._velocity, 'Pressure' : self._pressure, 'Lengths': lengths, 'Diameters': pipes_diameters, 'Degrees': degrees, 'Main valves': self._main_valves, 'Valves per sector': self._valves_per_sect}


    ###################################################################################################################
    # PRIVATE METHODS
    ###################################################################################################################
    def proj_distance(self, y1, x1, y2, x2, z1=0, z2=0):
        '''
        Methode for calculating distance from projected coordinates
        '''
        a = (y1-y2)**2 + (x1-x2)**2 + (z1-z2)**2
        return(math.sqrt(a))


    def tooclose(self, G, node, tolerance):
        '''
        Finds the nodes that are within `tolerance' distance the node.
        '''
        cluster = []
        if G.has_node(node):
            x1 = G.nodes[node]['x']
            y1 = G.nodes[node]['y']
            for u in G.nodes():
                x2 = G.nodes[u]['x']
                y2 = G.nodes[u]['y']
                if self.proj_distance(y1, x1, y2, x2) <= tolerance:
                    cluster.append(u)             
            cluster.remove(node)
            cluster.append(node)
        return cluster


    def zone_position(self, G, zone):
        '''
        Finds zone position in relation to the geographic center of the graph
        '''
        center = self.graph_center(G)
        zone_center = zone['center']
        if self.proj_distance(center[1], center[0], zone_center[1], zone_center[0]) <= self._distance/2:
            return 'middle'
        elif (zone_center[0] > center[0] and abs(zone_center[1]-center[1]) <= self._distance/2):
            return 'middle east'
        elif (zone_center[0] < center[0] and abs(zone_center[1]-center[1]) <= self._distance/2):
            return 'middle west'
        elif (zone_center[1] < center[1] and abs(zone_center[0]-center[0]) <= self._distance/2):
            return 'middle south'
        elif (zone_center[1] > center[1] and abs(zone_center[0]-center[0]) <= self._distance/2):
            return 'middle north'
        elif (zone_center[0] > center[0] and abs(zone_center[1]-center[1]) > self._distance/2):
            return 'north east'
        elif (zone_center[0] < center[0] and abs(zone_center[1]-center[1]) > self._distance/2):
            return 'south west'
        elif (zone_center[1] < center[1] and abs(zone_center[0]-center[0]) > self._distance/2):
            return 'south east'
        elif (zone_center[1] > center[1] and abs(zone_center[0]-center[0]) > self._distance/2):
            return 'north east'
    
        
    def find_points(self, G, zone, sizing = True):
        '''
        Generates a list of points for the zone depending oh it's position relative to the center
        of the graph G
        Sizing = True if we are trying to figure out the sizes of the zones (not trying to find extremity points
        depending on the position of the zone)
        '''         
        zone_center = zone['center']
        dx = zone['dx']
        dy = zone['dy']
        points = []
        if sizing == True:
            for node in list(G.nodes):
                if (zone_center[0]-dx <= G.nodes[node]['x'] <= zone_center[0]+dx and zone_center[1]-dy <= G.nodes[node]['y'] <= zone_center[1]+dy): points.append(node)        
        else:
            position = self.zone_position(G, zone)
            if position == 'middle':
                for node in list(G.nodes):
                    if (zone_center[0]-dx <= G.nodes[node]['x'] <= zone_center[0]+dx and zone_center[1]-dy <= G.nodes[node]['y'] <= zone_center[1]+dy): points.append(node)
            elif position == 'middle east':
                for node in list(G.nodes):
                    if (zone_center[0] <= G.nodes[node]['x'] <= zone_center[0]+dx and zone_center[1]-dy <= G.nodes[node]['y'] <= zone_center[1]+dy): points.append(node)
            elif position == 'middle west':
                for node in list(G.nodes):
                    if (zone_center[0]-dx <= G.nodes[node]['x'] <= zone_center[0] and zone_center[1]-dy <= G.nodes[node]['y'] <= zone_center[1]+dy): points.append(node)
            elif position == 'middle south':
                for node in list(G.nodes):
                    if (zone_center[0]-dx <= G.nodes[node]['x'] <= zone_center[0]+dx and zone_center[1]-dy <= G.nodes[node]['y'] <= zone_center[1]): points.append(node)
            elif position == 'middle north':
                for node in list(G.nodes):
                    if (zone_center[0]-dx <= G.nodes[node]['x'] <= zone_center[0]+dx and zone_center[1] <= G.nodes[node]['y'] <= zone_center[1]+dy): points.append(node)
            elif position == 'north east':
                for node in list(G.nodes):
                    if (zone_center[0] <= G.nodes[node]['x'] <= zone_center[0]+dx and zone_center[1] <= G.nodes[node]['y'] <= zone_center[1]+dy): points.append(node)
            elif position == 'north west':
                for node in list(G.nodes):
                    if (zone_center[0]-dx <= G.nodes[node]['x'] <= zone_center[0] and zone_center[1] <= G.nodes[node]['y'] <= zone_center[1]+dy): points.append(node)
            elif position == 'south east':
                for node in list(G.nodes):
                    if (zone_center[0] <= G.nodes[node]['x'] <= zone_center[0]+dx and zone_center[1]-dy <= G.nodes[node]['y'] <= zone_center[1]): points.append(node)
            elif position == 'south west':
                for node in list(G.nodes):
                    if (zone_center[0]-dx <= G.nodes[node]['x'] <= zone_center[0] and zone_center[1]-dy <= G.nodes[node]['y'] <= zone_center[1]): points.append(node)
        return points


    def graph_center(self, projected_G):
        '''
        Methode for finding the coordinates of the center of the graph
        '''
        x, y = 0, 0
        for node in projected_G.nodes():
            x += projected_G.nodes[node]['x']
            y += projected_G.nodes[node]['y']
        return (x/len(projected_G), y/len(projected_G))


    def outsider_zones(self, zones):
        dx = self._distance
        dy = self._distance
        upper_left = self._center
        upper_middle = self._center
        upper_right = self._center
        middle_left = self._center
        middle_right = self._center
        down_left = self._center
        down_middle = self._center
        down_right = self._center
        out_zones={}
        for zone in zones:
            zone_center = zones[zone]['center']
            # the x middle section
            if self._center[0]-dx/2 <= zone_center[0] <= self._center[0]+dx/2:
                if zone_center[1] >= upper_middle[1]: upper_middle, out_zones['upper_middle'] = zone_center, zone
                if zone_center[1] <= down_middle[1]: down_middle, out_zones['down_middle'] = zone_center, zone
            # the y middle section
            if self._center[1]-dy/2 <= zone_center[1] <= self._center[1]+dy/2:
                if zone_center[0] <= middle_left[0]: middle_left, out_zones['middle_left'] = zone_center, zone
                if zone_center[0] >= middle_right[0]: middle_right, out_zones['middle_right'] = zone_center, zone
            # the upper section
            # the right
            if zone_center[0] >= upper_right[0] and zone_center[1] >= upper_right[1]:
                upper_right, out_zones['upper_right'] = zone_center, zone
            # the left
            if zone_center[0] <= upper_left[0] and zone_center[1] >= upper_left[1]:
                upper_left, out_zones['upper_left'] = zone_center, zone        
            # the down section
            # the right
            if zone_center[0] >= down_right[0] and zone_center[1] <= down_right[1]:
                down_right, out_zones['down_right'] = zone_center, zone
            # the left
            if zone_center[0] <= down_left[0] and zone_center[1] <= down_left[1]:
                down_left, out_zones['down_left'] = zone_center, zone                        
        return out_zones


    def zonage(self, G, max_nodes=random.randint(100, 120), Reservoirs=False):
        '''
        Divides the graph to zones have at most max_nodes
        '''
        self._center = self.graph_center(G)
        center = self.graph_center(G)
        dx = self._distance
        dy = self._distance
        npoints = len(G)
        zones = {'0': 0}
        axis = 'x'
        biggest = '0'
        z = 0
        if Reservoirs== False: 
            max_nodes = min([int(len(G)/8), max_nodes])
            # max_nodes = len(G)/5
        while (npoints >= max_nodes):
            # print('npoints {}' .format(npoints))
            del zones[biggest]
            if axis == 'x':
                center1 = (center[0] + dx/2, center[1]) 
                zone1 = {'center': center1, 'dx': dx/2, 'dy': dy, 'axis': 'y'}
                zone1['points'] = self.find_points(G, zone1)
                center2 = (center[0] - dx/2, center[1])
                zone2 = {'center': center2, 'dx': dx/2, 'dy': dy, 'axis': 'y'}
                zone2['points'] = self.find_points(G, zone2)
                zones['{}'.format(z)] = zone1
                zones['{}'.format(z+1)] = zone2
                z += 2
                axis = 'y'
            else:
                center1 = (center[0], center[1] + dy/2) 
                zone1 = {'center': center1, 'dx': dx, 'dy': dy/2, 'axis': 'x'}
                zone1['points'] = self.find_points(G, zone1)
                center2 = (center[0], center[1] - dy/2)
                zone2 = {'center': center2, 'dx': dx, 'dy': dy/2, 'axis': 'x'}
                zone2['points'] = self.find_points(G, zone2)
                zones['{}'.format(z)] = zone1
                zones['{}'.format(z+1)] = zone2
                z += 2
                axis = 'x'            
            current = zones[list(zones)[0]]
            npoints = len(current['points'])
            center = current['center']
            dx = current['dx']
            dy = current['dy']
            axis = current['axis']
            biggest = list(zones)[0]                        
            for zone in zones:
                if len(zones[zone]['points']) > npoints:
                    npoints = len(zones[zone]['points'])
                    center = zones[zone]['center']
                    dx = zones[zone]['dx']
                    dy = zones[zone]['dy']
                    axis = zones[zone]['axis']
                    biggest = zone
        return zones


    def main_connected(self, G):
        cur_graph = G # whatever graph you're working with
        if not nx.is_connected(cur_graph):
            # get a list of unconnected networks
            sub_graphs = list(nx.connected_components(cur_graph))
            main_graph = list(sub_graphs[0])
            # find the largest network in that list
            for sg in sub_graphs:
                if len(list(sg)) > len(main_graph):
                    main_graph = list(sg)
            return G.subgraph(main_graph).copy()    
        else: return G


    def direction(self, ux, uy, vx, vy):
        return {'ex': vx-ux, 'ey': vy-uy}


    def in_way(self, G, u, v, cluster):
        nodes_in_way = {}
        udata = G.nodes[u]
        vdata = G.nodes[v]
        ux, uy = udata['x'], udata['y']
        vx, vy = vdata['x'], vdata['y']
        dir_uv = self.direction(ux, uy, vx, vy)
        ortho_uv = {'ex': -dir_uv['ey'], 'ey': dir_uv['ex']}
        for n in cluster:
            if n!=u:
                ndata = G.nodes[n]
                nx, ny = ndata['x']-ux, ndata['y']-uy
                projn_uv = (nx*dir_uv['ex']+ny*dir_uv['ey'])/(math.sqrt(dir_uv['ex']**2+dir_uv['ey']**2))
                projn_ortho_uv = (nx*ortho_uv['ex']+ny*ortho_uv['ey'])/(math.sqrt(ortho_uv['ex']**2+ortho_uv['ey']**2)) 
                if (0 < projn_uv < self.proj_distance(uy, vy, ux, vx) and abs(projn_ortho_uv) < 50 and n != v):
                    nodes_in_way['{}'.format(n)] = projn_uv
        nodes_in_way = {k: v for k, v in sorted(nodes_in_way.items(), key=lambda item: item[1])}
        return nodes_in_way


    def connect_zones(self, G, zones):
        '''
        Connect the zones with what will be the the edges of the main distribution
        network
        '''
        # Connect the centers of the zones
        edges_list = []
        neighbors = {}
        # We find the neighbors of each zone
        for zone in zones:
            neighbors['{}'.format(zone)]=[]
            center = zones[zone]['center']
            axis = zones[zone]['axis']
            zones[zone]['points'] = self.find_points(G, zones[zone], False)
            if axis == 'x': radius = 2*zones[zone]['dx'] 
            else: radius = 2*zones[zone]['dy']
            while len(neighbors['{}'.format(zone)]) < 1:
                for n in zones:
                    if n!=zone:
                        n_center = zones[n]['center']
                        distance = self.proj_distance(center[1], center[0], n_center[1], n_center[0])
                        if distance <= radius+10: neighbors['{}'.format(zone)].append(n)
                radius += radius/2
        p_list = []
        for z in neighbors:
            indicator = True
            if len(zones['{}'.format(z)]['points']) > 0:
                for p in p_list:
                    if int(p[1]) == int(z):
                        p1 = p[0]
                        indicator = False
                if indicator:
                    p1 = random.choice(zones['{}'.format(z)]['points'])
                p_list.append((p1, z))
            else: continue
            for n in neighbors['{}'.format(z)]:
                indicator = True
                if len(zones['{}'.format(n)]['points']) > 0:
                    for p in p_list:
                        if int(p[1]) == int(n):
                            p2 = p[0]
                            indicator = False
                    if indicator:
                        p2 = random.choice(zones['{}'.format(n)]['points'])
                    p_list.append((p2, n))
                try:
                    path = nx.shortest_path(G, p1, p2)
                except Exception:
                    continue
                edges = [(path[j], path[j+1]) for j in range(len(path)-1)]
                edges_list.append(edges)

        out_zones = self.outsider_zones(zones)
        # print(out_zones)
        ul, ml, dl, dm, dr, mr, ur, um = None, None, None, None, None, None, None, None
        for zone in out_zones:
            z = out_zones[zone]
            coords = [(p, G.nodes[p]['x'], G.nodes[p]['y']) for p in zones[z]['points']]
            # print(coords)
            if len(coords) > 0:
                if zone == 'upper_left':
                    sort = sorted(coords, key=lambda x: x[2], reverse=True)[0:min(9, len(coords))]
                    ul = min(sort, key=lambda x: x[1])[0] 
                if zone == 'middle_left':
                    ml = min(coords, key=lambda x: x[1])[0]
                if zone == 'down_left':
                    sort = sorted(coords, key=lambda x: x[2])[0:min(9, len(coords))]
                    dl = min(sort, key=lambda x: x[1])[0]
                if zone == 'down_middle':
                    dm = min(coords, key=lambda x: x[2])[0]
                if zone == 'down_right':
                    sort = sorted(coords, key=lambda x: x[2])[0:min(9, len(coords))]
                    dr = max(sort, key=lambda x: x[1])[0]
                if zone == 'middle_right':
                    mr = max(coords, key=lambda x: x[1])[0]
                if zone == 'upper_right':
                    sort = sorted(coords, key=lambda x: x[2], reverse=True)[0:min(9, len(coords))]
                    ur = max(sort, key=lambda x: x[2])[0]
                if zone == 'upper_middle':
                    um = max(coords, key=lambda x: x[2])[0]

        self._out_points = {
            'upper_left': ul,
            'middle_left': ml,
            'down_left': dl,
            'down_middle': dm,
            'down_right': dr,
            'middle_right': mr,
            'upper_right': ur,
            'upper_middle': um
        }

        points1, points2 = self._out_points.copy(), self._out_points.copy()
        del points1['upper_middle']
        points2['upper_middle'] = um
        del points2['upper_left']
        points2['upper_left'] = ul
        for (p1,p2) in list(zip(points1, points2)):
            # if p1 != None and p2 != None:
            try:
                path = nx.shortest_path(G, out_points[p1], out_points[p2])
                edges = [(path[j], path[j+1]) for j in range(len(path)-1)]
                edges_list.append(edges)
            except Exception:
                continue            

        edges_list = [u for v in edges_list for u in v]
        return edges_list            


    def reservoir_zone(self, G, local_highest, head):
        reservoir = G.nodes[local_highest]
        sf = self.sf
        xr = reservoir['x']
        yr = reservoir['y']
        zr = reservoir['elevation']+head
        reservoir_zone = []
        for node in G.nodes():
            x = G.nodes[node]['x']
            y = G.nodes[node]['y']
            z = G.nodes[node]['elevation']
            l = self.proj_distance(y, x, yr, xr, z, zr)
            if l*sf+20 < zr-z:
                try:
                    l = nx.algorithms.shortest_path_length(G, local_highest, node, 'length')
                except Exception:
                    continue
            # Estimate pressure at node
            p = zr-z-sf*l
            if p > 20 :
                reservoir_zone.append(node)
            else :
                continue
        return reservoir_zone


    def highest_elev(self, G):
        highest_elevation = 0
        for node in G.nodes():
            z = G.nodes[node]['elevation']
            if z > highest_elevation:
                highest_elevation = z
                highest_node = node
        return highest_node

