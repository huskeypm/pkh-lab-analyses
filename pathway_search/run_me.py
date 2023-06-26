#!/usr/bin/env python
# coding: utf-8

import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
from functools import reduce
from tabulate import tabulate
from operator import itemgetter, attrgetter
import graphviz
import pygraphviz as pgv
import yaml
import os

##### Workflow:
##1. Load networks. In this case, 2 manually curated networks: networkExtended-looped4.txt and networkExtended_edit3.txt. Combine these 2 networks.
##2. Initialize the network to put data into nodes. Edges are 'interactions' that are part of the network and converted to 1 (activate) or -1 (inhibit).
##3. Load RNA data.  
##4. Load antibody (Ab) data, if available. (To do)
##5. Add edge weights. 1 for activating interactions and 100 for inhibiting interactions.
##6. Calculate weights based on nodes and edge weights. Node weights are calculated using rna data (fold change) and antibody (ab) data, if available.
##7. Specify the start and end for the pathway search using the weighted search - dijkstra_path. In this case, the start are the receptors in the ligandlist and end is M1_polarization.
##8. Sort by weight and write an output.


### Load data
# current directory path
#input_dir = "C:/Users/Geraldine Desktop/OneDrive - Loyola University Chicago/Documents/PKH22/trial"
input_dir = os.getcwd()

# load network data from file
def list_edges(file):
    sif_file = open(input_dir + "/" + file, 'r')
    line_list = []
    for line in sif_file:
        #print(line)
        line_item = line.strip().split('\t')
        line_list.append(line_item)
        #print(line_list)
    #print(line_list)
    a_line = np.array(line_list)
    #a_line
    a_line[0:,1] = np.where(a_line[0:,1]=='up-regulates',1,-1)
    #print(a_line)
    #create list of 3-tuple edges from network data 
    network = [(a_line[i,0], a_line[i,2], {'interaction':a_line[i,1]}) for i in range(len(line_list))]
    #print(network)
    return network

# convert to agraph and change color & shape of arrow
def agraph_color(network):
    # convert to pgv graph
    G = nx.nx_agraph.to_agraph(network)
    for u,v in G.edges():
        z = G.get_edge(u,v)
        if z.attr["interaction"] == '-1':
            z.attr["color"] = 'red'
            z.attr["arrowhead"] ="tee"
        else:
            z.attr["color"] = 'black'
    return G

def pos_neg(path):
    sg = newDG2.subgraph(path)
    z = [int(e[2]['interaction']) for e in sg.edges(data=True)]
    result = reduce((lambda a, b:a*b), z)
    inh = 'inhibit' if result==-1 else 'promote'
    return inh

def get_inh(path, network):
    inter_list = [int(network[path[i]][path[i+1]]['interaction']) for i in range(len(p)-1)]
    result = reduce((lambda a, b:a*b), inter_list)
    inh = 'inhibit' if result==-1 else 'promote'
    return inh

#load network
x = list_edges('networkExtended-looped4.txt')
#x

#Create directional graph using network  
X = nx.DiGraph(x)

list(X.edges(data=True))

#show graph
aX = agraph_color(X)
aX.layout('dot')
aX

#load network 
y = list_edges('networkExtended_edit3.txt')
#y

#Create directional graph using network  
Y = nx.DiGraph(y)

list(Y.edges(data=True))

#show graph
aY = agraph_color(Y)
aY.layout('dot')
aY
# make an image file
#G.draw('G.png')

len(aY)

# create pic of network
#aY.draw('aY_23.png')

X.number_of_nodes()
X.number_of_edges()

Y.number_of_nodes()
Y.number_of_edges()

# combine 2 networks
z = x + y
z

Z = nx.DiGraph(z)

Z.add_edges_from(Y.edges)

list(Z.edges(data=True))

Z.number_of_nodes()
Z.number_of_edges()

aZ = agraph_color(Z)
aZ.layout('dot')
aZ
# make an image file
#aZ.draw('aZ.png')

list(Z.edges(data=True))

# Create new network using y network above to put node data

newDG2 = nx.DiGraph(Z)

newDG2 = nx.DiGraph(Z)
nx.set_node_attributes(newDG2, 0, 'fc')
nx.set_node_attributes(newDG2, 0, 'pval')
nx.set_node_attributes(newDG2, 0, 'ab')
nx.set_node_attributes(newDG2, 0, 'pab')
list(newDG2.nodes(data=True))

# Load new RNA data M1 vs M2 fdr 10%

rna_m1m2fdr10 = pd.read_csv(input_dir + "/" + 'm1m2fdr10.csv')
rna_m1m2fdr10

rna_m1m2fdr10['logFC']

# make list of 3-tuple from rna data
rna_m1m2fdr10_list = [(rna_m1m2fdr10['Gene.symbol'][i], rna_m1m2fdr10['adj.P.Val'][i], rna_m1m2fdr10['logFC'][i])  for i in range(len(rna_m1m2fdr10))]
#rna_m1m2fdr10_list

#create an rna dictionary to put the gene name, adj. p-value and fold change value
rna_dict3 = {}
for rna_tuple in rna_m1m2fdr10_list:
    rna_dict3[rna_tuple[0]] = {'pval':rna_tuple[1], 'fc':rna_tuple[2]}
#rna_dict3

# put rna data into nodes in network
nx.set_node_attributes(newDG2, rna_dict3)
list(newDG2.nodes(data=True))

#create a ranked list

#1st, reverse order
rev_rna_list = list(reversed(rna_m1m2fdr10_list))
#rev_rna_list

scores = [item[1] for item in rev_rna_list]
#scores

ranks = []
for s in scores:
    ranks.append(scores.index(s)+1)
#ranks

len(scores)
len(ranks)

#create an rna dictionary to put the gene name and rank
gene_rank = {}
for i in range(len(rna_m1m2fdr10_list)):
    gene_rank[rna_m1m2fdr10_list[i][0]] = {'rank':ranks[i]}
gene_rank


rna_m1m2fdr10_list[0][0]


#initialize network before putting rank as node data
nx.set_node_attributes(newDG2, 0, 'rank')
list(newDG2.nodes(data=True))


# add rank in the network nodes
nx.set_node_attributes(newDG2, gene_rank)
list(newDG2.nodes(data=True))


newDG2_nlist = list(newDG2.nodes(data=True))
#newDG2_nlist


list(newDG2.edges(data=True))

Zedgelist = list(Z.edges)
Zedgelist

#create function to calculate edge weight using gene_rank

def calc_wt(edge, gene_rank):
    node1 = edge[0]
    node2 = edge[1]
    node1_wt = gene_rank[node1]['rank'] if node1 in gene_rank else 0
    node2_wt = gene_rank[node2]['rank'] if node2 in gene_rank else 0
    #print(node1_wt, node2_wt)
    wt = (node1_wt + node2_wt)/2
    return wt

#add activating interactions weight

def get_inh(path, network):
    inter_list = [int(network[path[i]][path[i+1]]['interaction']) for i in range(len(p)-1)]
    result = reduce((lambda a, b:a*b), inter_list)
    inh = 'inhibit' if result==-1 else 'promote'
    return inh

def calc_inter_wt(edge, network):
    edge_inter = nx.get_edge_attributes(network, 'interaction')
    #edge_inter[edge]
    if int(edge_inter[edge]) > 0:
        i_wt = 1
    else:
        i_wt = 100
    return i_wt

calc_inter_wt(('CCL5', 'CCR5'), newDG2)


#initialize network to add 'i_wt' interactions weight

iwt_attr = {}
for edge in Zedgelist:
    #print(edge)
    iwt_attr[edge] = {'i_wt':calc_inter_wt(edge, newDG2)}
    #wt_attr.update = {(edge): {'wt':0}}
#print(iwt_attr)

#add activating interactions weight to network
nx.set_edge_attributes(newDG2, iwt_attr)
list(newDG2.edges(data=True))

# Add node weights

#rna only
def node_wt_rna(u, v):
    rna_node1 = 1 if newDG2.nodes[u].get("fc", 0) > 0 else -1 if newDG2.nodes[u].get("fc", 0) < 0 else 0
    rna_node2 = 1 if newDG2.nodes[v].get("fc", 0) > 0 else -1 if newDG2.nodes[v].get("fc", 0) < 0 else 0
    node_u_wt = rna_node1
    node_v_wt = rna_node2
    edge_wt = newDG2.edges[u, v].get("i_wt", 1)
    return edge_wt - node_u_wt / 2 - node_v_wt / 2 


#rna+ab
def node_wt_rna_ab(u, v):
    rna_node1 = 1 if newDG2.nodes[u].get("fc", 0) > 0 else -1 if newDG2.nodes[u].get("fc", 0) < 0 else 0
    rna_node2 = 1 if newDG2.nodes[v].get("fc", 0) > 0 else -1 if newDG2.nodes[v].get("fc", 0) < 0 else 0
    ab_node1 = 1 if 1 > newDG2.nodes[u].get("ab", 0) > 0 else -1 if newDG2.nodes[u].get("ab", 0) > 1 else 0
    ab_node2 = 1 if 1 > newDG2.nodes[v].get("ab", 0) > 0 else -1 if newDG2.nodes[v].get("ab", 0) > 1 else 0
    node_u_wt = (rna_node1 + ab_node1) / 2
    node_v_wt = (rna_node2 + ab_node2) / 2
    edge_wt = newDG2.edges[u, v].get("i_wt", 1)
    return edge_wt - node_u_wt / 2 - node_v_wt / 2


node_wt_rna('CCL5', 'CCR5')

node_wt_rna_ab('CCL5', 'CCR5')


#make dictionary to add edge+node(rna only) weights
iwt_rna_attr = {}
for edge in Zedgelist:
    iwt_rna_attr[edge] = {'iwt_rna':node_wt_rna(edge[0], edge[1])}
#print(iwt_rna_attr)

#add edge+node(rna only) weights
nx.set_edge_attributes(newDG2, iwt_rna_attr)
list(newDG2.edges(data=True))

#make dictionary to add edge+node(rna+ab) weights
iwt_rna_ab_attr = {}
for edge in Zedgelist:
    iwt_rna_ab_attr[edge] = {'iwt_rna_ab':node_wt_rna_ab(edge[0], edge[1])}
#print(iwt_rna_ab_attr)


#add edge+node(rna +ab) weights
nx.set_edge_attributes(newDG2, iwt_rna_ab_attr)
list(newDG2.edges(data=True))


with open('ligandlist2.yaml', 'r') as file:
    ligandlist = yaml.safe_load(file)


#ligandlist

r_list = [r for r in ligandlist.keys()]
r_list

end = ['M1_polarization']

# paths from receptor to end (M1_polarization)
r_paths = []
for s in r_list:
    for t in end:
        #check if node is in the network (s in DG)
        #check if path exist in the network (nx.has_path(DG, s, t)
        #if s in network and t in network and nx.has_path(network, s, t)
        if s in newDG2 and t in newDG2 and nx.has_path(newDG2, s, t):
            r_paths.append(nx.dijkstra_path(newDG2, s, t, weight='iwt_rna_ab'))
print(r_paths)


# pathways of all the receptors in the ligandlist
r_path_data = []
for p in r_paths:
    p_i_wt = nx.path_weight(newDG2, p, 'iwt_rna_ab')
    l = ligandlist[p[0]]
    r_path_data.append([p[0], get_inh(p, newDG2), p[-1], p_i_wt, p, l])
head = ['receptor', 'action', 'end', 'weight', 'pathway', 'ligand']
print(tabulate(sorted(r_path_data, key=itemgetter(3)), headers=head))


### Output to Chris' paper SI Table 2

# get only target receptors and show only 'receptor', 'weight', 'pathway', 'ligand'
select_rec = ['TLR1','TLR2', 'TLR4', 'CCR1', 'CCR3', 'CCR5', 'TNFRSF1A', 'TNFRSF1B', 'ITGAL', 'ITGB2']


# select the end point of pathway
end = ['M1_polarization']


# get the pathways from receptor to end / outcome
paths = []
for s in select_rec:
    for t in end:
        #check if node is in the network (s in DG)
        #check if path exist in the network (nx.has_path(DG, s, t)
        #if s in network and t in network and nx.has_path(network, s, t)
        if s in newDG2 and t in newDG2 and nx.has_path(newDG2, s, t):
            paths.append(nx.dijkstra_path(newDG2, s, t, weight='iwt_rna_ab'))
print(paths)


# get the scores
path_data = []
for p in paths:
    p_i_wt = nx.path_weight(newDG2, p, 'iwt_rna_ab')
    l = ligandlist[p[0]]
    path_data.append([p[0], p_i_wt, p, l])
head = ['receptor', 'weight', 'pathway', 'ligand']
print(tabulate(sorted(path_data, key=itemgetter(1)), headers=head))

# create an output file
outfile = open(input_dir+"/"+"rec_lig_paths.txt", 'w')
#head = ['receptor', 'action', 'end', 'weight', 'pathway', 'ligand']
#print(tabulate(sorted(r_path_data, key=itemgetter(3)), headers=head))
outfile.write('receptor, weight, pathway, ligand\n')
sorted_paths = sorted(path_data, key=itemgetter(1))
for path in sorted_paths:
    #outfile.write("%s\n" % path)
    outfile.write("%s\n" % path)
#outfile.write('\t'.join(str(path)))
#close output file object
outfile.close()





