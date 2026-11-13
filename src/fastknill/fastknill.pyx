import networkx as nx
from scipy.special import binom
from itertools import combinations
from math import dist
import numpy as np
import re
from scipy.spatial.distance import cdist
from sklearn.neighbors import kneighbors_graph


def knn_graph(W, k):
    """Filter symmetric adjacency matrix to keep k strongest neighbors per node, return NetworkX graph."""
    W = W.copy()
    n = W.shape[0]
    
    # Set self-connections to -inf to avoid selecting them
    np.fill_diagonal(W, -np.inf)
    
    # For each row, find indices of top k entries
    rows = np.arange(n)[:, None]
    top_k_idx = np.argpartition(-W, kth=k, axis=1)[:, :k]
    
    # Create a new zero matrix
    W_k = np.zeros_like(W)

    # Populate top-k entries per row
    W_k[rows, top_k_idx] = W[rows, top_k_idx]

    # Make symmetric: max(W_k, W_k.T) to preserve strongest edges in either direction
    W_k = np.maximum(W_k, W_k.T)
    
    # Create a NetworkX graph from the filtered adjacency matrix
    G = nx.from_numpy_array(W_k)
    
    return G



# Global neighborhood dictionary
cdef dict Neigh

def parse_graphml_to_networkx(filename):
    file = open(filename)
    Graph = nx.Graph()
    Dict = {}
    counter = -1
    weight_key = None
    key_info = "None"
    quoted_strings = ""

    for line in file:
        if "<key id" in line:
            quoted_strings = re.findall(r'"(.*?)"', line)
        if "weight" in quoted_strings:
            weight_key = quoted_strings[0]
            key_info = "<data key=" + '"' + weight_key + '"'
        if "<node id" in line:
            quoted_strings = re.findall(r'"(.*?)"', line)
            node = quoted_strings[0]
            counter += 1
            Dict[node] = counter
        if "<edge id" in line:
            quoted_strings = re.findall(r'"(.*?)"', line)
            Id, x, y = quoted_strings
            if x != y:
                Graph.add_edge(*[Dict[x], Dict[y]])
        if key_info in line:
            quoted_strings = re.findall(r'>(.*?)<', line)
            if x != y:
                Graph[Dict[x]][Dict[y]]["weight"] = float(quoted_strings[0])

    Graph.remove_edges_from(nx.selfloop_edges(Graph))
    Dict = {Dict[i]: i for i in Dict.keys()}
    return Graph, Dict    

def filter_graph(Obj, cutoff=np.inf,n_neighbors=None,metric="ëuclidean"):
    if isinstance(Obj, nx.Graph):
        L = list(Obj.nodes())
        mapping = {i: L[i] for i in range(len(L))}
        M = nx.to_numpy_array(Obj, weight="weight")
        M = np.where(M <= cutoff, M, 0)
        H = nx.from_numpy_array(M)
        if n_neighbors!=None:
            H=knn_graph(nx.to_numpy_array(H),k=n_neighbors)
        return H, mapping
    elif isinstance(Obj, dict):
        L = list(Obj.keys())
        mapping = {i: L[i] for i in range(len(L))}
        points = np.array(list(Obj.values()))
        M = cdist(points, points,metric)
        M = np.where(M <= cutoff, M, 0)
        H = nx.from_numpy_array(M)
        if n_neighbors!=None:
            H=knn_graph(nx.to_numpy_array(H),k=n_neighbors)
        return H, mapping
    elif isinstance(Obj, np.ndarray):
        mapping = {i: i for i in range(len(Obj))}
        M = np.where(Obj <= cutoff, Obj, 0)
        H = nx.from_numpy_array(M)
        if n_neighbors!=None:
            H=knn_graph(nx.to_numpy_array(H),k=n_neighbors)
        return H, mapping
    elif isinstance(Obj, str):
        H, mapping = parse_graphml_to_networkx(Obj)
        H = nx.to_numpy_array(H)
        H = np.where(H <= cutoff, H, 0)
        H = nx.from_numpy_array(H)
        if n_neighbors!=None:
            H=knn_graph(nx.to_numpy_array(H),k=n_neighbors)
        return H, mapping
    else:
        raise ValueError(f"Object type ({type(Obj)}) is not allowed to provide a network input: Object types allowed: nx.Graph, dict, np.array and graphml (str)")

# Manual set intersection counter for speed
cdef int intersect_size(set a, set b):
    cdef int count = 0
    for x in a:
        if x in b:
            count += 1
    return count



def compute_euler(object D, float cutoff=np.inf, n_neighbors=None, int max_dim=1, mapped_nodes=False,metric="euclidean"):
    """Computes the Euler Characteristics up to clique dimension (dim) from the object provided.
    
     Input: D, an object that can be:
    
        - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
        - A nx.Graph object. The edge's attibutes `weight` will be considered.
        - A np.ndarray object that represents a symmetric matrix of float numbers
        cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
        (if D is a nx.Graph) of float values (if D is a np.matrix);
    -n_neighbors: Interger, the number of the nearest neighbours to consider. If None, it is not applied.

    - A string, the file name of the graphml/xml file.
    - mapped_nodes: True or False - if the mapping of the original node labels are kept
    - metric: string, the metric in which the pairwise distance is computed. If D is not a dictionary, this parameter is ignored.
    
    dim: integer, the maximum simplex dimension allowed to the computation.
    
    Output:
    An integer, the Euler characteristics of the network"""
    global Neigh
    G, mapping = filter_graph(D, cutoff, n_neighbors,metric)

    cdef int i, j, n, d, x, y, old_node, w
    cdef int num_nodes = G.number_of_nodes()
    cdef tuple old_cell,edges
    cdef dict N = {d: {n: set() for n in G.nodes()} for d in range(1, max_dim + 1)}
    
    # Use a mutable container for E
    cdef list E = [G.number_of_nodes()]
   # cdef dict E={0:G.number_of_nodes()}

    def func(tuple cell):
        d = len(cell)
        E[0] += (-1) ** (d-1)
        #print(E[0])
        if d < max_dim + 1:
            if d == 2:
                x, y = cell
                N[d][x].update({y})
                N[d][y].update({x})
            else:
                old_cell = cell[:-1]
                old_node = cell[-1]
                for x in old_cell:
                    N[d][x].update({old_node})

            W = set.intersection(*[N[d][node] for node in cell])
            for w in W:
                new_cell = (*cell, w)
                func(new_cell)

    # Start the recursion with edges (1-simplices)
    for edge in G.edges():
        func(edge)

    return E[0]

def compute_knill(object D, float cutoff=np.inf, n_neighbors=None, int max_dim=1, mapped_nodes=False,metric="euclidean"):
    """Computes the Knill curvature of each node up to dimension (dim) from the object provided.
    
     Input: D, an object that can be:
    
        - A dictionary in which the keys are enumerated from 0 to len(D)-1 and the values are the N-dimensional points. 
        - A nx.Graph object. The edge's attibutes `weight` will be considered.
        - A np.ndarray object that represents a symmetric matrix of float numbers
        cutoff: A float number for the threshold values for distance (if D is a dictionary), weights 
        (if D is a nx.Graph) of float values (if D is a np.matrix);
    -n_neighbors: Interger, the number of the nearest neighbours to consider. If None, it is not applied.
    - A string, the file name of the graphml/xml file.
    - mapped_nodes: True or False - if the mapping of the original node labels are kept
    
    - max_dim: integer, the maximum simplex dimension allowed to the computation.
    - metric: string, the metric in which the pairwise distance is computed. If D is not a dictionary, this parameter is ignored.
    
    
    Output:
    a dictionary whose keys are the noded and the values are the Knill curvature of each node."""
    global Neigh
    G, mapping = filter_graph(D, cutoff, n_neighbors,metric=metric)

    cdef int i, j, n, d, x, y, old_node, w,node
    cdef int num_nodes = G.number_of_nodes()
    cdef tuple old_cell,edges
    cdef dict N = {d: {n: set() for n in G.nodes()} for d in range(1, max_dim + 1)}
    
    # Use a mutable container for E
    cdef dict K = {i:1. for i in range(G.number_of_nodes())}
   # cdef dict E={0:G.number_of_nodes()}

    def func(tuple cell):
        d = len(cell)
        for node in cell:
            K[node] += ((-1) ** (d-1))/(d)
        #print(E[0])
        if d < max_dim + 1:
            if d == 2:
                x, y = cell
                N[d][x].update({y})
                N[d][y].update({x})
            else:
                old_cell = cell[:-1]
                old_node = cell[-1]
                for x in old_cell:
                    N[d][x].update({old_node})

            W = set.intersection(*[N[d][node] for node in cell])
            for w in W:
                new_cell = (*cell, w)
                func(new_cell)

    # Start the recursion with edges (1-simplices)
    for edge in G.edges():
        func(edge)

    return K
