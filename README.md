# FastKnill - An alternative fast computation of Euler charactetistics and Knill curvature from networks

***
   **Author:** Danillo Barros de Souza
   -- 
   **ORCID: [0000-0002-7792-8862](https://orcid.org/0000-0002-7762-8862)**

***
This work is inspired by the proposed new methods for computations of cliques in Vietoris-Rips complexes [[1]](https://arxiv.org/abs/2502.14593). We provide a set-theoretical Python coding for efficiently computing the Knill curvature and the Euler characteristics for simplicial complexes. For alternative geometric computations, see also [[2]](https://github.com/danillodbs16/fastforman)

## Current version:
  0.1.0


## Content:

- `fastknill.pyx`
- `compiler.py`
- `setup.py`

## Python version: 
    3.8.5
## Package requirement:

- `numpy`
- `networkx`
- `cython`
- `scikit-learn`


## Installation:

The files `fastknill.pyx`, `compiler.py` and `setup.py` must be placed at the same directory.

 ### Local installation:
To this purpose, the user will need to compile the file `compiler.py` in cython by running the command:

```
python3 compiler.py build_ext --inplace
```

### Global installation:

Run `setup.py` by executing the command:

```
pip install .
```

After successful compilation, the generated `fastknill` file can be imported as a package, for instance:

```python
import fastforman as fk
```
## Functions:

- ``compute_knill``:

  ***Input:*** 
        
    - `D`, can be:
     - A `dictionary` in which the keys are integer numbers from 0 to len(D)-1 and the values are the N-dimensional points;

     - A symmetric `numpy.matrix` of float numbers;

     - A simple undirected `nx.Graph` that may include the feature `weights`on edges attributes;
     
     - A `string` of the `graph/graphml` file in a directory.

    - `cutoff`: Float number for the threshold distance between points.
    - `n_neighbors`: Interger, the number of the nearest neighbours to consider. If None, it is not applied.
    - `dim`: Integer, the maximum simplex dimension to compute FRC.
    - `mapped_nodes`: True or False - if the mapping of the original node labels is kept.
   
    
   ***Output:*** 
     
   A dictionary whose keys are the noded and the values are the Knill curvature of each node.
    
- ``compute_euler``:

  ***Input:*** 
        
    - `D`, can be:
     - A `dictionary` in which the keys are integer numbers from 0 to len(D)-1 and the values are the N-dimensional points;

     - A symmetric `numpy.matrix` of float numbers;

     - A simple undirected `nx.Graph` that may include the feature `weights` on edges attributes;
     
     - A `string` of the `graph/graphml` file in a directory.

    - `cutoff`: Float number for the threshold distance between points.
    - `n_neighbors`: Integer, the number of the nearest neighbours to consider. If None, it is not applied.
    - `dim`: Integer, the maximum simplex dimension to compute FRC.
    - `mapped_nodes`: True or False - if the mapping of the original node labels is kept.
   
    
   ***Output:*** 
     
   An integer, the Euler characteristics of the input network.
    

## Contact and Support:

danillo.dbs16@gmail.com, dbarros@bcamath.org

## References: 

[1] Alternative set-theoretical algorithms for efficient computations of cliques in Vietoris-Rips complexes, Barros de Souza, Danillo; da Cunha, Jontatas, A.N. Santos, Fernando; Desroches, Mathieu; Jost, Juergen & Rodigues, Serafim; [link](https://arxiv.org/abs/2502.14593)
[2] Efficient set-theoretic algorithms for computing high-order Forman-Ricci curvature on abstract simplicial complexes; Barros de Souza, Danillo; da Cunha, Jontatas, A.N. Santos, Fernando; Jost, Juergen & Rodigues, Serafim; [Link](https://royalsocietypublishing.org/doi/10.1098/rspa.2024.0364)

