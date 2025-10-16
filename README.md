Name: Simon Ghyselincks
Number: 12145124

## Running Via GUI

The mesh modification controller is loaded with a mesh connectivity:

```c++
  mohecore::Mesh_connectivity mesh;
  mohecore::Mesh_io io(mesh);
  mohecore::Mesh_modifier_edge_collapse modi(mesh);
  modi.initialize();
```  

Initialize computes all the required $Q$ matrices for the collapse and sets up the priority queue.

In practice the main functions that are used are either popping a valid edge collapse from the heap, or peeking at the top n candidates without popping them.
```c++
int n = 10; // number of candidates to evaluate
std::vector<int> top_k_he = modi.get_top_n_candidates(n); // get the top n half edges integer ids

mohecore::Mesh_modifier_edge_collapse::MergeCandidate candidate{0.0, {0, 0}, Eigen::Vector3d::Zero(), 0};

// Perform edge collapses up to the requested number
int collapses = 0;
while (collapses < n && modi.get_min_pair(candidate))
{
    if (globalvars::modi.collapse_edge(candidate))
    {
        printf("Collapsing edge (%d, %d) with error %.9f\n",
                candidate.pair.v1, candidate.pair.v2, candidate.error);
        collapses++;
    }
}
```

The merge candidate is a struct that holds all the tracking information for a pair and sorts it by error in the heap. It will return a boolean if the collapse was successful. The get_min_pair function will pop the valid pair from the heap or return false if the heap is empty.

#### Data Structures
`VertexPair` gives a cannonical way to specify an undirected edge
`MergeCandidate` Stores all the needed data in the heap for the merge with expiration versioning
`SymQuadric` is a 4x4 symmetric matrix with accumulation and error computation methods

#### Full API
- *initializaiton* Forms all the Q matrices attached to $v$ and starts the heap
- `add_or_update_pair` is a method to add or update a vertex pair as valid, the new x_opt and error will auto compute and go into the heap
- `get_min_pair` pops the valid pair from the heap or returns false once it is empty, also cleans out stale entries
- `invalidate_pair` will invalidate the version of the pair that is in the heap
- `get_all_pairs_from_vertex` gives the set of all edges around a vertex, these can be used to invalidate or update vertices
- `get_top_n_candidates` returns the half edges that have the lowest error, also cleans the heap/queue, no pop
- `get_all_neighbors_from_vetex` collects neighbors into a set to do topology evaluation
- `is_legal_collapse` is a catch all function to evaluate legality of a collapse
- `collapse_edge` does the deed itself.