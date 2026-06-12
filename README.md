# Labyrinth Solver (A* Pathfinding with Custom Memory Allocator)

A high-performance C++ maze solver that finds the shortest path through a labyrinth using the **A*** search algorithm combined with a custom hash table and fixed-size memory allocator.

The project demonstrates several advanced systems-programming techniques:

* A* heuristic search
* Custom object allocation
* Hash-based vertex storage
* Priority queue optimization
* Path reconstruction
* Low-overhead memory management

---

## Overview

The application reads a text-based labyrinth from a file and searches for the shortest route from the entrance at the top of the maze to the exit at the bottom.

After a solution is found, the path is marked with `+` characters and the completed maze is printed to standard output.

Example:

```text
######## #######
#++++++#       #
# ####+# ##### #
# #   +# #   # #
# # ### # # # # #
# #     # # # # #
# ####### # # # #
#         #   # #
############### #
```

---

## Algorithm

The solver uses the **A*** pathfinding algorithm.

### Heuristic

The search heuristic is the Manhattan distance:

```text
h(x, y) = |x - goal_x| + |y - goal_y|
```

This heuristic is admissible for grid-based movement in four directions and guarantees an optimal solution.

### Cost Function

The algorithm evaluates vertices using:

```text
f(n) = g(n) + h(n)
```

Where:

* `g(n)` = distance traveled from the start
* `h(n)` = Manhattan distance to the exit

Vertices with the lowest estimated total cost are explored first.

---

## Architecture

### Vertex

Represents a maze cell explored by the algorithm.

```cpp
struct Vertex
{
    VertexKey key;
    Vertex* next;
    ...
};
```

Stores:

* Cell coordinates
* Path length
* Parent pointer for reconstruction
* Hash-chain linkage

---

### VertexKey

Represents maze coordinates.

```cpp
struct VertexKey
{
    int x;
    int y;
};
```

Provides:

* Coordinate storage
* Manhattan distance calculation
* Hash generation

---

### Hash Table

The solver stores explored cells in a fixed-size hash table:

```cpp
HASH_BITS = 12
HASH_SIZE = 4096
```

Collisions are resolved using linked lists.

Benefits:

* Constant-time average lookup
* Fast duplicate detection
* Reduced search overhead

---

### Priority Queue

The open set is implemented using:

```cpp
std::priority_queue<Priority>
```

Vertices are ordered by their estimated total path cost.

This allows the solver to always expand the most promising candidate first.

---

## Custom Memory Allocator

A major feature of this project is the custom fixed-size allocator.

### FixedAlloc

Inspired by MFC's `CFixedAlloc`.

Features:

* Block-based allocation
* Free-list management
* Constant-time allocation
* Constant-time deallocation
* Reduced heap fragmentation

```cpp
FixedAlloc Vertex::s_alloc(sizeof(Vertex), 64);
```

Instead of allocating every vertex individually through the system heap, vertices are obtained from preallocated blocks.

---

### Plex

Inspired by MFC's `CPlex`.

Used to:

* Allocate large memory blocks
* Chain allocated blocks together
* Release all memory at once

This significantly reduces allocator overhead during large searches.

---

## Search Process

1. Load maze from file.
2. Locate entrance and exit.
3. Insert initial vertex into the priority queue.
4. Repeatedly:

   * Remove best candidate.
   * Expand neighbors.
   * Evaluate heuristic.
   * Store/update vertices.
5. Stop when exit is reached.
6. Reconstruct path using parent pointers.
7. Print solved maze.

---

## Input Format

The labyrinth is provided as a text file.

Example:

```text
####### #######
#             #
# ##### ##### #
# #   #     # #
# # # ##### # #
# # #     # # #
# # ##### # # #
#       #     #
############# #
```

Rules:

* `#` = wall
* `' '` (space) = walkable cell
* Opening in first row = entrance
* Opening in last row = exit

---

## Building

### Visual Studio

```bash
cl labyrinth.cpp
```

### GCC / Clang

```bash
g++ -O2 -std=c++11 labyrinth.cpp -o labyrinth
```

---

## Usage

```bash
labyrinth maze.txt
```

Output:

```bash
++++++++++++
+##########+
+          +
++++++++++++
```

---

## Performance Characteristics

### Time Complexity

Average-case A* performance:

```text
O(E log V)
```

where:

* `V` = number of reachable cells
* `E` = number of connections

### Memory Complexity

```text
O(V)
```

The custom allocator minimizes heap allocation overhead and improves cache locality.

---

## Notable Techniques

* A* pathfinding
* Manhattan heuristic
* Custom hash table
* Separate chaining collision resolution
* Fixed-size memory pools
* Object free lists
* Path reconstruction
* Priority queue search optimization

---

## Educational Value

This project serves as a practical example of:

* Graph search algorithms
* Heuristic optimization
* Memory pool implementation
* Game AI pathfinding
* Systems-level C++ programming
* Performance-oriented data structures

It is particularly useful for studying how custom allocators and hash-based storage can significantly improve the efficiency of graph traversal algorithms.
