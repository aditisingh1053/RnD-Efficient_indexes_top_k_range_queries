# Constrained Density Top-k (CDT) Problem

The **Constrained Density Top-k (CDT)** problem is defined as follows:

Given a set \(P\) of weighted points in \(\mathbb{R}^d\), a parameter \(\delta > 0\), and an integer \(k\), the goal is to select a subset \(S \subseteq P\) of at most \(k\) points such that:

- All pairs of points in \(S\) are at least \(\delta\) apart (i.e., \(S\) is \(\delta\)-diverse),
- The total weight of the selected points is maximized,
- Optionally, queries may specify a range from which the output points should be chosen.

This problem is important in scenarios where you want a diverse selection of the top-k most relevant (or heaviest) results.

---

## File Descriptions

### 🔹 `CDT_basic.cpp`
Basic greedy implementation of the CDT problem as described in the first paragraph of the paper.

- Picks the heaviest remaining point,
- Removes all points within distance \(\delta\) of the chosen point,
- Continues until \(k\) points are selected or no candidates remain.

### 🔹 `CDT_bbd.cpp`
An optimized version of the basic CDT using spatial data structures:

- **BBD-tree**: For efficient spatial partitioning and pruning.
- **Range trees**: To quickly find the heaviest point within a query range.

This version significantly improves performance over large datasets and higher dimensions.

### 🔹 `CDT_index1d.cpp`
Implements a \((1 - \varepsilon)\)-approximation algorithm for **1D data**:

- Uses **shifted grid technique** to discretize the space,
- Builds dynamic programming tables to store optimal sub-solutions,
- Combines them efficiently to produce a near-optimal, \(\delta\)-diverse top-k result.

### 🔹 `CDT_index2d.cpp`
Extends the index-based solution to **2D data**:

- Constructs 2D shifted grids \(G_0, ..., G_{s-1}\),
- Uses 2D range trees with secondary 1D structures,
- Partitions the query into center and boundary regions and processes them via dynamic programming.

Achieves fast and accurate CDT query processing in two dimensions.

### 🔹 `CDT_index_d.cpp`
Generalized index-based CDT solution for **\(d\)-dimensional data**:

- Builds a \(d\)-dimensional shifted grid,
- Constructs a \(d\)-dimensional range tree with secondary \((d-1)\)-dimensional trees,
- Handles both center and boundary region solutions,
- Merges them using efficient dynamic programming.

Offers high-quality approximations with polynomial preprocessing and logarithmic query time in higher dimensions.

---

## 🛠️ Compilation

Each file is self-contained and can be compiled using a standard C++11 or later compiler:

```bash
g++ -std=c++11 CDT_basic.cpp -o CDT_basic
g++ -std=c++11 CDT_bbd.cpp -o CDT_bbd
g++ -std=c++11 CDT_index1d.cpp -o CDT_index1d
g++ -std=c++11 CDT_index2d.cpp -o CDT_index2d
g++ -std=c++11 CDT_index_d.cpp -o CDT_index_d
