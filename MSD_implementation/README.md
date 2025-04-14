## Max-Sum Diversification (MSD) Problem

The **Max-Sum Diversification (MSD)** problem is defined as follows:

Given a set \( P \) of weighted points in \( \mathbb{R}^d \), an integer \( k \), and a trade-off parameter \( \lambda \geq 0 \), the goal is to select a subset \( S \subseteq P \) of size exactly \( k \) such that the **total diversity and relevance** of the points is maximized.

### Objective Function

For a subset \( S \subseteq P \), the MSD objective function is defined as:

\[
f(S) = \sum_{\substack{p_i, p_j \in S \\ i < j}} \|p_i - p_j\| + \lambda \sum_{p_i \in S} w(p_i)
\]

Where:
- \( \|p_i - p_j\| \) is the Euclidean distance between points \( p_i \) and \( p_j \), capturing **diversity**.
- \( w(p_i) \) is the **weight** (or relevance) of point \( p_i \).
- \( \lambda \) balances **diversity** vs **total weight**.

The goal is to find:

\[
\text{MS}(P, k) = \arg\max_{\substack{S \subseteq P \\ |S| = k}} f(S)
\]

This formulation selects the most **diverse** and **relevant** subset of exactly \( k \) points from \( P \).

### Use Case: Range Queries

In practical scenarios, we may be interested in computing MSD over a **range-constrained subset** \( X \subseteq P \), i.e., find:

\[
\text{MS}(X, k) = \arg\max_{\substack{S \subseteq X \\ |S| = k}} f(S)
\]

This is useful for spatial databases, recommender systems, or visual summarization, where relevance and diversity both matter and selections need to be localized or filtered.


---

## File Descriptions

### 🔹 `MSD_explicit.cpp`
Implementation of the MSD problem as described in the first paragraph of the paper.
Input structure is directly given as stream of n offline points

### 🔹 `MSD_implicit.cpp`
Implementation of same problem but with implicit representation of points in the form of boxes.

---

## Testcases
The directory structure is:
Testcases_exp
|
|--input
    |--input1
    |--input2
|-output
    |--output1
    |--output2
Similarly goes on for Testcases_imp
Note that the output for small testcases have been generated using brute force while that of large files using some approximation heurisics.
Input structure of Explicit is
|--number of points(n)    dimension(d)   size of subset(k)  lambda(l)
|--n lines each of length (d+1) where last input is weight of the point

Input structure of Implicit is
|--number of boxes(s)    dimension(d)   size of subset(k)  lambda(l)
|--s lines each containing number of points in ith box(n)
    |--n lines each of length (d+1) where last input is weight of the point

Output structure of both implicit and explicit is:
|--k lines each of length (d+1) where last input is weight of the point

---

## Evaluation
Bash script `evaluate_MSD.sh` reads in inputs from `Testcases_exp` directory and compare the MSD objective function of the subset obtained from `MSD_explicit.cpp` with the correct outputs and gives the score `min(1,objective(MSD_explicit)/objective(correct))`. The intention behind using this script and scoring is that since the algorithm is approximation, we cant directly compare in case of large input sizes and we wanted to know the accuracy of the algorithm. Currently this only supports input types for explicit implementations.

---
## Usage of script
The script automatically runs the inputs.If you want to change the directory of Testcases, you can modify the scipt. It takes the cpp file name as argument.
./evaluate_MSD.sh <filename.cpp>

## Compilation

Each file is self-contained and can be compiled using a standard C++11 or later compiler:

```bash
g++ -std=c++11 MSD_explicit.cpp -o MSD_exp
g++ -std=c++11 MSD_implicit.cpp -o MSD_imp
