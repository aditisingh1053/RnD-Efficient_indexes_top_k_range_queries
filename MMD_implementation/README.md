## Max-Min Diversity (MMD) Problem

The **Max-Min Diversity (MMD)** problem is defined as follows:

Given a set \( P \) of weighted points in \( \mathbb{R}^d \), an integer \( k \), and a trade-off parameter \( \lambda \geq 0 \), the goal is to select a subset \( S \subseteq P \) of size exactly \( k \) such that the **minimum pairwise distance** among selected points and their **total weight** are jointly maximized.

### Objective Function

For a subset \( S \subseteq P \), define the MMD objective as:

\[
f(S) = \min_{\substack{p_i, p_j \in S \\ i < j}} \|p_i - p_j\| + \lambda \cdot \min_{p_i \in S} w(p_i)
\]

Where:
- \( \min_{\substack{p_i, p_j \in S}} \|p_i - p_j\| \) captures **diversity** by maximizing the **minimum distance** between any pair of points in \( S \),
- \( \min_{p_i \in S} w(p_i) \) captures the **minimum weight** among the selected points, promoting relevance across all elements,
- \( \lambda \geq 0 \) is a parameter that balances the two components.

The objective is to compute:

\[
\text{MMD}(P, k) = \arg\max_{\substack{S \subseteq P \\ |S| = k}} f(S)
\]

This ensures that all selected points are **well-separated** and each is **individually relevant**, making it suitable for applications needing robust and diverse summarization.

### Use Case: Implicit Range Queries

The MMD problem is especially useful when the point set \( P \) is too large to store or must be accessed implicitly through bounding regions or streaming queries. In such cases, MMD is evaluated over a dynamically generated subset \( X \subseteq P \) defined by the query, i.e.,

\[
\text{MMD}(X, k) = \arg\max_{\substack{S \subseteq X \\ |S| = k}} f(S)
\]

This variant is applicable in spatial databases, clustering with fairness constraints, and anytime selection tasks where diversity and representativeness both matter.


---

## File Descriptions

### 🔹 `MMD_explicit.cpp`
Implementation of the MMD problem as described in the first paragraph of the paper.
Input structure is directly given as stream of n offline points

### 🔹 `MMD_implicit.cpp`
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
Bash script `evaluate_MMD.sh` reads in inputs from `Testcases_exp` directory and compare the MSD objective function of the subset obtained from `MMD_explicit.cpp` with the correct outputs and gives the score `min(1,objective(MMD_explicit)/objective(correct))`. The intention behind using this script and scoring is that since the algorithm is approximation, we cant directly compare in case of large input sizes and we wanted to know the accuracy of the algorithm. Currently this only supports input types for explicit implementations.

---
## Usage of script
The script automatically runs the inputs.If you want to change the directory of Testcases, you can modify the scipt. It takes the cpp file name as argument.
./evaluate_MMD.sh <filename.cpp>

## Compilation

Each file is self-contained and can be compiled using a standard C++11 or later compiler:

```bash
g++ -std=c++11 MMD_explicit.cpp -o MMD_exp
g++ -std=c++11 MMD_implicit.cpp -o MMD_imp
