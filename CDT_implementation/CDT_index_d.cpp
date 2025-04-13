#include <bits/stdc++.h>
using namespace std;

// -------------------- Data Structures --------------------

struct PointD {
    double weight;
    vector<double> coords; // d-dimensional coordinates
};

bool compareByWeightD(const PointD &a, const PointD &b) {
    return a.weight > b.weight;
}

void printPointD(const PointD &p) {
    cout << "(";
    for (int i = 0; i < p.coords.size(); i++) {
        cout << p.coords[i] << (i == p.coords.size()-1 ? "" : ", ");
    }
    cout << ") -> " << p.weight;
}

// DPEntry stores a candidate subset (solution) along with its total weight.
struct DPEntry {
    vector<PointD> subset;
    double totalWeight;
};

const double NEG_INF = -1e9;

// -------------------- Brute-force for a Single Cell (General d) --------------------
// Enumerate all subsets of cellPoints of size up to min(j_max, |cellPoints|).
// Only accept those subsets that are 1-diverse, i.e., any two points p,q satisfy
// EuclideanDistance(p,q) >= 1.
// Returns a DP table (vector of DPEntry) of size (j_max+1), where dp[j] is the best
// subset of exactly j points.
double euclideanDistance(const PointD &a, const PointD &b) {
    double sum = 0;
    for (int i = 0; i < a.coords.size(); i++) {
        double diff = a.coords[i] - b.coords[i];
        sum += diff*diff;
    }
    return sqrt(sum);
}

vector<DPEntry> computeCellDP_d(const vector<PointD>& cellPoints, int j_max) {
    int n = cellPoints.size();
    vector<DPEntry> best(j_max+1, { {}, NEG_INF });
    best[0] = { {}, 0.0 };

    int total = 1 << n;
    for (int mask = 1; mask < total; mask++){
        vector<PointD> subset;
        double sumW = 0;
        int count = 0;
        for (int i = 0; i < n; i++){
            if (mask & (1 << i)) {
                subset.push_back(cellPoints[i]);
                sumW += cellPoints[i].weight;
                count++;
            }
        }
        if (count > j_max) continue;
        bool valid = true;
        for (int i = 0; i < subset.size() && valid; i++){
            for (int j = i+1; j < subset.size(); j++){
                if (euclideanDistance(subset[i], subset[j]) < 1.0) {
                    valid = false;
                    break;
                }
            }
        }
        if (!valid) continue;
        if (sumW > best[count].totalWeight) {
            best[count] = { subset, sumW };
        }
    }
    return best;
}

// -------------------- Dynamic Programming Over Cells (General d) --------------------
// We assume we have ordered the grid cells lexicographically.
// cells: vector of cell point-sets (each cell is a vector<PointD>).
// cellDP: corresponding DP table for each cell, computed by computeCellDP_d.
// We combine them in order (like knapsack DP) to obtain a DP table for the union.
vector<DPEntry> combineCells_d(const vector<vector<PointD>>& cells, const vector<vector<DPEntry>>& cellDP, int k) {
    int m = cells.size();
    vector<vector<DPEntry>> dp(m+1, vector<DPEntry>(k+1, { {}, NEG_INF }));
    dp[0][0] = { {}, 0.0 };
    for (int i = 1; i <= m; i++){
        for (int j = 0; j <= k; j++){
            dp[i][j] = dp[i-1][j]; // Option: take 0 points from cell i.
            int maxPtsInCell = cellDP[i-1].size() - 1; // possible sizes in this cell.
            for (int p = 1; p <= min(j, maxPtsInCell); p++){
                if (dp[i-1][j-p].totalWeight == NEG_INF || cellDP[i-1][p].totalWeight == NEG_INF)
                    continue;
                double candidateWeight = dp[i-1][j-p].totalWeight + cellDP[i-1][p].totalWeight;
                if (candidateWeight > dp[i][j].totalWeight){
                    vector<PointD> candidateSubset = dp[i-1][j-p].subset;
                    candidateSubset.insert(candidateSubset.end(), cellDP[i-1][p].subset.begin(), cellDP[i-1][p].subset.end());
                    dp[i][j] = { candidateSubset, candidateWeight };
                }
            }
        }
    }
    cout << "DP over cells complete. Final DP values:" << endl;
    for (int j = 0; j <= k; j++){
        cout << "  For j = " << j << ", best weight = " << dp[m][j].totalWeight << endl;
    }
    return dp[m];
}

// -------------------- Generalized d-Dimensional Grid Construction --------------------
// We now generalize the shifted grid to d dimensions. Let r be a parameter (typically chosen as r = ceil(c/ε) for some constant c).
// For each grid, we have a shift vector s ∈ {0, 1, ..., r-1}^d.
// For each point p ∈ P, and for each dimension j, let a_j = floor((p[j] - s[j]) / r).
// Then p lies in cell with index vector a = (a_1, ..., a_d), and the left boundary for coordinate j is s[j] + r * a_j.
// We remove p if for any dimension j, p[j] - (s[j] + r*a_j) < 1.
struct GridD {
    vector<int> shift;      // length d, each in {0, 1, ..., r-1}
    double cellLength;      // = r
    // Map: key is vector<int> of length d (the cell index), value is vector of points in that cell.
    map<vector<int>, vector<PointD>> cells;
};

struct VecIntCompare {
    bool operator()(const vector<int> &a, const vector<int> &b) const {
        return lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
    }
};

GridD buildGridD(const vector<PointD>& P, const vector<int>& shift, double r) {
    int d = shift.size();
    GridD grid;
    grid.shift = shift;
    grid.cellLength = r;
    for (const PointD &p : P) {
        vector<int> cellIdx(d, 0);
        bool valid = true;
        for (int j = 0; j < d; j++){
            int a = floor((p.coords[j] - shift[j]) / r);
            cellIdx[j] = a;
            double leftBoundary = shift[j] + r * a;
            if (p.coords[j] - leftBoundary < 1.0) {
                valid = false;
                break;
            }
        }
        if (!valid) continue;
        grid.cells[cellIdx].push_back(p);
    }
    return grid;
}

void sortGridCellsD(GridD &grid) {
    // For each cell, sort points lexicographically by coordinates.
    for (auto &kv : grid.cells) {
        sort(kv.second.begin(), kv.second.end(), [&](const PointD &a, const PointD &b) {
            return lexicographical_compare(a.coords.begin(), a.coords.end(),
                                           b.coords.begin(), b.coords.end());
        });
    }
}

// -------------------- Offline Algorithm for a d-Dimensional Grid --------------------
// We assume the cells are ordered lexicographically (by their index vector).
// We then compute a DP table over the cells similar to the 2D version.
vector<DPEntry> computeGridDP_d(const GridD &grid, int k, double epsilon) {
    // Collect cell indices and sort them lex.
    vector<vector<int>> cellIndices;
    for (auto &kv : grid.cells) {
        cellIndices.push_back(kv.first);
    }
    sort(cellIndices.begin(), cellIndices.end(), VecIntCompare());
    
    vector<vector<DPEntry>> cellDP;
    for (auto &idx : cellIndices) {
        const vector<PointD>& cellPoints = grid.cells.at(idx);
        vector<DPEntry> dpCell = computeCellDP_d(cellPoints, k);
        cout << "Cell [";
        for (int x : idx) cout << x << " ";
        cout << "] with " << cellPoints.size() << " points:" << endl;
        for (int j = 0; j < dpCell.size(); j++){
            cout << "  Best subset of size " << j << " has weight " << dpCell[j].totalWeight << endl;
        }
        cellDP.push_back(dpCell);
    }
    vector<vector<PointD>> cells;
    for (auto &idx : cellIndices) {
        cells.push_back(grid.cells.at(idx));
    }
    vector<DPEntry> dpResult = combineCells_d(cells, cellDP, k);
    return dpResult;
}

// -------------------- Offline Index for d Dimensions --------------------
// For a given d, let r = ceil(c/ε). (For instance, in R2, c was 2; for general d, one can choose an appropriate constant.)
// Then there are s = r^d grids. For each grid, the shift vector is chosen from {0,1,...,r-1}^d.
// We build each grid, compute its offline DP solution T(P(i)) and choose the best.
pair<vector<PointD>, double> offlineCDT_d(const vector<PointD>& P, int k, double epsilon, int d) {
    int r = ceil(2.0 / epsilon); // using the R2 value; this can be adjusted for higher d if needed.
    int s = pow(r, d);
    cout << "Using s = " << s << " grids in R^" << d << " (r = " << r << ")" << endl;
    pair<vector<PointD>, double> bestSolution = { {}, -1e9 };
    
    // We need to enumerate all shift vectors of length d with each coordinate in [0, r-1].
    vector<int> current(d, 0);
    function<void(int)> rec = [&](int pos) {
        if (pos == d) {
            cout << "\n=== Processing Grid with shift = [";
            for (int v : current) cout << v << " ";
            cout << "] ===" << endl;
            GridD grid = buildGridD(P, current, r);
            sortGridCellsD(grid);
            cout << "Non-empty cells in this grid:" << endl;
            for (auto &kv : grid.cells) {
                cout << "  Cell [";
                for (int v : kv.first) cout << v << " ";
                cout << "] has " << kv.second.size() << " points." << endl;
            }
            vector<DPEntry> dpResult = computeGridDP_d(grid, k, epsilon);
            // Choose the best among dp[0] ... dp[k]
            int best_j = 0;
            for (int j = 0; j <= k; j++){
                if (dpResult[j].totalWeight != NEG_INF && dpResult[j].totalWeight > dpResult[best_j].totalWeight)
                    best_j = j;
            }
            double candidateWeight = dpResult[best_j].totalWeight;
            cout << "Grid with shift [";
            for (int v : current) cout << v << " ";
            cout << "]: Best weight for up to k = " << k << " is " << candidateWeight << " (using " << best_j << " points)" << endl;
            if (candidateWeight > bestSolution.second) {
                bestSolution = { dpResult[best_j].subset, candidateWeight };
            }
            return;
        }
        for (int i = 0; i < r; i++){
            current[pos] = i;
            rec(pos+1);
        }
    };
    rec(0);
    return bestSolution;
}

// -------------------- Query Procedure for d Dimensions --------------------
// Given a query rectangle ρ (specified by two vectors: lower bounds and upper bounds),
// restrict P to those points in ρ and run offlineCDT_d on that subset.
pair<vector<PointD>, double> queryProcedure_d(const vector<PointD>& P, int k, double epsilon,
                                                const vector<double>& qLower, const vector<double>& qUpper) {
    cout << "\n[Query] Using query rectangle with lower bounds: [";
    for (double x : qLower) cout << x << " ";
    cout << "] and upper bounds: [";
    for (double x : qUpper) cout << x << " ";
    cout << "]" << endl;
    vector<PointD> Pq;
    for (const PointD &p : P){
        bool inside = true;
        for (int j = 0; j < qLower.size(); j++){
            if (p.coords[j] < qLower[j] || p.coords[j] > qUpper[j]){
                inside = false;
                break;
            }
        }
        if (inside)
            Pq.push_back(p);
    }
    cout << "[Query] " << Pq.size() << " points in query rectangle." << endl;
    return offlineCDT_d(Pq, k, epsilon, qLower.size());
}

// -------------------- Main Function --------------------

int main(){
    int n, k, d;
    double delta, epsilon;
    // Input: n k δ d, then n lines: weight c1 ... cd, then query rectangle: q1_left ... qd_left q1_right ... qd_right
    cin >> n >> k >> delta >> d;
    cin >> epsilon; // read ε
    vector<PointD> P(n);
    for (int i = 0; i < n; i++){
        cin >> P[i].weight;
        P[i].coords.resize(d);
        for (int j = 0; j < d; j++){
            cin >> P[i].coords[j];
        }
    }
    vector<double> qLower(d), qUpper(d);
    for (int j = 0; j < d; j++){
        cin >> qLower[j];
    }
    for (int j = 0; j < d; j++){
        cin >> qUpper[j];
    }
    
    cout << "Input parameters: n = " << n << ", k = " << k << ", δ = " << delta << ", d = " << d << endl;
    cout << "ε = " << epsilon << endl;
    cout << "Query Rectangle: lower bounds = [";
    for (double x : qLower) cout << x << " ";
    cout << "], upper bounds = [";
    for (double x : qUpper) cout << x << " ";
    cout << "]" << endl;
    
    // Scale coordinates if needed so that δ becomes 1.
    if (fabs(delta - 1.0) > 1e-9){
        cout << "Scaling coordinates by factor " << 1.0/delta << " so that δ becomes 1." << endl;
        for (auto &p : P){
            for (int j = 0; j < d; j++){
                p.coords[j] *= (1.0 / delta);
            }
        }
        for (int j = 0; j < d; j++){
            qLower[j] *= (1.0 / delta);
            qUpper[j] *= (1.0 / delta);
        }
    }
    
    cout << "\n=== Preprocessing (Offline Algorithm for R^d) ===" << endl;
    pair<vector<PointD>, double> sol = offlineCDT_d(P, k, epsilon, d);
    cout << "\n=== Offline CDT Solution over entire P in R^" << d << " ===" << endl;
    // Scale back before printing.
    if (fabs(delta - 1.0) > 1e-9){
        for (auto &p : sol.first){
            for (int j = 0; j < d; j++){
                p.coords[j] *= delta;
            }
        }
    }
    cout << "Best total weight: " << sol.second << "\nSelected points:" << endl;
    for (const auto &p : sol.first){
        cout << p.weight << " at (";
        for (double x : p.coords) cout << x << " ";
        cout << ")" << endl;
    }
    
    cout << "\n=== Query Phase for R^" << d << " ===" << endl;
    pair<vector<PointD>, double> qsol = queryProcedure_d(P, k, epsilon, qLower, qUpper);
    cout << "\n=== Query Result in R^" << d << " ===" << endl;
    // Scale back query results.
    if (fabs(delta - 1.0) > 1e-9){
        for (auto &p : qsol.first){
            for (int j = 0; j < d; j++){
                p.coords[j] *= delta;
            }
        }
    }
    cout << "Best total weight: " << qsol.second << "\nSelected points:" << endl;
    for (const auto &p : qsol.first){
        cout << p.weight << " at (";
        for (double x : p.coords) cout << x << " ";
        cout << ")" << endl;
    }
    
    return 0;
}
