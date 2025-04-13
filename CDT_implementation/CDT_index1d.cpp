#include <bits/stdc++.h>
using namespace std;

// ---------- Data Structures and Constants ----------

const double NEG_INF = -1e9;

struct Point {
    double weight;
    double coord; // 1D coordinate (used when d == 1)
    vector<double> coords; // General coordinate vector (d dimensions)
};

bool compareByWeight(const Point &a, const Point &b) {
    return a.weight > b.weight;
}

// For debug printing of a point.
void printPoint(const Point &p) {
    cout << "(" << p.coord << ") -> " << p.weight;
}

// DPEntry stores a candidate subset (solution) along with its total weight.
struct DPEntry {
    vector<Point> subset;
    double totalWeight;
};

DPEntry makeEmptyEntry() {
    return { vector<Point>(), 0.0 };
}

// ---------- Brute-force for a single cell ----------
// Enumerate all subsets of cellPoints of size up to min(j_max, |cellPoints|)
// Only accept those subsets that are 1-diverse (i.e., pairwise distance >= 1)
// Return a DP table of size (j_max+1), where dp[j] is the best subset of exactly j points.
vector<DPEntry> computeCellDP(const vector<Point>& cellPoints, int j_max) {
    int n = cellPoints.size();
    vector<DPEntry> best(j_max+1, { {}, NEG_INF });
    best[0] = { {}, 0.0 };

    int total = 1 << n;
    for (int mask = 1; mask < total; mask++) {
        vector<Point> subset;
        double sumW = 0;
        int count = 0;
        for (int i = 0; i < n; i++) {
            if (mask & (1 << i)) {
                subset.push_back(cellPoints[i]);
                sumW += cellPoints[i].weight;
                count++;
            }
        }
        if (count > j_max) continue;
        bool valid = true;
        for (int i = 0; i < subset.size() && valid; i++) {
            for (int j = i+1; j < subset.size(); j++) {
                if (fabs(subset[i].coord - subset[j].coord) < 1.0) {
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

// ---------- Dynamic Programming over cells ----------
// Combine cells (left-to-right) to compute T(P(i)) for a grid.
vector<DPEntry> combineCells(const vector<vector<Point>>& cells, const vector<vector<DPEntry>>& cellDP, int k) {
    int m = cells.size();
    vector<vector<DPEntry>> dp(m+1, vector<DPEntry>(k+1, { {}, NEG_INF }));
    dp[0][0] = { {}, 0.0 };
    
    for (int i = 1; i <= m; i++) {
        for (int j = 0; j <= k; j++) {
            dp[i][j] = dp[i-1][j];
            int maxPointsInCell = cellDP[i-1].size()-1;
            for (int p = 1; p <= min(j, maxPointsInCell); p++) {
                if (dp[i-1][j-p].totalWeight == NEG_INF || cellDP[i-1][p].totalWeight == NEG_INF)
                    continue;
                double candidateWeight = dp[i-1][j-p].totalWeight + cellDP[i-1][p].totalWeight;
                if (candidateWeight > dp[i][j].totalWeight) {
                    vector<Point> candidateSubset = dp[i-1][j-p].subset;
                    candidateSubset.insert(candidateSubset.end(), cellDP[i-1][p].subset.begin(), cellDP[i-1][p].subset.end());
                    dp[i][j] = { candidateSubset, candidateWeight };
                }
            }
        }
    }
    cout << "DP over cells complete. Final values:" << endl;
    for (int j = 0; j <= k; j++) {
        cout << "  For j = " << j << ", best weight = " << dp[m][j].totalWeight << endl;
    }
    return dp[m];
}

// ---------- Grid Construction ----------
// For 1D, we use shifted grids. Let r = ceil(1/ε). For each grid Gi (i = 0, 1, ..., r-1),
// grid cell intervals are [L, L + r), where L = shift + r*a for some a ∈ Z.
// For grid Gi, define P(i) as follows: for each point p in P, if p lies in a cell,
// remove it if its distance from the left endpoint is < 1.
struct Grid {
    int shift;         // i in 0 ... r-1
    double cellLength; // = r
    map<int, vector<Point>> cells;
};

Grid buildGrid(const vector<Point>& P, int shift, double r) {
    Grid grid;
    grid.shift = shift;
    grid.cellLength = r;
    for (const Point& p : P) {
        int a = floor((p.coord - shift) / r);
        double leftEndpoint = shift + r * a;
        if (p.coord - leftEndpoint < 1.0) continue;
        grid.cells[a].push_back(p);
    }
    return grid;
}

void sortGridCells(Grid &grid) {
    for (auto &kv : grid.cells) {
        sort(kv.second.begin(), kv.second.end(), [](const Point &a, const Point &b) {
            return a.coord < b.coord;
        });
    }
}

// ---------- Offline Algorithm to Compute CD(P,k) for a Grid ----------
// Partition grid P(i) into non-empty cells C1,...,Cm (ordered by cell index)
// and compute T(P(i)) via DP.
vector<DPEntry> computeGridDP(const Grid &grid, int k, double epsilon) {
    vector<vector<DPEntry>> cellDP;
    vector<int> cellIndices;
    for (auto &kv : grid.cells) {
        cellIndices.push_back(kv.first);
    }
    sort(cellIndices.begin(), cellIndices.end());
    
    for (int idx : cellIndices) {
        const vector<Point>& cellPoints = grid.cells.at(idx);
        vector<DPEntry> dpCell = computeCellDP(cellPoints, k);
        cout << "Cell " << idx << " with " << cellPoints.size() << " points:" << endl;
        for (int j = 0; j < dpCell.size(); j++) {
            cout << "  Best subset of size " << j << " has weight " << dpCell[j].totalWeight << endl;
        }
        cellDP.push_back(dpCell);
    }
    vector<vector<Point>> cells; // only for passing to combineCells (not used further)
    for (int idx : cellIndices) {
        cells.push_back(grid.cells.at(idx));
    }
    vector<DPEntry> dpResult = combineCells(cells, cellDP, k);
    return dpResult;
}

// ---------- Offline Algorithm: Try All Grids and Pick Best ----------
// For r = ceil(1/ε), build grids G0,...,G_{r-1}. For each grid, compute T(P(i), k)
// and pick the best (highest total weight).
pair<vector<Point>, double> offlineCDT(const vector<Point>& P, int k, double epsilon) {
    int r = ceil(1.0 / epsilon);
    cout << "Using r = " << r << " grids (shift = 0 to " << r-1 << ")" << endl;
    pair<vector<Point>, double> bestSolution = { {}, NEG_INF };
    for (int i = 0; i < r; i++) {
        cout << "\n=== Processing Grid G" << i << " ===" << endl;
        Grid grid = buildGrid(P, i, r);
        sortGridCells(grid);
        cout << "Non-empty cells in G" << i << ":" << endl;
        for (auto &kv : grid.cells) {
            cout << "  Cell " << kv.first << " has " << kv.second.size() << " points." << endl;
        }
        vector<DPEntry> dpResult = computeGridDP(grid, k, epsilon);
        double weightCandidate = dpResult[k].totalWeight;
        cout << "Grid G" << i << ": Best weight for k = " << k << " is " << weightCandidate << endl;
        if (weightCandidate > bestSolution.second) {
            bestSolution = { dpResult[k].subset, weightCandidate };
        }
    }
    return bestSolution;
}

// ---------- Query Procedure ----------
// In the index version, one builds an index on each grid. Here, for simplicity,
// we assume the query interval is the entire range (or extend it to F-aligned).
// We then run the offline algorithm on P ∩ queryInterval.
pair<vector<Point>, double> queryProcedure(const vector<Point>& P, int k, double epsilon, const pair<double,double>& queryInterval) {
    cout << "\n[Query] Using query interval [" << queryInterval.first << ", " << queryInterval.second << "]" << endl;
    vector<Point> Pq;
    for (const Point &p : P) {
        if (p.coord >= queryInterval.first && p.coord <= queryInterval.second)
            Pq.push_back(p);
    }
    cout << "[Query] " << Pq.size() << " points in query interval." << endl;
    return offlineCDT(Pq, k, epsilon);
}

// ---------- Main Function ----------

int main(){
    int n, k, d;
    d=1;
    double delta, epsilon;
    // Input: n k δ d, then n lines: weight coordinate, then query interval (left right)
    cin >> n >> k >> delta;
    int storedelta=delta;
    vector<Point> P(n);
    for (int i = 0; i < n; i++){
        cin >> P[i].weight;
        P[i].coords.resize(d);
        for (int j = 0; j < d; j++){
            cin >> P[i].coords[j];
        }
        if(d == 1) P[i].coord = P[i].coords[0];
    }
    epsilon=0.1;
    cout << "Input parameters: n = " << n << ", k = " << k << ", δ = " << delta << ", d = " << d << endl;
    pair<double, double> queryInterval;
    cin >> queryInterval.first >> queryInterval.second;
    
    // Scale coordinates if needed so that δ becomes 1.
    if (fabs(delta - 1.0) > 1e-9) {
        cout << "Scaling coordinates by factor " << 1.0/delta << " so that δ becomes 1." << endl;
        for (auto &p : P) {
            p.coord *= (1.0 / delta);
            for (auto &x : p.coords) {
                x *= (1.0 / delta);
            }
        }
        queryInterval.first *= (1.0 / delta);
        queryInterval.second *= (1.0 / delta);
    }
    
    cout << "\n=== Preprocessing (Offline Algorithm) ===" << endl;
    pair<vector<Point>, double> sol = offlineCDT(P, k, epsilon);
    cout << "\n=== Offline CDT Solution over entire P ===" << endl;
    cout << "Best total weight: " << sol.second << "\nSelected points:" << endl;
    for (const auto &p : sol.first) {
        cout << p.weight << " at " << p.coord << endl;
    }
    
    cout << "\n=== Query Phase ===" << endl;
    pair<vector<Point>, double> qsol = queryProcedure(P, k, epsilon, queryInterval);
    cout << "\n=== Query Result ===" << endl;
    cout << "Best total weight: " << qsol.second << "\nSelected points:" << endl;
    for (const auto &p : qsol.first) {
        cout << p.weight << " at " << p.coord*storedelta << endl;
    }
    
    return 0;
}
