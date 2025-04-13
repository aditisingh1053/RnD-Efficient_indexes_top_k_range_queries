#include <bits/stdc++.h>
using namespace std;

// -------------------- Data Structures --------------------

// In 2D, each point has a weight, and x,y coordinates.
struct Point2D {
    double weight;
    double x, y;
};

bool compareByWeight2D(const Point2D &a, const Point2D &b) {
    return a.weight > b.weight;
}

void printPoint2D(const Point2D &p) {
    cout << "(" << p.x << ", " << p.y << ") -> " << p.weight;
}

// DPEntry stores a candidate subset (solution) along with its total weight.
struct DPEntry {
    vector<Point2D> subset;
    double totalWeight;
};

const double NEG_INF = -1e9;

// -------------------- Brute-force for a Single Cell in 2D --------------------
// Enumerate all subsets of cellPoints of size up to min(j_max, |cellPoints|).
// Only accept those subsets that are 1-diverse, i.e., any two points in the subset
// have Euclidean distance at least 1 (recall that we scale so that δ = 1).
// Returns a DP table (vector of DPEntry) of size (j_max+1), where dp[j] is the best subset of exactly j points.
vector<DPEntry> computeCellDP2D(const vector<Point2D>& cellPoints, int j_max) {
    int n = cellPoints.size();
    vector<DPEntry> best(j_max+1, { {}, NEG_INF });
    best[0] = { {}, 0.0 };

    int total = 1 << n;
    for (int mask = 1; mask < total; mask++){
        vector<Point2D> subset;
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
        for (int i = 0; i < (int)subset.size() && valid; i++){
            for (int j = i+1; j < (int)subset.size(); j++){
                double dx = subset[i].x - subset[j].x;
                double dy = subset[i].y - subset[j].y;
                double dist = sqrt(dx*dx + dy*dy);
                if (dist < 1.0) {
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

// -------------------- Dynamic Programming Over Cells in 2D --------------------
// Given a sequence of cells (each cell is a vector<Point2D>) and their DP tables (computed via computeCellDP2D),
// combine them (in lexicographic order) to compute T(P(i)) for a grid.
// This is analogous to a knapsack DP: dp[i][j] is the best solution using the first i cells with exactly j points.
vector<DPEntry> combineCells2D(const vector<vector<Point2D>>& cells, const vector<vector<DPEntry>>& cellDP, int k) {
    int m = cells.size();
    vector<vector<DPEntry>> dp(m+1, vector<DPEntry>(k+1, { {}, NEG_INF }));
    dp[0][0] = { {}, 0.0 };
    for (int i = 1; i <= m; i++){
        for (int j = 0; j <= k; j++){
            dp[i][j] = dp[i-1][j]; // Option: take 0 points from current cell.
            int maxPtsInCell = cellDP[i-1].size()-1; // p from 0 to maxPtsInCell.
            for (int p = 1; p <= min(j, maxPtsInCell); p++){
                if (dp[i-1][j-p].totalWeight == NEG_INF || cellDP[i-1][p].totalWeight == NEG_INF)
                    continue;
                double candidateWeight = dp[i-1][j-p].totalWeight + cellDP[i-1][p].totalWeight;
                if (candidateWeight > dp[i][j].totalWeight){
                    vector<Point2D> candidateSubset = dp[i-1][j-p].subset;
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

// -------------------- 2D Grid Construction --------------------
// We extend the 1D shifted grid to 2D. Let r = ceil(2/ε) and s = r^2.
// For each grid index i (0 ≤ i < s), if i = α * r + β (with α, β ∈ [0, r-1]),
// then grid G_i is defined by shifts: shiftX = α, shiftY = β, and its cells are
// { (α + r*a, β + r*b) : a, b ∈ Z }.
// Each grid cell is indexed by its bottom-left vertex.
// For each cell, remove points that lie within distance 1 from its bottom or left boundary.
struct Grid2D {
    int shiftX, shiftY; // determined by α and β.
    double cellLength;  // = r
    // Map: cell index (a,b) -> vector of points in that cell.
    map<pair<int,int>, vector<Point2D>> cells;
};

Grid2D buildGrid2D(const vector<Point2D>& P, int shiftX, int shiftY, double r) {
    Grid2D grid;
    grid.shiftX = shiftX;
    grid.shiftY = shiftY;
    grid.cellLength = r;
    for (const Point2D &p : P){
        int a = floor((p.x - shiftX) / r);
        int b = floor((p.y - shiftY) / r);
        double leftX = shiftX + r * a;
        double bottomY = shiftY + r * b;
        // Remove p if its distance from the left boundary is < 1 or from bottom boundary is < 1.
        if ((p.x - leftX) < 1.0 || (p.y - bottomY) < 1.0)
            continue;
        grid.cells[{a, b}].push_back(p);
    }
    return grid;
}

void sortGridCells2D(Grid2D &grid) {
    for (auto &kv : grid.cells){
        sort(kv.second.begin(), kv.second.end(), [](const Point2D &a, const Point2D &b){
            if (a.x != b.x) return a.x < b.x;
            return a.y < b.y;
        });
    }
}

// -------------------- Offline Algorithm for a 2D Grid --------------------
// Given a 2D grid (one of the grids G_i), partition the grid into non-empty cells
// (ordered lexicographically by cell index) and compute T(P(i)) via dynamic programming.
vector<DPEntry> computeGridDP2D(const Grid2D &grid, int k, double epsilon) {
    vector<pair<int,int>> cellIndices;
    for (auto &kv : grid.cells){
        cellIndices.push_back(kv.first);
    }
    // Sort lexicographically: first by a, then by b.
    sort(cellIndices.begin(), cellIndices.end(), [](const pair<int,int> &p1, const pair<int,int> &p2){
        if (p1.first != p2.first) return p1.first < p2.first;
        return p1.second < p2.second;
    });
    
    vector<vector<DPEntry>> cellDP;
    for (auto &idx : cellIndices){
        const vector<Point2D>& cellPoints = grid.cells.at(idx);
        vector<DPEntry> dpCell = computeCellDP2D(cellPoints, k);
        cout << "Cell (" << idx.first << "," << idx.second << ") with " << cellPoints.size() << " points:" << endl;
        for (int j = 0; j < dpCell.size(); j++){
            cout << "  Best subset of size " << j << " has weight " << dpCell[j].totalWeight << endl;
        }
        cellDP.push_back(dpCell);
    }
    vector<vector<Point2D>> cells;
    for (auto &idx : cellIndices){
        cells.push_back(grid.cells.at(idx));
    }
    vector<DPEntry> dpResult = combineCells2D(cells, cellDP, k);
    return dpResult;
}

// -------------------- Offline Index for 2D --------------------
// Let r = ceil(2/ε) and s = r^2. For each grid G_i, build the grid using shifts:
// if i = α*r + β then shiftX = α, shiftY = β. Compute T(P(i)) for that grid and pick the best.
pair<vector<Point2D>, double> offlineCDT2D(const vector<Point2D>& P, int k, double epsilon) {
    int r = ceil(2.0 / epsilon);
    int s = r * r;
    cout << "Using s = " << s << " grids in R2 (r = " << r << ")" << endl;
    pair<vector<Point2D>, double> bestSolution = { {}, -1e9 };
    for (int i = 0; i < s; i++){
        int alpha = i / r, beta = i % r;
        cout << "\n=== Processing Grid G" << i << " (shiftX = " << alpha << ", shiftY = " << beta << ") ===" << endl;
        Grid2D grid = buildGrid2D(P, alpha, beta, r);
        sortGridCells2D(grid);
        cout << "Non-empty cells in G" << i << ":" << endl;
        for (auto &kv : grid.cells){
            cout << "  Cell (" << kv.first.first << "," << kv.first.second << ") has " << kv.second.size() << " points." << endl;
        }
        vector<DPEntry> dpResult = computeGridDP2D(grid, k, epsilon);
        // Instead of simply reading dpResult[k], choose the best over j = 0...k.
        int best_j = 0;
        for (int j = 0; j <= k; j++){
            if (dpResult[j].totalWeight != NEG_INF && dpResult[j].totalWeight > dpResult[best_j].totalWeight)
                best_j = j;
        }
        double weightCandidate = dpResult[best_j].totalWeight;
        cout << "Grid G" << i << ": Best weight for up to k = " << k << " is " << weightCandidate << " (using " << best_j << " points)" << endl;
        if (weightCandidate > bestSolution.second){
            bestSolution = { dpResult[best_j].subset, weightCandidate };
        }
    }
    return bestSolution;
}

// -------------------- Query Procedure for 2D --------------------
// Given a query rectangle ρ = [x_left, y_bottom, x_right, y_top],
// extend it to the smallest F-aligned rectangle (F: fine grid of cell size 1/r),
// then restrict P to those points in ρ and run the offline algorithm.
pair<vector<Point2D>, double> queryProcedure2D(const vector<Point2D>& P, int k, double epsilon, const tuple<double,double,double,double>& queryRect) {
    double xL, yL, xR, yR;
    tie(xL, yL, xR, yR) = queryRect;
    cout << "\n[Query] Using query rectangle [(" << xL << ", " << yL << "), (" << xR << ", " << yR << ")]" << endl;
    vector<Point2D> Pq;
    for (const Point2D &p : P){
        if (p.x >= xL && p.x <= xR && p.y >= yL && p.y <= yR)
            Pq.push_back(p);
    }
    cout << "[Query] " << Pq.size() << " points in query rectangle." << endl;
    return offlineCDT2D(Pq, k, epsilon);
}

// -------------------- Main Function --------------------

int main(){
    int n, k, d;
    double delta, epsilon;
    // Input: n k δ d, then n lines: weight x y, then query rectangle: x_left y_bottom x_right y_top
    cin >> n >> k >> delta >> d;
    // cin >> epsilon; // read ε
    epsilon = 0.1; // default ε
    vector<Point2D> P(n);
    for (int i = 0; i < n; i++){
        cin >> P[i].weight >> P[i].x >> P[i].y;
    }
    double qxL, qyL, qxR, qyR;
    cin >> qxL  >> qxR >> qyL >> qyR;
    tuple<double,double,double,double> queryRect = make_tuple(qxL, qyL, qxR, qyR);
    
    cout << "Input parameters: n = " << n << ", k = " << k << ", δ = " << delta << ", d = " << d << endl;
    cout << "ε = " << epsilon << ", Query Rectangle = [(" << qxL << ", " << qyL << "), (" << qxR << ", " << qyR << ")]" << endl;
    
    // Scale coordinates if needed so that δ becomes 1.
    if (fabs(delta - 1.0) > 1e-9){
        cout << "Scaling coordinates by factor " << 1.0/delta << " so that δ becomes 1." << endl;
        for (auto &p : P){
            p.x *= (1.0 / delta);
            p.y *= (1.0 / delta);
        }
        qxL *= (1.0 / delta); qyL *= (1.0 / delta);
        qxR *= (1.0 / delta); qyR *= (1.0 / delta);
        queryRect = make_tuple(qxL, qyL, qxR, qyR);
    }
    
    cout << "\n=== Preprocessing (Offline Algorithm for R2) ===" << endl;
    pair<vector<Point2D>, double> sol = offlineCDT2D(P, k, epsilon);
    cout << "\n=== Offline CDT Solution over entire P in R2 ===" << endl;
    // Scale back before printing
    for (auto &p : sol.first){
        p.x *= delta;
        p.y *= delta;
    }
    cout << "Best total weight: " << sol.second << "\nSelected points:" << endl;
    for (const auto &p : sol.first){
        cout << p.weight << " at (" << p.x << ", " << p.y << ")" << endl;
    }
    
    cout << "\n=== Query Phase for R2 ===" << endl;
    pair<vector<Point2D>, double> qsol = queryProcedure2D(P, k, epsilon, queryRect);
    cout << "\n=== Query Result in R2 ===" << endl;
    // Scale back query results.
    for (auto &p : qsol.first){
        p.x *= delta;
        p.y *= delta;
    }
    cout << "Best total weight: " << qsol.second << "\nSelected points:" << endl;
    for (const auto &p : qsol.first){
        cout << p.weight << " at (" << p.x << ", " << p.y << ")" << endl;
    }
    
    return 0;
}
