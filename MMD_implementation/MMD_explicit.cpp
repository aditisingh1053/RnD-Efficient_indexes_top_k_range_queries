#include <bits/stdc++.h>
#include <cmath>
using namespace std;

// --- Data structures and helper functions ---

// Structure to represent a point.
struct Point {
    vector<double> coord;  // d-dimensional coordinates
    double weight;
    int id;
};

// Compute Euclidean distance between two points.
double euclideanDistance(const Point &a, const Point &b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.coord.size(); i++){
        double diff = a.coord[i] - b.coord[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}

// Brute-force computation of the closest pair distance (μ) in a set of points.
double closestPairDistance(const vector<Point>& pts) {
    double best = numeric_limits<double>::max();
    for (size_t i = 0; i < pts.size(); i++){
        for (size_t j = i + 1; j < pts.size(); j++){
            best = min(best, euclideanDistance(pts[i], pts[j]));
        }
    }
    // If only one point is available, define the value as 0.
    return (pts.size() < 2 ? 0.0 : best);
}

// Compute the minimum weight among points in a set.
double minWeight(const vector<Point>& pts) {
    double m = numeric_limits<double>::max();
    for (auto &p: pts) {
        m = min(m, p.weight);
    }
    return m;
}

// Packing structure for each radius level.
struct Packing {
   double r;                // current radius value for the packing
   vector<Point> A;         // maintained subset A
   vector<Point> N;         // recent insertions buffer N
   vector<Point> D;         // deleted points (from pruning)
};

// --- Main Algorithm ---

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    // Read input: number of points n, dimension d, desired subset size k,
    // and λ (the weight multiplier in the utility function).
    int n, d, k;
    double lambda; 
    double epsilon = 0.1; // constant ε (can be adjusted if needed)
    cin >> n >> d >> k >> lambda;
    
    vector<Point> stream(n);
    for (int i = 0; i < n; i++){
        stream[i].coord.resize(d);
        for (int j = 0; j < d; j++){
            cin >> stream[i].coord[j];
        }
        cin >> stream[i].weight;
        stream[i].id = i + 1;
    }
    
    // Sort the points in non-increasing order by weight.
    sort(stream.begin(), stream.end(), [](const Point &a, const Point &b) {
        return a.weight > b.weight;
    });
    
    // --- Compute the number of packing levels J and compute alpha ---
    // We need the smallest integer J such that (1+ε)^J >= 1 + 1/ε.
    int J = 0;
    while(pow(1 + epsilon, J) < 1 + 1/epsilon) J++;
    // In the algorithm, alpha is defined as (1+ε)^J.
    double alpha = pow(1 + epsilon, J);
    
    // --- Initialization using the first k points ---
    // Assume the stream has at least k points.
    vector<Point> initialPk(stream.begin(), stream.begin() + k);
    double init_mu = closestPairDistance(initialPk);
    
    // Prepare an array of packings, one for each radius level.
    vector<Packing> packings(J);
    for (int i = 0; i < J; i++){
        // Initialize the radius as: r_i = (1+ε)^i * μ(Pk) / alpha.
        packings[i].r = pow(1 + epsilon, i) * init_mu / alpha;
        // A, N, D remain empty initially.
    }
    
    // Our current best candidate solution (subset of k points)
    vector<Point> bestSolution = initialPk;
    // Compute utility f(T)=μ(T)+λ·min{w(p): p in T}
    double bestUtility = closestPairDistance(initialPk) + lambda * minWeight(initialPk);
    
    // --- Process the remaining points in the stream ---
    for (int t = k; t < n; t++){
        Point p = stream[t];
        // Update each packing level with the new point.
        for (int i = 0; i < J; i++){
            // Check if there exists some point in A∪N within distance r.
            bool found = false;
            for (const auto &q : packings[i].A) {
                if(euclideanDistance(p, q) <= packings[i].r) { 
                    found = true; 
                    break; 
                }
            }
            if(!found) {
                for (const auto &q : packings[i].N) {
                    if(euclideanDistance(p, q) <= packings[i].r) { 
                        found = true; 
                        break; 
                    }
                }
            }
            // If p is already "covered" by A∪N, no update is needed.
            if(found) continue;
            
            // Otherwise, add p to the buffer N.
            packings[i].N.push_back(p);
            
            // If the union A∪N now has exactly k points, update the radius.
            if(packings[i].A.size() + packings[i].N.size() == (size_t)k){
                // Build the combined set from A and N.
                vector<Point> combined = packings[i].A;
                combined.insert(combined.end(), packings[i].N.begin(), packings[i].N.end());
                double current_mu = closestPairDistance(combined);
                
                // Increase the current radius r until the invariant holds:
                // find the smallest integer m such that α^m * r > μ(A∪N)
                int m = 0;
                while(pow(alpha, m) * packings[i].r <= current_mu) 
                    m++;
                packings[i].r = pow(alpha, m) * packings[i].r;
                
                // Absorb all recent insertions into A and clear the buffers.
                packings[i].A.insert(packings[i].A.end(), packings[i].N.begin(), packings[i].N.end());
                packings[i].N.clear();
                packings[i].D.clear();
                
                // --- Pruning process ---
                // Remove points from A until its closest pair distance is at least r.
                bool pruning = true;
                while(pruning && packings[i].A.size() >= 2){
                    double cp = closestPairDistance(packings[i].A);
                    if(cp >= packings[i].r){
                        pruning = false;
                    } else {
                        // Find a pair (p, q) with distance < r and remove one of them.
                        bool removed = false;
                        for (size_t a = 0; a < packings[i].A.size() && !removed; a++){
                            for (size_t b = a + 1; b < packings[i].A.size() && !removed; b++){
                                if(euclideanDistance(packings[i].A[a], packings[i].A[b]) < packings[i].r){
                                    // Remove the second point and record it in D.
                                    packings[i].D.push_back(packings[i].A[b]);
                                    packings[i].A.erase(packings[i].A.begin() + b);
                                    removed = true;
                                }
                            }
                        }
                        if(!removed) break; // Safety break
                    }
                }
            } // end: if (|A∪N| == k)
        } // end: for each packing level
        
        // --- Select the candidate solution from the current packings ---
        // Pick the packing level with the largest radius.
        int bestIndex = 0;
        double max_r = -1.0;
        for (int i = 0; i < J; i++){
            if(packings[i].r > max_r){
                max_r = packings[i].r;
                bestIndex = i;
            }
        }
        // Form the candidate T = A ∪ D.
        vector<Point> candidate = packings[bestIndex].A;
        candidate.insert(candidate.end(), packings[bestIndex].D.begin(), packings[bestIndex].D.end());
        // If candidate size does not equal k, fall back to the initial prefix.
        if(candidate.size() != (size_t)k)
            candidate = initialPk;
        
        // Compute candidate utility: f(candidate)=μ(candidate)+λ·minWeight(candidate).
        double candidate_mu = closestPairDistance(candidate);
        double candidate_min_weight = minWeight(candidate);
        double candidate_utility = candidate_mu + lambda * candidate_min_weight;
        
        if(candidate_utility > bestUtility) {
            bestUtility = candidate_utility;
            bestSolution = candidate;
        }
    }
    
    // --- Output the best solution subset ---
    // (Here, we output the coordinates and weight of each point in the best subset.)
    for (const auto &pt : bestSolution) {
        for (size_t j = 0; j < pt.coord.size(); j++){
            cout << pt.coord[j] << " ";
        }
        cout << pt.weight << "\n";
    }
    
    return 0;
}
