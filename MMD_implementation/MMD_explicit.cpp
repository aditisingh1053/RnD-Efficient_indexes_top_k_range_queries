#include <bits/stdc++.h>
#include <cmath>
// #include "ECP.h"
using namespace std;

// ---------------------------------------------------------------------------
// Data structures and helper functions

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

// Brute-force computation of the closest pair distance (O(n^2)).
double bruteForceClosestPair(const vector<Point>& pts) {
    double best = numeric_limits<double>::max();
    int n = pts.size();
    if(n < 2) return 0.0;
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++){
            best = min(best, euclideanDistance(pts[i], pts[j]));
        }
    }
    return best;
}

// Compute the minimum weight among points in a set.
double minWeight(const vector<Point>& pts) {
    double m = numeric_limits<double>::max();
    for (const auto &p: pts) {
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

// ---------------------------------------------------------------------------
// Main Algorithm
int main(){
    // Read input: number of points n, dimension d, desired subset size k,
    // and lambda (the weight multiplier in the utility function).
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
    double init_mu = bruteForceClosestPair(initialPk);
    
    // Prepare an array of packings, one for each radius level.
    vector<Packing> packings(J);
    // Also prepare a parallel vector of dynamic ECPs for each packing level.
    // vector<DynamicECP> ecp_packings(J);
    
    for (int i = 0; i < J; i++){
        // Initialize the radius as: r_i = (1+ε)^i * μ(Pk) / alpha.
        packings[i].r = pow(1 + epsilon, i) * init_mu / alpha;
        // A, N, D remain empty initially.
    }
    
    // Our current best candidate solution (subset of k points)
    vector<Point> bestSolution = initialPk;
    // Compute utility f(T)=μ(T)+λ·minWeight(T)
    double bestUtility = bruteForceClosestPair(initialPk) + lambda * minWeight(initialPk);
    
    // --- Process the remaining points in the stream ---
    for (int t = k; t < n; t++){
        Point p = stream[t];
        // Update each packing level with the new point.
        for (int i = 0; i < J; i++){
            bool found = false;
            // Check in A.
            for (const auto &q : packings[i].A) {
                if(euclideanDistance(p, q) <= packings[i].r) { 
                    found = true; 
                    break; 
                }
            }
            // Check in N if not found.
            if (!found) {
                for (const auto &q : packings[i].N) {
                    if(euclideanDistance(p, q) <= packings[i].r) { 
                        found = true; 
                        break; 
                    }
                }
            }
            if(found) continue;  // p is covered by A∪N
            
            // Otherwise, add p to the buffer N.
            packings[i].N.push_back(p);
            
            // If A∪N now has exactly k points, update the radius.
            if(packings[i].A.size() + packings[i].N.size() == (size_t)k){
                // Build the combined set from A and N.
                vector<Point> combined = packings[i].A;
                combined.insert(combined.end(), packings[i].N.begin(), packings[i].N.end());
                
                // Instead of calling bruteForceClosestPair, use the dynamic ECP:
                // First, insert all points in A into the dynamic structure.
                // (For simplicity we assume that when a packing level is updated for the first time,
                //  its dynamic ECP is empty.)
                // if(ecp_packings[i].getPoints().empty()){
                //     for (const auto &q : packings[i].A) {
                //         ecp_packings[i].insert(q);
                //     }
                // }
                // // Also insert any new points from N.
                // for (const auto &q : packings[i].N) {
                //     ecp_packings[i].insert(q);
                // }
                // Query the current closest pair distance.
                // double current_mu = ecp_packings[i].getClosestPair();
                double current_mu=bruteForceClosestPair(combined);
                
                // Increase r until the invariant holds: find smallest integer m such that α^m * r > current_mu.
                int m = 0;
                while(pow(alpha, m) * packings[i].r <= current_mu) 
                    m++;
                packings[i].r = pow(alpha, m) * packings[i].r;
                
                // Absorb all points from N into A and clear buffers N and D.
                packings[i].A.insert(packings[i].A.end(), packings[i].N.begin(), packings[i].N.end());
                packings[i].N.clear();
                packings[i].D.clear();
                // (Our dynamic ECP already has these points.)
                
                // --- Pruning process using dynamic ECP ---
                bool pruning = true;
                while(pruning && packings[i].A.size() >= 2){
                    double cp = bruteForceClosestPair(packings[i].A);
                    if(cp >= packings[i].r) {
                        pruning = false;
                    } else {
                        // Find a pair (p, q) with distance < r and remove one of them.
                        // Here, we scan the current set A brute-force.
                        bool removed = false;
                        for (size_t a = 0; a < packings[i].A.size() && !removed; a++){
                            for (size_t b = a + 1; b < packings[i].A.size() && !removed; b++){
                                if(euclideanDistance(packings[i].A[a], packings[i].A[b]) < packings[i].r){
                                    // Remove the second point.
                                    Point toRemove = packings[i].A[b];
                                    packings[i].D.push_back(toRemove);
                                    // Erase from A.
                                    packings[i].A.erase(packings[i].A.begin() + b);
                                    // Also remove from the dynamic ECP.
                                    // ecp_packings[i].remove(toRemove);
                                    removed = true;
                                }
                            }
                        }
                        if(!removed) break; // Safety break in case no removal occurs.
                    }
                }
            } // end if(A∪N == k)
        } // end for each packing level
        
        // --- Select the candidate solution from the current packings ---
        // Choose the packing level with the largest radius.
        int bestIndex = 0;
        double max_r = -1.0;
        for (int i = 0; i < J; i++){
            if(packings[i].r > max_r){
                max_r = packings[i].r;
                bestIndex = i;
            }
        }
        // Form candidate = A ∪ D.
        vector<Point> candidate = packings[bestIndex].A;
        candidate.insert(candidate.end(), packings[bestIndex].D.begin(), packings[bestIndex].D.end());
        // If candidate size is not equal to k, fall back to the initial prefix.
        if(candidate.size() != (size_t)k)
            candidate = initialPk;
        
        // Compute candidate utility.
        double candidate_mu = bruteForceClosestPair(candidate);
        double candidate_min_weight = minWeight(candidate);
        double candidate_utility = candidate_mu + lambda * candidate_min_weight;
        
        if(candidate_utility > bestUtility) {
            bestUtility = candidate_utility;
            bestSolution = candidate;
        }
    }
    
    // --- Output the best solution subset ---
    // Output the coordinates and weight of each point in the best subset.
    for (const auto &pt : bestSolution) {
        for (size_t j = 0; j < pt.coord.size(); j++){
            cout << pt.coord[j] << " ";
        }
        cout << pt.weight << "\n";
    }
    
    return 0;
}
