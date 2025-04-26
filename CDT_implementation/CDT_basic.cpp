#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <map>
#include <algorithm>

using namespace std;

typedef pair<double, double> Point;

double euclidean_distance(const Point &p1, const Point &p2) {
    return sqrt(pow(p1.first - p2.first, 2) + pow(p1.second - p2.second, 2));
}
void printSet(const set<pair<int, Point>> &S) {
    cout << "Selected points:\n";
    for (const auto &[w, p] : S) {
        cout << "(" << p.first << ", " << p.second << ") -> weight " << w << "\n";
    }
}

vector<pair<double,Point>> greedy_cdt(vector<pair<double, Point>> &w, int k, double delta) {
    vector<pair<double,Point>> S; // Set of selected points
    set<pair<double,Point>> point_set(w.begin(), w.end()); // Convert to set for efficient removal
    // printSet(point_set);
    for (int i = 0; i < k; ++i) {

        if (point_set.empty()) break;

        auto last_element = --point_set.end();
        // cout<<(*last_element).first<<endl;
        S.push_back(*last_element);
        point_set.erase(*last_element);

        
         // Remove all points within distance delta from p
        for (auto it = point_set.begin(); it != point_set.end(); ) {

            if (euclidean_distance((*it).second,(*last_element).second) < delta) {
                it = point_set.erase(it);
            } else {
                ++it;
            }
        }
    }
    return S;
}

int main() {


    int N;
    int k;
    double delta;
    int d;
    vector<pair<double,Point>> w; 
 
    cin>>N>>k>>delta>>d;
    for(int i=0;i<N;i++){
        double weight;
        double x,y;
        cin>>weight>>x>>y;
        w.push_back({weight,{x,y}});
    }
    

    auto selected_points = greedy_cdt(w, k, delta);
    
    cout << "Selected points:" << endl;
    for (const auto &q : selected_points) {
        auto p=q.second;
        cout << "(" << p.first << ", " << p.second << ") ->" << q.first<< endl;
    }
    return 0;
}
