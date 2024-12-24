#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

class Solution {
public:
    int minRefuelStops(int target, int startFuel, vector<vector<int>>& stations) {
        stations.push_back({target, 0});
        priority_queue<int> pq;
        int fuel = startFuel;
        int stops = 0;
        int prevPosition = 0;

        for (const auto& station : stations) {
            int position = station[0];
            int fuelAtStation = station[1];
            fuel -= (position - prevPosition);
            while (fuel < 0) {  
                if (pq.empty()) return -1;
                fuel += pq.top();
                pq.pop();
                stops++;
            }
            pq.push(fuelAtStation);
            prevPosition = position;
        }
        return stops;
    }
};
