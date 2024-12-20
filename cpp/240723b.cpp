#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace std;

int min_operations(int n, vector<int>& a) {
    // 目标排列是 [1, 2, ..., n]
    vector<int> target(n);
    for (int i = 0; i < n; ++i) {
        target[i] = i + 1;
    }
    
    // 找到当前排列与目标排列的位置偏移
    unordered_map<int, int> index_map;
    for (int i = 0; i < n; ++i) {
        index_map[a[i]] = i;
    }
    
    vector<int> distances(n);
    for (int i = 0; i < n; ++i) {
        distances[i] = index_map[target[i]] - i;
    }
    
    // 计算左移和右移的最小操作次数
    int left_moves = 0, right_moves = 0;
    for (int i = 0; i < n; ++i) {
        left_moves += (distances[i] + n) % n;
        right_moves += (-distances[i] + n) % n;
    }
    
    return min(left_moves, right_moves);
}

int main() {
    int n;
    cin >> n;
    vector<int> a(n);
    for (int i = 0; i < n; ++i) {
        cin >> a[i];
    }
    
    cout << min_operations(n, a) << endl;
    
    return 0;
}
