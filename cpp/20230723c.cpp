#include <iostream>
#include <vector>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <algorithm> // 包含algorithm头文件以使用reverse
#include <tuple>
#include <sstream>

using namespace std;

// 将状态转换为字符串以便存储在集合和字典中
string state_to_string(const vector<int>& state) {
    string result;
    for (int num : state) {
        result += to_string(num) + ",";
    }
    return result;
}

// 回溯路径以构造操作步骤
vector<string> construct_path(unordered_map<string, pair<string, string>>& parents, string& state) {
    vector<string> path;
    while (parents[state].first != "") {
        path.push_back(parents[state].second);
        state = parents[state].first;
    }
    reverse(path.begin(), path.end());
    return path;
}

// 计算最小操作次数并记录每一步操作
pair<int, vector<string>> min_operations_to_target(int n, vector<int>& a) {
    vector<int> target(n);
    for (int i = 0; i < n; ++i) {
        target[i] = i;
    }
    string target_state = state_to_string(target);

    string initial_state = state_to_string(a);
    if (initial_state == target_state) {
        return {0, {}};
    }

    queue<pair<vector<int>, int>> q;
    q.push({a, 0});
    unordered_set<string> visited;
    visited.insert(initial_state);

    unordered_map<string, pair<string, string>> parents;
    parents[initial_state] = {"", ""};

    while (!q.empty()) {
        auto [current_state, steps] = q.front();
        q.pop();
        string current_str = state_to_string(current_state);

        // 操作1: 丢弃最高位，整体左移，最低位补充一个合适的数字
        for (int new_digit = 0; new_digit < n; ++new_digit) {
            vector<int> new_state(current_state.begin() + 1, current_state.end());
            new_state.push_back(new_digit);
            string new_state_str = state_to_string(new_state);
            if (new_state_str == target_state) {
                parents[new_state_str] = {current_str, "丢弃最高位，补充" + to_string(new_digit)};
                return {steps + 1, construct_path(parents, new_state_str)};
            }
            if (visited.find(new_state_str) == visited.end()) {
                visited.insert(new_state_str);
                q.push({new_state, steps + 1});
                parents[new_state_str] = {current_str, "丢弃最高位，补充" + to_string(new_digit)};
            }
        }

        // 操作2: 丢弃最低位，整体右移，最高位补充一个合适的数字
        for (int new_digit = 0; new_digit < n; ++new_digit) {
            vector<int> new_state = {new_digit};
            new_state.insert(new_state.end(), current_state.begin(), current_state.end() - 1);
            string new_state_str = state_to_string(new_state);
            if (new_state_str == target_state) {
                parents[new_state_str] = {current_str, "丢弃最低位，补充" + to_string(new_digit)};
                return {steps + 1, construct_path(parents, new_state_str)};
            }
            if (visited.find(new_state_str) == visited.end()) {
                visited.insert(new_state_str);
                q.push({new_state, steps + 1});
                parents[new_state_str] = {current_str, "丢弃最低位，补充" + to_string(new_digit)};
            }
        }
    }

    return {-1, {}};
}

int main() {
    int n;
    cin >> n;
    cin.ignore(); // 忽略换行符

    string line;
    getline(cin, line);
    stringstream ss(line);
    vector<int> a(n);
    for (int i = 0; i < n; ++i) {
        ss >> a[i];
        a[i]--; // 每个数字减去1
    }

    auto [result, steps] = min_operations_to_target(n, a);

    cout << result << endl;
    for (const string& step : steps) {
        cout << step << endl;
    }

    return 0;
}
