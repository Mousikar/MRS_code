// #include <iostream>
// using namespace std;

// int main(){
//     cout<<"Hello, VScode!"<<endl;
//     cout<<"Hello, VScode!"<<endl;
//     cout<<"Hello, VScode!"<<endl;
//     cout<<"Hello, VScode!"<<endl;
//     cout<<"Hello, VScode!"<<endl;
//     cout<<"Hello, VScode!"<<endl;
//     cout<<"Hello, VScode!"<<endl;
//     cout<<"Hello, VScode!"<<endl;
//     return 0;
// }
#include <iostream>
#include <vector>

using namespace std;

int main() {
    int n;
    cin >> n;

    vector<int> heights(n);
    for (int i = 0; i < n; ++i) {
        cin >> heights[i];
    }

    int peaks = 0;

    // 检查左边界
    if (heights[0] > heights[1]) {
        peaks++;
    }

    // 检查中间部分
    for (int i = 1; i < n - 1; ++i) {
        if (heights[i] > heights[i - 1] && heights[i] > heights[i + 1]) {
            peaks++;
        }
    }

    // 检查右边界
    if (heights[n - 1] > heights[n - 2]) {
        peaks++;
    }

    cout << peaks << endl;

    return 0;
}
