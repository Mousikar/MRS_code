#include <iostream>
#include <vector>
#include <algorithm>
#define MOD 1000000007

using namespace std;

int countSequences(int n, int m) {
    // 创建一个 DP 表，dp[i][j][k] 表示前 i 个数字异或结果为 k 且最大值为 j 的序列数
    vector<vector<vector<int>>> dp(n + 1, vector<vector<int>>(m + 1, vector<int>(m + 1, 0)));

    // 初始化 base case
    dp[0][0][0] = 1;  // 长度为0,异或结果为0,最大值为0的序列数量为1(空序列)

    // 填充 DP 表
    for (int i = 1; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            for (int k = 0; k <= m; ++k) {
                for (int x = 0; x <= j; ++x) {
                    dp[i][j][k] = (dp[i][j][k] + dp[i-1][x][k ^ j]) % MOD;
                }
            }
        }
    }

    // 计算结果
    int result = 0;
    for (int j = 0; j <= m; ++j) {
        result = (result + dp[n][j][m]) % MOD;
    }

    return result;
}

int main() {
    int n, m;
    cin >> n >> m;
    cout << countSequences(n, m) << endl;
    return 0;
}
