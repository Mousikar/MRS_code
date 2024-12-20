#include <iostream>
#include <deque>
#include <vector>
#include <algorithm>

int min_operations(int n, std::vector<int>& a) {
    std::deque<int> dq(a.begin(), a.end());
    int operations = 0;

    // 创建目标序列
    std::vector<int> target(n);
    for (int i = 0; i < n; ++i) {
        target[i] = i + 1;
    }

    while (dq != std::deque<int>(target.begin(), target.end())) {
        // 尝试从左边进行操作
        std::deque<int> dq_left = dq;
        int left_operations = 0;
        while (dq_left != std::deque<int>(target.begin(), target.end()) && left_operations <= n) {
            dq_left.pop_front();
            dq_left.push_back(0);
            left_operations++;
        }

        // 尝试从右边进行操作
        std::deque<int> dq_right = dq;
        int right_operations = 0;
        while (dq_right != std::deque<int>(target.begin(), target.end()) && right_operations <= n) {
            dq_right.pop_back();
            dq_right.push_front(0);
            right_operations++;
        }

        // 更新操作次数
        operations += std::min(left_operations, right_operations);

        // 更新队列状态
        if (left_operations <= right_operations) {
            dq = dq_left;
        } else {
            dq = dq_right;
        }
    }

    return operations;
}

int main() {
    int n;
    std::cin >> n;
    std::vector<int> a(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> a[i];
    }

    int result = min_operations(n, a);
    std::cout << result << std::endl;

    return 0;
}
