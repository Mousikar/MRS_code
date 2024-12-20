// ------------------------------------------------------------------------------
// #include <iostream>
// #include <stack>
// #include <string>
// using namespace std;

// bool isValid(std::string s) {
//     std::stack<char> stack;

//     for (char c : s) {
//         if (c == '(' || c == '{' || c == '[') {
//             stack.push(c);
//         } else {
//             if (stack.empty()) return false;
//             char top = stack.top();
//             stack.pop();
//             if ((c == ')' && top != '(') || 
//                 (c == '}' && top != '{') || 
//                 (c == ']' && top != '[')) {
//                 return false;
//             }
//         }
//     }

//     return stack.empty();
// }

// int main() {
//     string s;
//     cin>>s;
//     std::cout << isValid(s);

//     return 0;
// }
// ------------------------------------------------------------------------------
// 示例 1:
// 输入: "()"
// 输出: true
// 示例 2:
// 输入: "()[]{}"
// 输出: true
// 示例 3:
// 输入: "(]"
// 输出: false
// 示例 4:
// 输入: "([)]"
// 输出: false
// 示例 5:
// 输入: "{[]}"
// 输出: true