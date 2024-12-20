#include <iostream>
#include <string.h>
#include <stdio.h>
using namespace std;
int main()
{
    char a[13] = "Hello world!";
    char *b = "Hello world!";
    cout << sizeof(a) << endl;
    cout << sizeof(b) << endl;
    cout << strlen(a) << endl;
    cout << strlen(b) << endl;
    return 0;
}