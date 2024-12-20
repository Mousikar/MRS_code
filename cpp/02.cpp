#include <stdio.h>
#include <iostream>
using namespace std;

class B
{
    int a;
    public:
        B(int i=2):a(i){num++;std::cout<<a<<"_"<<num;}
        static int num;
};
int B::num =0;
void f3()
{
    B o1;
    if(o1.num>3)
        throw 1.0;
    else
        throw -1;
}
void f2()
{
    B *o2=new B;
    try
    {
        f3();
    }
    catch(double)
    {
        std::cout << "#" ;
    }
}
void f1()
{
    try
    {
        f2();
        throw 1;
    }
    catch(int)
    {
        std::cout << "&" << endl;
    }
    
}
int main()
{
    B o3(5);
    B &o4=o3;
    f1();
    return 0;
}