//function to get random value from 0 to "a"
#include <iostream>

double getrand(int a) {
    double dec,rvalue;
    dec = (rand() % 100)*0.01;
    rvalue = (rand() % a) + dec;
    return rvalue;
}