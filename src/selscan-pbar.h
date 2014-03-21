#ifndef __XP_IHH_PBAR_H__
#define __XP_IHH_PBAR_H__
#include <iostream>

using namespace std;

struct Bar
{
    double total;
    double current;
    int totalTicks;
    int currentTick;
};

void advanceBar(Bar &bar, double inc);
void barInit(Bar &bar, double total, int totalTicks);

#endif
