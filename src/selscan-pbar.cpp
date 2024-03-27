/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014  Zachary A Szpiech

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/

#include "selscan-pbar.h"
#include <pthread.h>

pthread_mutex_t mutex_progress = PTHREAD_MUTEX_INITIALIZER;

void advanceBar(Bar &bar, double inc)
{
    pthread_mutex_lock(&mutex_progress);
    bar.current += inc;
    if (bar.current / bar.total > double(bar.currentTick) / double(bar.totalTicks))
    {
        bar.currentTick++;
        for (int i = 0; i < bar.totalTicks + 2; i++) cerr << '\b';
        cerr << "|";
        for (int i = 1; i <= bar.totalTicks; i++)
        {
            if (i < bar.currentTick) cerr << "=";
            else if (i == bar.currentTick) cerr << ">";
            else cerr << " ";
        }
        cerr << "|";
        cerr.flush();
    }
    pthread_mutex_unlock(&mutex_progress);
    return;
}

void barInit(Bar &bar, double total, int totalTicks)
{
    bar.total = total;
    bar.current = 0;
    bar.totalTicks = totalTicks;
    bar.currentTick = 0;
    return;
}

#include <iostream>
#include <ctime>
#include <cmath>

double static readTimer() {
    return clock() / (double) CLOCKS_PER_SEC;
}