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
