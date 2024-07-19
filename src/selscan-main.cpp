//selscan-main.cpp
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
#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <pthread.h>
#include <map>
#include "selscan-cli.h"
#include "selscan-maintools.h"
#include "selscan-data.h"
#include "selscan-pbar.h"
#include "param_t.h"

#define IS_TRUE(x) { if (!(x)) std::cout << __FUNCTION__ << " failed on line " << __LINE__ << std::endl; }

using namespace std;

int main(int argc, char *argv[])
{
    cout<<"Max threads supoorted: "<<std::thread::hardware_concurrency()<<endl;
    double time_start = (clock() / (double) CLOCKS_PER_SEC);

    auto start = std::chrono::high_resolution_clock::now();

    cerr << "selscan v" + VERSION + "\n";
// #ifdef PTW32_STATIC_LIB
//     pthread_win32_process_attach_np();
// #endif

    param_t params;
    bool ERROR = initalizeParameters(params,argc,argv);
    if (ERROR) return 1;
    ERROR = checkParameters(params,argc,argv);
    if (ERROR) return 1;

    MainTools mt(params); // pass command line parameters to MainTools
    mt.run(argc, argv); // run the main program

// #ifdef PTW32_STATIC_LIB
//     pthread_win32_process_detach_np();
// #endif

    return 0;
}
