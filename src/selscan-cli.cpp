/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014-2024  Zachary A Szpiech
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

#include "selscan-cli.h"

//Returns 1 if error
bool initalizeParameters(param_t &params,int argc, char *argv[]){
	params.setPreamble(PREAMBLE);
    params.addFlag(ARG_THREAD, DEFAULT_THREAD, "", HELP_THREAD);
    params.addFlag(ARG_FILENAME_POP1, DEFAULT_FILENAME_POP1, "", HELP_FILENAME_POP1);
    params.addFlag(ARG_FILENAME_POP2, DEFAULT_FILENAME_POP2, "", HELP_FILENAME_POP2);
    params.addFlag(ARG_FILENAME_POP1_TPED, DEFAULT_FILENAME_POP1_TPED, "", HELP_FILENAME_POP1_TPED);
    params.addFlag(ARG_FILENAME_POP2_TPED, DEFAULT_FILENAME_POP2_TPED, "", HELP_FILENAME_POP2_TPED);
    params.addFlag(ARG_FILENAME_POP1_VCF, DEFAULT_FILENAME_POP1_VCF, "", HELP_FILENAME_POP1_VCF);
    params.addFlag(ARG_FILENAME_POP2_VCF, DEFAULT_FILENAME_POP2_VCF, "", HELP_FILENAME_POP2_VCF);
    params.addFlag(ARG_FILENAME_MAP, DEFAULT_FILENAME_MAP, "", HELP_FILENAME_MAP);
    params.addFlag(ARG_PMAP, DEFAULT_PMAP, "", HELP_PMAP);
    params.addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
    params.addFlag(ARG_CUTOFF, DEFAULT_CUTOFF, "", HELP_CUTOFF);
    params.addFlag(ARG_MAX_GAP, DEFAULT_MAX_GAP, "", HELP_MAX_GAP);
    params.addFlag(ARG_GAP_SCALE, DEFAULT_GAP_SCALE, "", HELP_GAP_SCALE);
    params.addFlag(ARG_IHS, DEFAULT_IHS, "", HELP_IHS);
    params.addFlag(ARG_XPNSL, DEFAULT_XPNSL, "", HELP_XPNSL);
    params.addFlag(ARG_UNPHASED, DEFAULT_UNPHASED, "", HELP_UNPHASED);
    params.addFlag(ARG_NSL, DEFAULT_NSL, "", HELP_NSL);
    params.addFlag(ARG_IHS_DETAILED, DEFAULT_IHS_DETAILED, "", HELP_IHS_DETAILED);
    params.addFlag(ARG_SOFT, DEFAULT_SOFT, "", HELP_SOFT);
    params.addFlag(ARG_XP, DEFAULT_XP, "", HELP_XP);
    params.addFlag(ARG_ALT, DEFAULT_ALT, "", HELP_ALT);
    params.addFlag(ARG_MAF, DEFAULT_MAF, "", HELP_MAF);
    params.addFlag(ARG_EHH, DEFAULT_EHH, "", HELP_EHH);
    params.addFlag(ARG_QWIN, DEFAULT_QWIN, "", HELP_QWIN);
    params.addFlag(ARG_SOFT_K, DEFAULT_SOFT_K, "SILENT", HELP_SOFT_K);
    params.addFlag(ARG_MAX_EXTEND, DEFAULT_MAX_EXTEND, "", HELP_MAX_EXTEND);
    params.addFlag(ARG_MAX_EXTEND_NSL, DEFAULT_MAX_EXTEND_NSL, "", HELP_MAX_EXTEND_NSL);
    params.addFlag(ARG_SKIP, DEFAULT_SKIP, "", HELP_SKIP);
    params.addFlag(ARG_KEEP, DEFAULT_KEEP, "", HELP_KEEP);
    params.addFlag(ARG_TRUNC, DEFAULT_TRUNC, "", HELP_TRUNC);
    params.addFlag(ARG_PI, DEFAULT_PI, "", HELP_PI);
    params.addFlag(ARG_PI_WIN, DEFAULT_PI_WIN, "", HELP_PI_WIN);
    params.addFlag(ARG_WAGH, DEFAULT_WAGH, "", HELP_WAGH);

    try
    {
        params.parseCommandLine(argc, argv);
    }
    catch (...)
    {
        return 1;
    }

    return 0;
}
