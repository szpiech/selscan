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

#include <cmath> //isnan and isinf
#include "selscan-maintools.h"
pthread_mutex_t mutex_log = PTHREAD_MUTEX_INITIALIZER;

void calculatePi(HaplotypeData *hapData, MapData *mapData, int winsize, string outFilename)
{
    ofstream fout;
    fout.open(outFilename.c_str());
    if (fout.fail())
    {
        cerr << "ERROR: Failed to open " << outFilename << " for writing.\n";
        throw 1;
    }

    int start = 1;
    //int end = mapData->physicalPos[mapData->nloci - 1];
    int startLocus = 0;
    int endLocus = 0;
    //Identify the start and end indicies in the window
    int winStart = start;
    int winEnd = winStart + winsize - 1;
    int pos;
    int length = 0;
    double pi = 0;
    double denominator = (hapData->nhaps) * (hapData->nhaps - 1) * 0.5;
    int locus;
    for (locus = 0; locus < mapData->nloci; locus++)
    {
        pos = mapData->physicalPos[locus];
        while (pos > winEnd)
        {
            endLocus = locus - 1;
            length = endLocus - startLocus + 1;
            fout << winStart << " " << winEnd << " ";
            //do stuff

            if (length == 0)
            {
                pi = 0;
            }
            else
            {
                for (int i = 0; i < hapData->nhaps; i++)
                {
                    for (int j = i + 1; j < hapData->nhaps; j++)
                    {
                        pi += hamming_dist_ptr(hapData->data[i] + startLocus, hapData->data[j] + startLocus, length);
                    }
                }
            }

            fout << pi / denominator << endl;

            winStart += winsize;
            winEnd += winsize;

            startLocus = locus;
            pi = 0;
        }
    }
    //final window
    endLocus = locus - 1;
    length = endLocus - startLocus + 1;
    fout << winStart << " " << winEnd << " ";
    //do stuff

    if (length == 0)
    {
        pi = 0;
    }
    else
    {
        for (int i = 0; i < hapData->nhaps; i++)
        {
            for (int j = i + 1; j < hapData->nhaps; j++)
            {
                pi += hamming_dist_ptr(hapData->data[i] + startLocus, hapData->data[j] + startLocus, length);
            }
        }
    }

    fout << pi / denominator << endl;



    fout.close();
    return;
}

int queryFound(MapData *mapData, string query)
{
    for (int locus = 0; locus < mapData->nloci; locus++)
    {
        if (mapData->locusName[locus].compare(query) == 0) return locus;
        //This is only set when compiling windows binaries
        //And I needed to use itoa for some reason under windows
        #ifdef PTW32_STATIC_LIB
        char buffer[100];
        itoa(mapData->physicalPos[locus], buffer, 10);
        if (string(buffer).compare(query) == 0) return locus;
        #else
        if (to_string(mapData->physicalPos[locus]).compare(query) == 0) return locus;
        #endif
    }

    return -1;
}

double calcFreq(HaplotypeData *hapData, int locus, bool unphased)
{
    double total = 0;
    double freq = 0;

    for (int hap = 0; hap < hapData->nhaps; hap++)
    {
        if (hapData->data[hap][locus] != MISSING_CHAR)
        {
            if (hapData->data[hap][locus] == '1') freq += 1;
            else if (hapData->data[hap][locus] == '2') freq += 2;
            
            if (unphased) total+=2;
            else total++;
        }
    }
    return (freq / total);
}
double getMaxFreq(map<string, int> &count){
    map<string, int>::iterator i;
    double total = 0;
    double max = 0;

    for (i = count.begin(); i != count.end(); i++){
        total += i->second;
        if (i->second > max) max = i->second;
    }
    return max/total;
}

double getMaxFreq(double a, double b, double c){
    double total = a+b+c;
    double max = 0;

    if (a > b) max = a;
    if (c > a) max = c

    return max/total;
}

double scaleHomozygosity(double h, double M)
{
    if(M == 1) return 1;
    if(M == 0) return 0;

    double ceil_invM = ceil(1/M);
    double lowerBound = M*M;
    double upperBound = 1 - M*(ceil_invM - 1)*(2 - ceil_invM*M);
    
    return (h - lowerBound)/upperBound;
}

void query_locus(void *order)
{
    work_order_t *p = (work_order_t *)order;
    char **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nhaps = p->hapData->nhaps;
    int *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    //string *locusName = p->mapData->locusName;
    //ofstream *flog = p->flog;
    ofstream *fout = p->fout;
    string outFilename = p->filename;
    //int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    //int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    //double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    bool ALT = p->params->getBoolFlag(ARG_ALT);
    //bool WAGH = p->params->getBoolFlag(ARG_WAGH);
    double (*calc)(map<string, int> &, int, bool) = p->calc;
    bool unphased = p->params->getBoolFlag(ARG_UNPHASED);
    bool SCALE_HOM = p->params->getBoolFlag(ARG_RESCALE);

    int locus = p->queryLoc;
    int queryPad = p->params->getIntFlag(ARG_QWIN);
    int stopLeft = locus;
    for (int i = locus - 1; i >= 0; i--)
    {
        if (physicalPos[locus] - physicalPos[i] <= queryPad) stopLeft = i;
    }
    int stopRight = locus;
    for (int i = locus + 1; i < nloci; i++)
    {
        if (physicalPos[i] - physicalPos[locus] <= queryPad) stopRight = i;
    }


    //cerr << stopLeft << " " << physicalPos[stopLeft] << " " << stopRight << " " << physicalPos[stopRight] << endl;

    //EHH to the left of the core snp
    double current_derived_ehh = 1;
    double current_ancestral_ehh = 1;
    double current_ehh = 1;
    double derivedCount = 0;

    //used for unphased analyses
    double hetCount = 0;
    double ancestralCount = 0;
    double current_notDerived_ehh = 1;
    double current_notAncestral_ehh = 1;
    
    //A list of all the haplotypes
    //Starts with just the focal snp and grows outward
    string *haplotypeList = new string[nhaps];
    for (int hap = 0; hap < nhaps; hap++)
    {
        haplotypeList[hap] = data[hap][locus];
        if (unphased){
            derivedCount += ( data[hap][locus] == '2' ) ? 1 : 0;
            ancestralCount += ( data[hap][locus] == '0' ) ? 1 : 0;
            hetCount += ( data[hap][locus] == '1' ) ? 1 : 0;
        }
        else{
            derivedCount += ( data[hap][locus] == '1' ) ? 1 : 0;
            ancestralCount += ( data[hap][locus] == '0' ) ? 1 : 0;
        }
        
    }

    if (ALT){
        double h = derivedCount/nhaps;
        current_ehh = h;
        h = ancestralCount/nhaps;
        current_ehh += h;
        h = hetCount/nhaps;
        current_ehh += h;
    }
    else{
        current_ehh = (derivedCount > 1) ? (nCk(derivedCount, 2) / nCk(nhaps, 2)) : 0;
        current_ehh += (ancestralCount > 1) ? (nCk(ancestralCount, 2) / nCk(nhaps, 2)) : 0;
        current_ehh += (hetCount > 1) ? (nCk(hetCount, 2) / nCk(nhaps, 2)) : 0;
    }

    if(SCALE_HOM) current_ehh = scaleHomozygosity(current_ehh,getMaxFreq(derivedCount,ancestralCount,hetCount));
    

    
    /*
        if (derivedCount == 0 || derivedCount == nhaps)
        {
            cerr << "ERROR: " << locusName[locus] << " is monomorphic.\n";
            (*fout) << "ERROR: " << locusName[locus] << " is monomorphic.\n";
            return;
        }
    */
    //cerr << "numHaps: " << nhaps << "\nderivedCounts: " << derivedCount << endl;
    /*
    int **ancestralHapColor = new int *[int(nhaps - derivedCount)];
    for (int i = 0; i < nhaps - derivedCount; i++)
    {
        ancestralHapColor[i] = new int[stopRight - stopLeft + 1];
        ancestralHapColor[i][locus - stopLeft] = 0;
    }
    int **derivedHapColor = new int *[int(derivedCount)];
    for (int i = 0; i < derivedCount; i++)
    {
        derivedHapColor[i] = new int[stopRight - stopLeft + 1];
        derivedHapColor[i][locus - stopLeft] = 0;
    }
    */
    //cerr << "allocated hap color arrays.\n";

    bool isDerived;
    string hapStr;
    string *tempResults = new string[locus - stopLeft];
    int tempIndex = locus - stopLeft - 1;
    //int derivedCurrentColor = 0;
    //int ancestralCurrentColor = 0;

    for (int i = locus - 1; i >= stopLeft; i--)
    {
        int numDerived = 0;
        int numAncestral = 0;
        map<string, int> ancestralHapCount;
        map<string, int> derivedHapCount;
        map<string, int> hapCount;
        
        //for unphased
        int numHet = 0;
        map<string, int> notAncestralHapCount;
        map<string, int> notDerivedHapCount;

        for (int hap = 0; hap < nhaps; hap++)
        {
            if(unphased){
                haplotypeList[hap] += data[hap][i];
                hapStr = haplotypeList[hap];

                if (data[hap][locus] == '0'){
                    //count ancestral haplotype
                    if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                    else ancestralHapCount[hapStr]++;
                    numAncestral++;
                    //count non-derived hapoltype
                    if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                    else notDerivedHapCount[hapStr]++;
                }
                else if (data[hap][locus] == '1'){
                    //count non-derived hapoltype
                    if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                    else notDerivedHapCount[hapStr]++;
                    //count non-ancestral haplotype
                    if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                    else notAncestralHapCount[hapStr]++;
                    numHet++;
                }
                else{//data[hap][locus] == '2'
                    //count derived hapoltype
                    if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                    else derivedHapCount[hapStr]++;
                    numDerived++;
                    //count non-ancestral haplotype
                    if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                    else notAncestralHapCount[hapStr]++;
                }
            }
            else{
                isDerived = ( data[hap][locus] == '1') ? 1 : 0;
                //build haplotype string
                //char digit[2];
                //snprintf(digit, "%d", data[hap][i]);
                haplotypeList[hap] += data[hap][i];
                hapStr = haplotypeList[hap];

                if (isDerived)
                {
                    //count derived hapoltype
                    if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                    else derivedHapCount[hapStr]++;
                    numDerived++;
                }
                else
                {
                    //count ancestral haplotype
                    if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                    else ancestralHapCount[hapStr]++;
                    numAncestral++;
                }

                if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
                else hapCount[hapStr]++;
            }
        }
        //Write functions to fill in haplotype colors here
        //fillColors(derivedHapColor, derivedHapCount, haplotypeList, nhaps, tempIndex, derivedCurrentColor, true);
        //fillColors(ancestralHapColor, ancestralHapCount, haplotypeList, nhaps, tempIndex, ancestralCurrentColor, true);

        current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);
        current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);
        current_ehh = (*calc)(hapCount, numAncestral + numDerived + numHet, ALT);

        if(unphased){
            current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);
            current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);
        }

        if (SCALE_HOM){
            current_derived_ehh = scaleHomozygosity(current_derived_ehh,getMaxFreq(derivedHapCount));
            current_ancestral_ehh = scaleHomozygosity(current_ancestral_ehh,getMaxFreq(ancestralHapCount));
            current_ehh = scaleHomozygosity(current_ehh,getMaxFreq(hapCount));
            if(unphased){
                current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);
                current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);
            }
        }

    // use stringstreams  
        char tempStr[200];
        if(unphased){
            snprintf(tempStr,200, "%d\t%f\t%f\t%f\t%f\t%f\t%f", physicalPos[i] - physicalPos[locus], geneticPos[i] - geneticPos[locus], current_derived_ehh, current_ancestral_ehh, current_notDerived_ehh, current_notAncestral_ehh, current_ehh);
        }
        else{
            snprintf(tempStr,200, "%d\t%f\t%f\t%f\t%f", physicalPos[i] - physicalPos[locus], geneticPos[i] - geneticPos[locus], current_derived_ehh, current_ancestral_ehh, current_ehh);
        }
        tempResults[tempIndex] = string(tempStr);
        tempIndex--;
    }

    delete [] haplotypeList;

    (*fout) << "pdist\tgdist\tderEHH\tancEHH\t";
    if(unphased){
        (*fout) << "notAncEHH\t"
                << "notDerEHH\t";
    }
    (*fout) << "EHH" << endl;

    for (int i = 0; i < locus - stopLeft; i++)
    {
        (*fout) << tempResults[i] << "\n";
    }
    delete [] tempResults;

    //calculate EHH to the right
    current_derived_ehh = 1;
    current_ancestral_ehh = 1;

    if (ALT){
        double h = derivedCount/nhaps;
        current_ehh = h;
        h = ancestralCount/nhaps;
        current_ehh += h;
        h = hetCount/nhaps;
        current_ehh += h;
    }
    else{
        current_ehh = (derivedCount > 1) ? (nCk(derivedCount, 2) / nCk(nhaps, 2)) : 0;
        current_ehh += (ancestralCount > 1) ? (nCk(ancestralCount, 2) / nCk(nhaps, 2)) : 0;
        current_ehh += (hetCount > 1) ? (nCk(hetCount, 2) / nCk(nhaps, 2)) : 0;
    }

    fout->precision(6);
    (*fout) << std::fixed <<  physicalPos[locus] - physicalPos[locus]  << "\t"
            << geneticPos[locus] - geneticPos[locus] << "\t"
            << current_derived_ehh << "\t"
            << current_ancestral_ehh << "\t";
    if(unphased){
        (*fout) << current_notAncestral_ehh << "\t"
                << current_notDerived_ehh << "\t";
    }
    (*fout) << current_ehh << endl;

    //A list of all the haplotypes
    //Starts with just the focal snp and grows outward
    haplotypeList = new string[nhaps];
    for (int hap = 0; hap < nhaps; hap++)
    {
        //char digit[2];
        //snprintf(digit, "%d", data[hap][locus]);
        haplotypeList[hap] = data[hap][locus];
    }

    //derivedCurrentColor = 0;
    //ancestralCurrentColor = 0;

    //while(current_ancestral_ehh > EHH_CUTOFF || current_derived_ehh > EHH_CUTOFF)
    for (int i = locus + 1; i <= stopRight; i++)
    {
        int numDerived = 0;
        int numAncestral = 0;
        map<string, int> ancestralHapCount;
        map<string, int> derivedHapCount;
        map<string, int> hapCount;
        
        //for unphased
        int numHet = 0;
        map<string, int> notAncestralHapCount;
        map<string, int> notDerivedHapCount;

        for (int hap = 0; hap < nhaps; hap++)
        {
            if(unphased){
                haplotypeList[hap] += data[hap][i];
                hapStr = haplotypeList[hap];

                if (data[hap][locus] == '0'){
                    //count ancestral haplotype
                    if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                    else ancestralHapCount[hapStr]++;
                    numAncestral++;
                    //count non-derived hapoltype
                    if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                    else notDerivedHapCount[hapStr]++;
                }
                else if (data[hap][locus] == '1'){
                    //count non-derived hapoltype
                    if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                    else notDerivedHapCount[hapStr]++;
                    //count non-ancestral haplotype
                    if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                    else notAncestralHapCount[hapStr]++;
                    numHet++;
                }
                else{//data[hap][locus] == '2'
                    //count derived hapoltype
                    if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                    else derivedHapCount[hapStr]++;
                    numDerived++;
                    //count non-ancestral haplotype
                    if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                    else notAncestralHapCount[hapStr]++;
                }
            }
            else{
                isDerived = ( data[hap][locus] == '1' ) ? 1 : 0;
                //build haplotype string
                //char digit[2];
                //snprintf(digit, "%d", data[hap][i]);
                haplotypeList[hap] += data[hap][i];
                hapStr = haplotypeList[hap];

                if (isDerived)
                {
                    //count hapoltypes
                    if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                    else derivedHapCount[hapStr]++;
                    numDerived++;
                }
                else //ancestral
                {
                    if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                    else ancestralHapCount[hapStr]++;
                    numAncestral++;
                }

                if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
                else hapCount[hapStr]++;
            }
        }

        //Write functions to fill in haplotype colors here
        //fillColors(derivedHapColor, derivedHapCount, haplotypeList, nhaps, i - stopLeft, derivedCurrentColor, false);
        //fillColors(ancestralHapColor, ancestralHapCount, haplotypeList, nhaps, i - stopLeft, ancestralCurrentColor, false);

        current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);
        current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);
        current_ehh = (*calc)(hapCount, numAncestral + numDerived + numHet, ALT);

        if(unphased){
            current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);
            current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);
        }

        if (SCALE_HOM){
            current_derived_ehh = scaleHomozygosity(current_derived_ehh,getMaxFreq(derivedHapCount));
            current_ancestral_ehh = scaleHomozygosity(current_ancestral_ehh,getMaxFreq(ancestralHapCount));
            current_ehh = scaleHomozygosity(current_ehh,getMaxFreq(hapCount));
            if(unphased){
                current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);
                current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);
            }
        }

        (*fout) << std::fixed <<  physicalPos[i] - physicalPos[locus]  << "\t"
            << geneticPos[i] - geneticPos[locus] << "\t"
            << current_derived_ehh << "\t"
            << current_ancestral_ehh << "\t";
        if(unphased){
            (*fout) << current_notAncestral_ehh << "\t"
                    << current_notDerived_ehh << "\t";
        }
        (*fout) << current_ehh << endl;
    }

    delete [] haplotypeList;
    /*
    ofstream fout2;
    string outFilenameDer = outFilename + ".der.colormap";
    fout2.open(outFilenameDer.c_str());
    if (fout2.fail())
    {
        cerr << "ERROR: could not open " << outFilenameDer << " for writing.\n";
        throw 1;
    }

    for (int i = 0; i < derivedCount; i++)
    {
        for (int j = 0; j < stopRight - stopLeft + 1; j++)
        {
            fout2 << derivedHapColor[i][j] << " ";
        }
        fout2 << endl;
    }
    fout2.close();

    string outFilenameAnc = outFilename + ".anc.colormap";
    fout2.open(outFilenameAnc.c_str());
    if (fout2.fail())
    {
        cerr << "ERROR: could not open " << outFilenameAnc << " for writing.\n";
        throw 1;
    }

    for (int i = 0; i < nhaps - derivedCount; i++)
    {
        for (int j = 0; j < stopRight - stopLeft + 1; j++)
        {
            fout2 << ancestralHapColor[i][j] << " ";
        }
        fout2 << endl;
    }
    fout2.close();
    */
    return;
}

void query_locus_soft(void *order)
{
    work_order_t *p = (work_order_t *)order;
    char **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nhaps = p->hapData->nhaps;
    int *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    //string *locusName = p->mapData->locusName;

    //ofstream *flog = p->flog;
    ofstream *fout = p->fout;
    string outFilename = p->filename;
    //int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    //int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    //double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    //bool ALT = p->params->getBoolFlag(ARG_ALT);
    //bool WAGH = p->params->getBoolFlag(ARG_WAGH);

    int locus = p->queryLoc;
    int queryPad = p->params->getIntFlag(ARG_QWIN);
    int stopLeft = locus;
    for (int i = locus - 1; i >= 0; i--)
    {
        if (physicalPos[locus] - physicalPos[i] <= queryPad) stopLeft = i;
    }
    int stopRight = locus;
    for (int i = locus + 1; i < nloci; i++)
    {
        if (physicalPos[i] - physicalPos[locus] <= queryPad) stopRight = i;
    }

    //EHH to the left of the core snp

    double current_ehh1 = 1;
    double current_ehh2d1 = 1;
    double current_ehh12 = 1;

    double previous_ehh1 = 1;
    double previous_ehh2d1 = 1;
    double previous_ehh12 = 1;

    double derivedCount = 0;
    //A list of all the haplotypes
    //Starts with just the focal snp and grows outward
    map<string, int> tempHapCount;
    string *haplotypeList = new string[nhaps];
    for (int hap = 0; hap < nhaps; hap++)
    {
        derivedCount += data[hap][locus];
        char digit[2];
        snprintf(digit,2, "%d", data[hap][locus]);
        haplotypeList[hap] = digit;
        string hapStr = haplotypeList[hap];
        //count hapoltype freqs
        if (tempHapCount.count(hapStr) == 0) tempHapCount[hapStr] = 1;
        else tempHapCount[hapStr]++;
    }

    triplet_t res = calculateSoft(tempHapCount, nhaps);
    current_ehh1 = res.h1;
    current_ehh2d1 = res.h2dh1;
    current_ehh12 = res.h12;
    previous_ehh1 = res.h1;
    previous_ehh2d1 = res.h2dh1;
    previous_ehh12 = res.h12;

    //cerr << "numHaps: " << nhaps << "\nderivedCounts: " << derivedCount << endl;
    /*
    int **ancestralHapColor = new int*[int(nhaps-derivedCount)];
    for(int i = 0; i < nhaps-derivedCount; i++)
    {
      ancestralHapColor[i] = new int[stopRight-stopLeft+1];
      ancestralHapColor[i][locus-stopLeft] = 0;
    }
    int **derivedHapColor = new int*[int(derivedCount)];
    for(int i = 0; i < derivedCount; i++)
    {
      derivedHapColor[i] = new int[stopRight-stopLeft+1];
      derivedHapColor[i][locus-stopLeft] = 0;
    }
    */
    //cerr << "allocated hap color arrays.\n";

    string *tempResults = new string[locus - stopLeft];
    int tempIndex = locus - stopLeft - 1;
    //int derivedCurrentColor = 0;
    //int ancestralCurrentColor = 0;

    for (int i = locus - 1; i >= stopLeft; i--)
    {
        //int numDerived = 0;
        //int numAncestral = 0;
        map<string, int> hapCount;

        for (int hap = 0; hap < nhaps; hap++)
        {
            //build haplotype string
            char digit[2];
            snprintf(digit,2, "%d", data[hap][i]);
            haplotypeList[hap] += digit;
            string hapStr = haplotypeList[hap];

            //count hapoltype freqs
            if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
            else hapCount[hapStr]++;
        }

        //We've now counted all of the unique haplotypes extending out of the core SNP
        res = calculateSoft(hapCount, nhaps);
        current_ehh1 = res.h1;
        current_ehh2d1 = res.h2dh1;
        current_ehh12 = res.h12;
        //Write functions to fill in haplotype colors here
        /*
        fillColors(derivedHapColor, derivedHapCount,haplotypeList, nhaps,tempIndex,derivedCurrentColor,true);
        fillColors(ancestralHapColor, ancestralHapCount,haplotypeList, nhaps,tempIndex,ancestralCurrentColor,true);
        */

        char tempStr[100];
        snprintf(tempStr,100, "%d\t%f\t%f\t%f\t%f", physicalPos[i] - physicalPos[locus], geneticPos[i] - geneticPos[locus], current_ehh1, current_ehh12, current_ehh2d1);
        tempResults[tempIndex] = string(tempStr);
        tempIndex--;
    }

    delete [] haplotypeList;

    for (int i = 0; i < locus - stopLeft; i++)
    {
        (*fout) << tempResults[i] << "\n";
    }
    delete [] tempResults;

    //calculate EHH to the right
    current_ehh1 = 1;
    current_ehh2d1 = 1;
    current_ehh12 = 1;

    previous_ehh1 = 1;
    previous_ehh2d1 = 1;
    previous_ehh12 = 1;

    //A list of all the haplotypes
    //Starts with just the focal snp and grows outward
    tempHapCount.clear();
    haplotypeList = new string[nhaps];
    for (int hap = 0; hap < nhaps; hap++)
    {
        char digit[2];
        snprintf(digit,2, "%d", data[hap][locus]);
        haplotypeList[hap] = digit;
        string hapStr = haplotypeList[hap];
        //count hapoltype freqs
        if (tempHapCount.count(hapStr) == 0) tempHapCount[hapStr] = 1;
        else tempHapCount[hapStr]++;
    }

    res = calculateSoft(tempHapCount, nhaps);
    current_ehh1 = res.h1;
    current_ehh2d1 = res.h2dh1;
    current_ehh12 = res.h12;
    previous_ehh1 = res.h1;
    previous_ehh2d1 = res.h2dh1;
    previous_ehh12 = res.h12;

    fout->precision(6);
    (*fout) << std::fixed <<  physicalPos[locus] - physicalPos[locus]  << "\t"
            << geneticPos[locus] - geneticPos[locus] << "\t"
            << current_ehh1 << "\t"
            << current_ehh12 << "\t"
            << current_ehh2d1 << endl;

    //while(current_ancestral_ehh > EHH_CUTOFF || current_derived_ehh > EHH_CUTOFF)
    for (int i = locus + 1; i <= stopRight; i++)
    {
        map<string, int> hapCount;

        for (int hap = 0; hap < nhaps; hap++)
        {
            //build haplotype string
            char digit[2];
            snprintf(digit,2, "%d", data[hap][i]);
            haplotypeList[hap] += digit;
            string hapStr = haplotypeList[hap];

            //count hapoltype freqs
            if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
            else hapCount[hapStr]++;
        }

        //We've now counted all of the unique haplotypes extending out of the core SNP
        res = calculateSoft(hapCount, nhaps);
        current_ehh1 = res.h1;
        current_ehh2d1 = res.h2dh1;
        current_ehh12 = res.h12;

        //Write functions to fill in haplotype colors here
        /*
        fillColors(derivedHapColor, derivedHapCount,haplotypeList, nhaps,i-stopLeft,derivedCurrentColor,false);
        fillColors(ancestralHapColor, ancestralHapCount,haplotypeList, nhaps,i-stopLeft,ancestralCurrentColor,false);
        */

        (*fout) << physicalPos[i] - physicalPos[locus] << "\t"
                << geneticPos[i] - geneticPos[locus] << "\t"
                << current_ehh1 << "\t"
                << current_ehh12 << "\t"
                << current_ehh2d1 << endl;
    }

    delete [] haplotypeList;

    /*
    ofstream fout2;
    string outFilenameDer = outFilename + ".der.colormap";
    fout2.open(outFilenameDer.c_str());
    if(fout2.fail())
    {
      cerr << "ERROR: could not open " << outFilenameDer << " for writing.\n";
      throw 1;
    }
    for(int i = 0; i < derivedCount; i++)
    {
      for(int j = 0; j < stopRight-stopLeft+1; j++)
      {
        fout2 << derivedHapColor[i][j] << " ";
        }
      fout2 << endl;
    }
    fout2.close();
    string outFilenameAnc = outFilename + ".anc.colormap";
    fout2.open(outFilenameAnc.c_str());
    if(fout2.fail())
    {
      cerr << "ERROR: could not open " << outFilenameAnc << " for writing.\n";
      throw 1;
    }
    for(int i = 0; i < nhaps-derivedCount; i++)
    {
      for(int j = 0; j < stopRight-stopLeft+1; j++)
      {
        fout2 << ancestralHapColor[i][j] << " ";
        }
      fout2 << endl;
    }
    fout2.close();
    */
    return;
}

void fillColors(int **hapColor, map<string, int> &hapCount, string *haplotypeList, int hapListLength, int currentLoc, int &currentColor, bool left)
{
    map<string, int>::iterator i;
    //int numUniqueHaps = hapCount.size();
    //int mostCommonHapCount = 0;
    int nhaps = 0;

    for (i = hapCount.begin(); i != hapCount.end(); i++)
    {
        nhaps += i->second;
    }
    //holds colors for haplotypes that have already been seen in the search
    map<string, int> hapSeen;

    string mostCommonHap = "_NONE_";

    int colorIndex = 0;
    int previousLoc = (left) ? currentLoc + 1 : currentLoc - 1;
    for (int j = 0; j < hapListLength; j++)
    {
        string hapStr = haplotypeList[j];
        if (hapCount.count(hapStr) == 0) continue;

        /*
        for(int h = 0; h < hapListLength; h++)
        {
          if(h == j) cerr << ">>";
          else cerr << "  ";
          cerr << haplotypeList[h] << endl;
        }
        */

        if (hapCount[hapStr] == 1)
        {
            //cerr << "Unique hap\n";
            hapColor[colorIndex][currentLoc] = -1;
            colorIndex++;
            continue;
        }

        //If there was a split in the current haplotype family
        //let the most common continued haplotype keep the color
        //then increment the color for the less common one
        if (familyDidSplit(hapStr, hapCount[hapStr], hapColor, nhaps, colorIndex, previousLoc, mostCommonHap))
        {
            //cerr << "split\n";
            //cerr << hapStr << " " << mostCommonHap << endl;
            if (hapStr.compare(mostCommonHap) == 0)
            {
                //cerr << "common\n";
                hapColor[colorIndex][currentLoc] = hapColor[colorIndex][previousLoc];
            }
            else
            {
                //cerr << "not common ";
                if (hapSeen.count(hapStr) == 0) //not seen
                {
                    //cerr << "not seen\n";
                    hapColor[colorIndex][currentLoc] = ++currentColor;
                    hapSeen[hapStr] = currentColor;
                }
                else // seen
                {
                    //cerr << "seen\n";
                    hapColor[colorIndex][currentLoc] = hapSeen[hapStr];
                }
            }
        }
        //Family did not split, so keep the color it had
        else
        {
            //cerr << "no split\n";
            hapColor[colorIndex][currentLoc] = hapColor[colorIndex][previousLoc];
        }
        colorIndex++;

        //string junk;
        //cin >> junk;
    }
    return;
}

bool familyDidSplit(const string &hapStr, const int hapCount, int **hapColor, const int nhaps, const int colorIndex, const int previousLoc, string &mostCommonHap)
{
    //cerr << "most common hap: " << mostCommonHap << endl;
    int previousColor = hapColor[colorIndex][previousLoc];
    int numPrevColor = 0;

    for (int i = 0; i < nhaps; i++)
    {
        if (hapColor[i][previousLoc] == previousColor) numPrevColor++;
    }

    //cerr << "numPrevColor: " << numPrevColor << "\nhapCount: " << hapCount << endl;

    if (numPrevColor == hapCount) return false;

    if ( hapCount > double(numPrevColor) / 2.0  )
    {
        mostCommonHap = hapStr;
        return true;
    }
    else if (mostCommonHap.compare("_NONE_") == 0 && hapCount == double(numPrevColor) / 2.0)
    {
        mostCommonHap = hapStr;
        return true;
    }
    else return true;
}

void calc_ihs(void *order)
{
    work_order_t *p = (work_order_t *)order;
    char **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nhaps = p->hapData->nhaps;
    int *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    string *locusName = p->mapData->locusName;
    int id = p->id;
    double *ihs = p->ihs;
    double *ihh1 = p->ihh1;
    double *ihh2 = p->ihh2;
    double *ihhDerivedLeft    = p->ihhDerivedLeft;
    double *ihhDerivedRight   = p->ihhDerivedRight;
    double *ihhAncestralLeft  = p->ihhAncestralLeft;
    double *ihhAncestralRight = p->ihhAncestralRight;
    double *freq = p->freq;
    ofstream *flog = p->flog;
    Bar *pbar = p->bar;

    int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    bool ALT = p->params->getBoolFlag(ARG_ALT);
    bool WAGH = p->params->getBoolFlag(ARG_WAGH);
    bool TRUNC = p->params->getBoolFlag(ARG_TRUNC);
    double MAF = p->params->getDoubleFlag(ARG_MAF);
    int numThreads = p->params->getIntFlag(ARG_THREAD);
    int MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND) <= 0 ) ? physicalPos[nloci - 1] - physicalPos[0] : p->params->getIntFlag(ARG_MAX_EXTEND);
    //bool SKIP = p->params->getBoolFlag(ARG_SKIP);
    bool WRITE_DETAILED_IHS = p->params->getBoolFlag(ARG_IHS_DETAILED);
    double (*calc)(map<string, int> &, int, bool) = p->calc;

    bool unphased = p->params->getBoolFlag(ARG_UNPHASED);

    //int step = (stop - start) / (pbar->totalTicks);
    int step = (nloci / numThreads) / (pbar->totalTicks);
    if (step == 0) step = 1;

    bool isDerived;
    string hapStr;

    for (int locus = id; locus < nloci; locus += numThreads)
    {
        if (locus % step == 0) advanceBar(*pbar, double(step));

        ihs[locus] = 0;
        //freq[locus] = MISSING;
        ihh1[locus] = MISSING;
        ihh2[locus] = MISSING;
        if (WRITE_DETAILED_IHS) {
            ihhDerivedLeft[locus]    = MISSING;
            ihhDerivedRight[locus]   = MISSING;
            ihhAncestralLeft[locus]  = MISSING;
            ihhAncestralRight[locus] = MISSING;
        }
        bool skipLocus = 0;
        //If the focal snp has MAF < MAF, then skip this locus
        if (freq[locus] < MAF || freq[locus] > 1 - MAF)
        {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Locus " << locusName[locus]
                    << " has MAF < " << MAF << ". Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
            pthread_mutex_unlock(&mutex_log);
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        //EHH to the left of the core snp
        double current_derived_ehh = 1;
        double current_ancestral_ehh = 1;
        double previous_derived_ehh = 1;
        double previous_ancestral_ehh = 1;
        int currentLocus = locus;
        int nextLocus = locus - 1;
        double derived_ihh = 0;
        double ancestral_ihh = 0;

        double derived_ihh_left    = 0;
        double ancestral_ihh_left  = 0;
        double derived_ihh_right   = 0;
        double ancestral_ihh_right = 0;

        //used for unphased analyses
        double current_notDerived_ehh = 1;
        double current_notAncestral_ehh = 1;
        double previous_notDerived_ehh = 1;
        double previous_notAncestral_ehh = 1;
        double notDerived_ihh = 0;
        double notAncestral_ihh = 0;

        //left right for unphased
        double notDerived_ihh_left    = 0;
        double notAncestral_ihh_left  = 0;
        double notDerived_ihh_right   = 0;
        double notAncestral_ihh_right = 0;

        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        string *haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            //derivedCount += ( data[hap][locus] == '1' ) ? 1 : 0;
            haplotypeList[hap] = data[hap][locus];
        }

        while (current_derived_ehh > EHH_CUTOFF || current_ancestral_ehh > EHH_CUTOFF)
        {
            if (nextLocus < 0)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[currentLocus] - physicalPos[nextLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[currentLocus] - physicalPos[nextLocus]
                        << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            //Check to see if the gap between the markers is huge, if so, scale it in an ad hoc way as in
            //Voight et al. (2006)
            double scale = double(SCALE_PARAMETER) / double(physicalPos[currentLocus] - physicalPos[nextLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus--;

            int numDerived = 0;
            int numAncestral = 0;
            map<string, int> ancestralHapCount;
            map<string, int> derivedHapCount;

            //used for unphased analyses
            int numHet = 0;
            map<string, int> notAncestralHapCount;
            map<string, int> notDerivedHapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                if(unphased){
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];
                    if (data[hap][locus] == '0'){
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                    }
                    else if (data[hap][locus] == '1'){
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                        numHet++;
                    }
                    else{//data[hap][locus] == '2'
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                    }
                }
                else{
                    isDerived = ( data[hap][locus] == '1' ) ? 1 : 0;
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];

                    if (isDerived)
                    {
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                    }
                    else
                    {
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                    }
                }
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            //If locus is monomorphic, shoot a warning and skip locus
            //This probably isn't necessary any more
            if ( !unphased && (numDerived == 0 || numAncestral == 0) ) 
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: locus " << locusName[locus]
                        << " (number " << locus + 1 << ") is monomorphic. Skipping calculation at this locus.\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double freqHetGT = double(numHet) / double(numDerived + numAncestral + numHet);
            //double freqAncestralGT = double(numAncestral) / double(numDerived + numAncestral + numHet);

            if ( unphased && freqHetGT > 1-MAF ) 
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: locus " << locusName[locus]
                        << " (number " << locus + 1 << ") has too many hets. Skipping calculation at this locus. "
                        << "het: " << numHet << " hom0: " << numAncestral << " hom1: " << numDerived << ".\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            if (current_derived_ehh > EHH_CUTOFF)
            {
                current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                derived_ihh += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_derived_ehh + previous_derived_ehh);
                previous_derived_ehh = current_derived_ehh;
            }

            if (current_ancestral_ehh > EHH_CUTOFF)
            {
                current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                ancestral_ihh += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_ancestral_ehh + previous_ancestral_ehh);
                previous_ancestral_ehh = current_ancestral_ehh;
            }

            if(unphased){
                current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notDerived_ihh += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_notDerived_ehh + previous_notDerived_ehh);
                previous_notDerived_ehh = current_notDerived_ehh;

                current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notAncestral_ihh += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_notAncestral_ehh + previous_notAncestral_ehh);
                previous_notAncestral_ehh = current_notAncestral_ehh;
            }

            //check if currentLocus is beyond MAX_EXTEND
            if (physicalPos[locus] - physicalPos[currentLocus] >= MAX_EXTEND) break;
        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        derived_ihh_left   = derived_ihh;
        ancestral_ihh_left = ancestral_ihh;

        if(unphased){
            notDerived_ihh_left    = notDerived_ihh;    
            notAncestral_ihh_left  = notAncestral_ihh;
        }
/*
        if(unphased){
                ihh1[locus] = log10(derived_ihh / notDerived_ihh);
                ihh2[locus] = log10(ancestral_ihh / notAncestral_ihh);
                ihs[locus] = (ihh1[locus] > ihh2[locus]) ? ihh1[locus] : 0-ihh2[locus];
            }
*/
        //calculate EHH to the right
        current_derived_ehh = 1;
        current_ancestral_ehh = 1;
        previous_derived_ehh = 1;
        previous_ancestral_ehh = 1;
        currentLocus = locus;
        nextLocus = locus + 1;
        skipLocus = 0;

        //used for unphased analyses
        current_notDerived_ehh = 1;
        current_notAncestral_ehh = 1;
        previous_notDerived_ehh = 1;
        previous_notAncestral_ehh = 1;
        
        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            //char digit[2];
            //snprintf(digit, "%d", data[hap][locus]);
            haplotypeList[hap] = data[hap][locus];
        }

        while (current_ancestral_ehh > EHH_CUTOFF || current_derived_ehh > EHH_CUTOFF)
        {
            if (nextLocus > nloci - 1)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[nextLocus] - physicalPos[currentLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[nextLocus] - physicalPos[currentLocus]
                        << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double scale = double(SCALE_PARAMETER) / double(physicalPos[nextLocus] - physicalPos[currentLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus++;

            int numDerived = 0;
            int numAncestral = 0;
            map<string, int> ancestralHapCount;
            map<string, int> derivedHapCount;

            //used for unphased analyses
            int numHet = 0;
            map<string, int> notAncestralHapCount;
            map<string, int> notDerivedHapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                if(unphased){
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];
                    if (data[hap][locus] == '0'){
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                    }
                    else if (data[hap][locus] == '1'){
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                        numHet++;
                    }
                    else{//data[hap][locus] == '2'
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                    }
                }
                else{
                    isDerived = ( data[hap][locus] == '1') ? 1 : 0;
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];

                    if (isDerived)
                    {
                        //count hapoltypes
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                    }
                    else //ancestral
                    {
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                    }
                }
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            //If there are no derived alleles at a locus, shoot a warning and skip locus
            if (numDerived == 0 || numAncestral == 0)
            {
                //(*flog) << "WARNING: locus " << locusName[locus]
                //   << " (number " << locus+1 << ") is monomorphic.  Skipping calculation at this locus.\n";
                skipLocus = 1;
                break;
            }

            if (current_derived_ehh > EHH_CUTOFF)
            {
                current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                derived_ihh += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_derived_ehh + previous_derived_ehh);
                previous_derived_ehh = current_derived_ehh;
            }

            if (current_ancestral_ehh > EHH_CUTOFF)
            {
                current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                ancestral_ihh += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_ancestral_ehh + previous_ancestral_ehh);
                previous_ancestral_ehh = current_ancestral_ehh;
            }

            if(unphased){
                current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notDerived_ihh += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_notDerived_ehh + previous_notDerived_ehh);
                previous_notDerived_ehh = current_notDerived_ehh;

                current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notAncestral_ihh += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_notAncestral_ehh + previous_notAncestral_ehh);
                previous_notAncestral_ehh = current_notAncestral_ehh;
            }

            //check if currentLocus is beyond 1Mb
            if (physicalPos[currentLocus] - physicalPos[locus] >= MAX_EXTEND) break;

        }

        derived_ihh_right   = derived_ihh   - derived_ihh_left;
        ancestral_ihh_right = ancestral_ihh - ancestral_ihh_left;

        if(unphased){
            notDerived_ihh_right = notDerived_ihh - notDerived_ihh_left;
            notAncestral_ihh_right = notAncestral_ihh - notAncestral_ihh_left;
        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        if (ihs[locus] != MISSING)
        {
            if(unphased){
                ihh1[locus] = log10(derived_ihh / notDerived_ihh);
                ihh2[locus] = log10(ancestral_ihh / notAncestral_ihh);
                ihs[locus] = (ihh1[locus] > ihh2[locus]) ? ihh1[locus] : 0-ihh2[locus];
                /*
                if(isnan(ihs[locus])){
                    cerr << "NAN: " << derived_ihh << " " << notDerived_ihh << " " << ancestral_ihh << " " << notAncestral_ihh << endl;
                }
                */
            }
            else{
                ihh1[locus] = derived_ihh;
                ihh2[locus] = ancestral_ihh;
                ihs[locus] = log10(derived_ihh / ancestral_ihh);
            }
            if (WRITE_DETAILED_IHS) {
                if(unphased){//for the time being this is going to report "same" as phased
                    ihhDerivedLeft[locus]    = derived_ihh_left;
                    ihhDerivedRight[locus]   = derived_ihh_right;
                    ihhAncestralLeft[locus]  = ancestral_ihh_left;
                    ihhAncestralRight[locus] = ancestral_ihh_right;
                }
                else{
                    ihhDerivedLeft[locus]    = derived_ihh_left;
                    ihhDerivedRight[locus]   = derived_ihh_right;
                    ihhAncestralLeft[locus]  = ancestral_ihh_left;
                    ihhAncestralRight[locus] = ancestral_ihh_right;
                }
            }
            //freq[locus] = double(derivedCount) / double(nhaps);
        }
    }
}

void calc_nsl(void *order)
{
    work_order_t *p = (work_order_t *)order;
    char **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nhaps = p->hapData->nhaps;
    int *physicalPos = p->mapData->physicalPos;
    //double *geneticPos = p->mapData->geneticPos;
    string *locusName = p->mapData->locusName;
    int id = p->id;
    double *ihs = p->ihs;
    double *ihh1 = p->ihh1;
    double *ihh2 = p->ihh2;
    double *freq = p->freq;
    ofstream *flog = p->flog;
    Bar *pbar = p->bar;

    int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    //double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    double EHH_CUTOFF = 0;
    bool ALT = p->params->getBoolFlag(ARG_ALT);
    bool WAGH = p->params->getBoolFlag(ARG_WAGH);
    bool TRUNC = p->params->getBoolFlag(ARG_TRUNC);
    double MAF = p->params->getDoubleFlag(ARG_MAF);
    int numThreads = p->params->getIntFlag(ARG_THREAD);
    int MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND_NSL) <= 0 ) ? nloci : p->params->getIntFlag(ARG_MAX_EXTEND_NSL);

    bool unphased = p->params->getBoolFlag(ARG_UNPHASED);

    double (*calc)(map<string, int> &, int, bool) = p->calc;

    int step = (nloci / numThreads) / (pbar->totalTicks);
    if (step == 0) step = 1;

    bool isDerived;
    string hapStr;

    for (int locus = id; locus < nloci; locus += numThreads)
    {
        if (locus % step == 0) advanceBar(*pbar, double(step));

        ihs[locus] = 0;
        //freq[locus] = MISSING;
        ihh1[locus] = MISSING;
        ihh2[locus] = MISSING;
        bool skipLocus = 0;
        //If the focal snp has MAF < MAF, then skip this locus
        if (freq[locus] < MAF || freq[locus] > 1 - MAF)
        {
            pthread_mutex_lock(&mutex_log);
            (*flog) << "WARNING: Locus " << locusName[locus]
                    << " has MAF < " << MAF << ". Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
            pthread_mutex_unlock(&mutex_log);
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        //EHH to the left of the core snp
        double current_derived_ehh = 1;
        double current_ancestral_ehh = 1;
        double previous_derived_ehh = 1;
        double previous_ancestral_ehh = 1;
        int currentLocus = locus;
        int nextLocus = locus - 1;
        double derived_ihh = 0;
        double ancestral_ihh = 0;
        
        //used for unphased analyses
        double current_notDerived_ehh = 1;
        double current_notAncestral_ehh = 1;
        double previous_notDerived_ehh = 1;
        double previous_notAncestral_ehh = 1;
        double notDerived_ihh = 0;
        double notAncestral_ihh = 0;

        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        string *haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            //derivedCount += ( data[hap][locus] == '1' ) ? 1 : 0;
            haplotypeList[hap] = data[hap][locus];
        }

        while (current_derived_ehh > EHH_CUTOFF || current_ancestral_ehh > EHH_CUTOFF)
        {
            if (nextLocus < 0)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[currentLocus] - physicalPos[nextLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[currentLocus] - physicalPos[nextLocus]
                        << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            //Check to see if the gap between the markers is huge, if so, scale it in an ad hoc way as in
            //Voight et al. (2006)
            double scale = double(SCALE_PARAMETER) / double(physicalPos[currentLocus] - physicalPos[nextLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus--;

            int numDerived = 0;
            int numAncestral = 0;

            map<string, int> ancestralHapCount;
            map<string, int> derivedHapCount;

            //used for unphased analyses
            int numHet = 0;
            map<string, int> notAncestralHapCount;
            map<string, int> notDerivedHapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                if(unphased){
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];
                    if (data[hap][locus] == '0'){
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                    }
                    else if (data[hap][locus] == '1'){
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                        numHet++;
                    }
                    else{//data[hap][locus] == '2'
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                    }
                }
                else{
                    isDerived = ( data[hap][locus] == '1' ) ? 1 : 0;
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];

                    if (isDerived)
                    {
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                    }
                    else
                    {
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                    }
                }
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            //If locus is monomorphic, shoot a warning and skip locus
            //This probably isn't necessary any more
            if ( !unphased && (numDerived == 0 || numAncestral == 0) ) 
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: locus " << locusName[locus]
                        << " (number " << locus + 1 << ") is monomorphic. Skipping calculation at this locus.\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double freqHetGT = double(numHet) / double(numDerived + numAncestral + numHet);
            //double freqAncestralGT = double(numAncestral) / double(numDerived + numAncestral + numHet);

            if ( unphased && freqHetGT > 1-MAF ) 
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: locus " << locusName[locus]
                        << " (number " << locus + 1 << ") has too many hets. Skipping calculation at this locus. "
                        << "het: " << numHet << " hom0: " << numAncestral << " hom1: " << numDerived << ".\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            if (current_derived_ehh > EHH_CUTOFF)
            {
                current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                derived_ihh += 0.5 * scale * (current_derived_ehh + previous_derived_ehh);
                previous_derived_ehh = current_derived_ehh;
            }

            if (current_ancestral_ehh > EHH_CUTOFF)
            {
                current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                ancestral_ihh += 0.5 * scale * (current_ancestral_ehh + previous_ancestral_ehh);
                previous_ancestral_ehh = current_ancestral_ehh;
            }

            if(unphased){
                current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notDerived_ihh += 0.5 * scale * (current_notDerived_ehh + previous_notDerived_ehh);
                previous_notDerived_ehh = current_notDerived_ehh;

                current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notAncestral_ihh += 0.5 * scale * (current_notAncestral_ehh + previous_notAncestral_ehh);
                previous_notAncestral_ehh = current_notAncestral_ehh;
            }

            //check if currentLocus is beyond MAX_EXTEND
            if (locus - currentLocus >= MAX_EXTEND) break;
        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        //calculate EHH to the right
        current_derived_ehh = 1;
        current_ancestral_ehh = 1;
        previous_derived_ehh = 1;
        previous_ancestral_ehh = 1;
        currentLocus = locus;
        nextLocus = locus + 1;
        skipLocus = 0;

        //used for unphased analyses
        current_notDerived_ehh = 1;
        current_notAncestral_ehh = 1;
        previous_notDerived_ehh = 1;
        previous_notAncestral_ehh = 1;
        
        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            //char digit[2];
            //snprintf(digit, "%d", data[hap][locus]);
            haplotypeList[hap] = data[hap][locus];
        }

        while (current_ancestral_ehh > EHH_CUTOFF || current_derived_ehh > EHH_CUTOFF)
        {
            if (nextLocus > nloci - 1)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[nextLocus] - physicalPos[currentLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[nextLocus] - physicalPos[currentLocus]
                        << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double scale = double(SCALE_PARAMETER) / double(physicalPos[nextLocus] - physicalPos[currentLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus++;

            int numDerived = 0;
            int numAncestral = 0;

            map<string, int> ancestralHapCount;
            map<string, int> derivedHapCount;

            //used for unphased analyses
            int numHet = 0;
            map<string, int> notAncestralHapCount;
            map<string, int> notDerivedHapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                if(unphased){
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];
                    if (data[hap][locus] == '0'){
                        //count ancestral haplotype
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                    }
                    else if (data[hap][locus] == '1'){
                        //count non-derived hapoltype
                        if (notDerivedHapCount.count(hapStr) == 0) notDerivedHapCount[hapStr] = 1;
                        else notDerivedHapCount[hapStr]++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                        numHet++;
                    }
                    else{//data[hap][locus] == '2'
                        //count derived hapoltype
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                        //count non-ancestral haplotype
                        if (notAncestralHapCount.count(hapStr) == 0) notAncestralHapCount[hapStr] = 1;
                        else notAncestralHapCount[hapStr]++;
                    }
                }
                else{
                    isDerived = ( data[hap][locus] == '1') ? 1 : 0;
                    haplotypeList[hap] += data[hap][currentLocus];
                    hapStr = haplotypeList[hap];

                    if (isDerived)
                    {
                        //count hapoltypes
                        if (derivedHapCount.count(hapStr) == 0) derivedHapCount[hapStr] = 1;
                        else derivedHapCount[hapStr]++;
                        numDerived++;
                    }
                    else //ancestral
                    {
                        if (ancestralHapCount.count(hapStr) == 0) ancestralHapCount[hapStr] = 1;
                        else ancestralHapCount[hapStr]++;
                        numAncestral++;
                    }
                }
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            //If there are no derived alleles at a locus, shoot a warning and skip locus
            if (numDerived == 0 || numAncestral == 0)
            {
                //(*flog) << "WARNING: locus " << locusName[locus]
                //   << " (number " << locus+1 << ") is monomorphic.  Skipping calculation at this locus.\n";
                skipLocus = 1;
                break;
            }

            if (current_derived_ehh > EHH_CUTOFF)
            {
                current_derived_ehh = (*calc)(derivedHapCount, numDerived, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                derived_ihh += 0.5 * scale * (current_derived_ehh + previous_derived_ehh);
                previous_derived_ehh = current_derived_ehh;
            }

            if (current_ancestral_ehh > EHH_CUTOFF)
            {
                current_ancestral_ehh = (*calc)(ancestralHapCount, numAncestral, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                ancestral_ihh += 0.5 * scale * (current_ancestral_ehh + previous_ancestral_ehh);
                previous_ancestral_ehh = current_ancestral_ehh;
            }

            if(unphased){
                current_notDerived_ehh = (*calc)(notDerivedHapCount, numAncestral + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notDerived_ihh += 0.5 * scale * (current_notDerived_ehh + previous_notDerived_ehh);
                previous_notDerived_ehh = current_notDerived_ehh;

                current_notAncestral_ehh = (*calc)(notAncestralHapCount, numDerived + numHet, ALT);

                //directly calculate ihs, iteratively
                //Trapezoid rule
                notAncestral_ihh += 0.5 * scale * (current_notAncestral_ehh + previous_notAncestral_ehh);
                previous_notAncestral_ehh = current_notAncestral_ehh;
            }

            //check if currentLocus is beyond MAX_EXTEND
            if (currentLocus - locus >= MAX_EXTEND) break;

        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            ihs[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        if (ihs[locus] != MISSING)
        {
            if(unphased){
                ihh1[locus] = log10(derived_ihh / notDerived_ihh);
                ihh2[locus] = log10(ancestral_ihh / notAncestral_ihh);
                ihs[locus] = (ihh1[locus] > ihh2[locus]) ? ihh1[locus] : 0-ihh2[locus];
            }
            else{
                ihh1[locus] = derived_ihh;
                ihh2[locus] = ancestral_ihh;
                ihs[locus] = log10(derived_ihh / ancestral_ihh);
                //freq[locus] = double(derivedCount) / double(nhaps);
            }
        }
    }
}


void calc_soft_ihs(void *order)
{
    work_order_t *p = (work_order_t *)order;
    char **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nhaps = p->hapData->nhaps;
    int *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    string *locusName = p->mapData->locusName;
    //int start = p->first_index;
    //int stop = p->last_index;
    int id = p->id;
    double *h1 = p->ihh1;
    double *h2dh1 = p->ihh2;
    double *h12 = p->ihs;
    //double *freq = p->freq;
    ofstream *flog = p->flog;
    Bar *pbar = p->bar;

    int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    //bool ALT = p->params->getBoolFlag(ARG_ALT);
    //double MAF = p->params->getDoubleFlag(ARG_MAF);
    //bool WAGH = p->params->getBoolFlag(ARG_WAGH);
    int numThreads = p->params->getIntFlag(ARG_THREAD);
    bool TRUNC = p->params->getBoolFlag(ARG_TRUNC);
    int MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND) <= 0 ) ? physicalPos[nloci - 1] - physicalPos[0] : p->params->getIntFlag(ARG_MAX_EXTEND);

    int step = (nloci / numThreads) / (pbar->totalTicks);
    if (step == 0) step = 1;

    for (int locus = id; locus < nloci; locus += numThreads)
    {
        if (locus % step == 0) advanceBar(*pbar, double(step));

        //freq[locus] = MISSING;
        h1[locus] = -1;
        h2dh1[locus] = -1;
        h12[locus] = -1;

        //EHH to the left of the core snp
        double current_ehh1 = 1;
        double current_ehh2d1 = 1;
        double current_ehh12 = 1;

        double previous_ehh1 = 1;
        double previous_ehh2d1 = 1;
        double previous_ehh12 = 1;

        int currentLocus = locus;
        int nextLocus = locus - 1;
        bool skipLocus = 0;
        double ihh1 = 0;
        double ihh2d1 = 0;
        double ihh12 = 0;
        //double derivedCount = 0;
        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        map<string, int> tempHapCount;
        string *haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            //derivedCount += data[hap][locus];
            haplotypeList[hap] = data[hap][locus];
            string hapStr = haplotypeList[hap];
            //count hapoltype freqs
            if (tempHapCount.count(hapStr) == 0) tempHapCount[hapStr] = 1;
            else tempHapCount[hapStr]++;
        }

        triplet_t res = calculateSoft(tempHapCount, nhaps);
        current_ehh1 = res.h1;
        current_ehh2d1 = res.h2dh1;
        current_ehh12 = res.h12;
        previous_ehh1 = res.h1;
        previous_ehh2d1 = res.h2dh1;
        previous_ehh12 = res.h12;

        while (current_ehh1 > EHH_CUTOFF)
        {
            if (nextLocus < 0)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[currentLocus] - physicalPos[nextLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[currentLocus] - physicalPos[nextLocus] << "bp > " << MAX_GAP 
                << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            //Check to see if the gap between the markers is huge, if so, scale it in an ad hoc way as in
            //Voight et al. (2006)
            double scale = double(SCALE_PARAMETER) / double(physicalPos[currentLocus] - physicalPos[nextLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus--;

            map<string, int> hapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                //build haplotype string
                haplotypeList[hap] += data[hap][currentLocus];
                string hapStr = haplotypeList[hap];

                //count hapoltype freqs
                if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
                else hapCount[hapStr]++;
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            res = calculateSoft(hapCount, nhaps);
            current_ehh1 = res.h1;
            current_ehh2d1 = res.h2dh1;
            current_ehh12 = res.h12;

            //directly calculate integral, iteratively
            //Trapezoid rule
            ihh1 += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_ehh1 + previous_ehh1);
            previous_ehh1 = current_ehh1;
            ihh2d1 += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_ehh2d1 + previous_ehh2d1);
            previous_ehh2d1 = current_ehh2d1;
            ihh12 += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_ehh12 + previous_ehh12);
            previous_ehh12 = current_ehh12;

            //check if currentLocus is beyond 1Mb
            if (physicalPos[locus] - physicalPos[currentLocus] >= MAX_EXTEND) break;
        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            h1[locus] = MISSING;
            h2dh1[locus] = MISSING;
            h12[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        //calculate EHH to the right
        current_ehh1 = 1;
        current_ehh2d1 = 1;
        current_ehh12 = 1;

        previous_ehh1 = 1;
        previous_ehh2d1 = 1;
        previous_ehh12 = 1;

        currentLocus = locus;
        nextLocus = locus + 1;
        skipLocus = 0;
        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        tempHapCount.clear();
        haplotypeList = new string[nhaps];
        for (int hap = 0; hap < nhaps; hap++)
        {
            haplotypeList[hap] = data[hap][locus];
            string hapStr = haplotypeList[hap];
            //count hapoltype freqs
            if (tempHapCount.count(hapStr) == 0) tempHapCount[hapStr] = 1;
            else tempHapCount[hapStr]++;
        }

        res = calculateSoft(tempHapCount, nhaps);
        current_ehh1 = res.h1;
        current_ehh2d1 = res.h2dh1;
        current_ehh12 = res.h12;
        previous_ehh1 = res.h1;
        previous_ehh2d1 = res.h2dh1;
        previous_ehh12 = res.h12;

        while (current_ehh1 > EHH_CUTOFF)
        {
            if (nextLocus > nloci - 1)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[nextLocus] - physicalPos[currentLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[nextLocus] - physicalPos[currentLocus] << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double scale = double(SCALE_PARAMETER) / double(physicalPos[nextLocus] - physicalPos[currentLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus++;

            map<string, int> hapCount;

            for (int hap = 0; hap < nhaps; hap++)
            {
                //build haplotype string
                haplotypeList[hap] += data[hap][currentLocus];
                string hapStr = haplotypeList[hap];

                //count hapoltype freqs
                if (hapCount.count(hapStr) == 0) hapCount[hapStr] = 1;
                else hapCount[hapStr]++;
            }

            //We've now counted all of the unique haplotypes extending out of the core SNP
            res = calculateSoft(hapCount, nhaps);
            current_ehh1 = res.h1;
            current_ehh2d1 = res.h2dh1;
            current_ehh12 = res.h12;

            //directly calculate integral, iteratively
            //Trapezoid rule
            ihh1 += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_ehh1 + previous_ehh1);
            previous_ehh1 = current_ehh1;
            ihh2d1 += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_ehh2d1 + previous_ehh2d1);
            previous_ehh2d1 = current_ehh2d1;
            ihh12 += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_ehh12 + previous_ehh12);
            previous_ehh12 = current_ehh12;

            //check if currentLocus is beyond 1Mb
            if (physicalPos[currentLocus] - physicalPos[locus] >= MAX_EXTEND) break;
        }

        delete [] haplotypeList;

        if (skipLocus == 1)
        {
            h1[locus] = MISSING;
            h2dh1[locus] = MISSING;
            h12[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        if (h12[locus] != MISSING)
        {
            h1[locus] = ihh1;
            h2dh1[locus] = ihh2d1;
            h12[locus] = ihh12;
            //freq[locus] = double(derivedCount) / double(nhaps);
        }
    }
}


void calc_xpihh(void *order)
{
    work_order_t *p = (work_order_t *)order;

    char **data1 = p->hapData1->data;
    int nhaps1 = p->hapData1->nhaps;
    double *ihh1 = p->ihh1;
    //double *freq1 = p->freq1;

    char **data2 = p->hapData2->data;
    int nhaps2 = p->hapData2->nhaps;
    double *ihh2 = p->ihh2;
    //double *freq2 = p->freq2;

    int nloci = p->mapData->nloci;
    int *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    string *locusName = p->mapData->locusName;

    int id = p->id;

    ofstream *flog = p->flog;
    Bar *pbar = p->bar;

    bool TRUNC = p->params->getBoolFlag(ARG_TRUNC);
    int SCALE_PARAMETER = p->params->getIntFlag(ARG_GAP_SCALE);
    int MAX_GAP = p->params->getIntFlag(ARG_MAX_GAP);
    double EHH_CUTOFF = p->params->getDoubleFlag(ARG_CUTOFF);
    bool ALT = p->params->getBoolFlag(ARG_ALT);
    bool WAGH = p->params->getBoolFlag(ARG_WAGH);
    int numThreads = p->params->getIntFlag(ARG_THREAD);
    bool CALC_XPNSL = p->params->getBoolFlag(ARG_XPNSL);

    bool unphased = p->params->getBoolFlag(ARG_UNPHASED);

    int MAX_EXTEND;
    if (!CALC_XPNSL){
        MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND) <= 0 ) ? physicalPos[nloci - 1] - physicalPos[0] : p->params->getIntFlag(ARG_MAX_EXTEND);
    }
    else{
        MAX_EXTEND = ( p->params->getIntFlag(ARG_MAX_EXTEND_NSL) <= 0 ) ? geneticPos[nloci - 1] - geneticPos[0] : p->params->getIntFlag(ARG_MAX_EXTEND_NSL);
    }
    int step = (nloci / numThreads) / (pbar->totalTicks);
    if (step == 0) step = 1;

    string hapStr;

    for (int locus = id; locus < nloci; locus += numThreads)
    {
        if (locus % step == 0) advanceBar(*pbar, double(step));

        ihh1[locus] = 0;
        ihh2[locus] = 0;

        //freq1[locus] = MISSING;
        //freq2[locus] = MISSING;

        //EHH to the left of the core snp
        double current_pooled_ehh = 1;
        double previous_pooled_ehh = 1;
        double derivedCountPooled = 0;

        double current_pop1_ehh = 1;
        double previous_pop1_ehh = 1;
        double ihhPop1 = 0;
        double derivedCount1 = 0;

        double current_pop2_ehh = 1;
        double previous_pop2_ehh = 1;
        double ihhPop2 = 0;
        double derivedCount2 = 0;

        //For unphased analyses
        double ancestralCount1 = 0;
        double ancestralCount2 = 0;
        double ancestralCountPooled = 0;
        double hetCount1 = 0;
        double hetCount2 = 0;
        double hetCountPooled = 0;

        int currentLocus = locus;
        int nextLocus = locus - 1;
        bool skipLocus = 0;

        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        string *haplotypeList1, *haplotypeList2, *haplotypeListPooled;
        haplotypeList1 = new string[nhaps1];
        haplotypeList2 = new string[nhaps2];
        haplotypeListPooled = new string[nhaps1 + nhaps2];
        for (int hap = 0; hap < nhaps1 + nhaps2; hap++)
        {
            //char digit[2];
            //Pop1
            if (hap < nhaps1)
            {
                //snprintf(digit, "%d", data1[hap][locus]);
                haplotypeList1[hap] = data1[hap][locus];
                haplotypeListPooled[hap] = data1[hap][locus];
                if(unphased){
                    derivedCount1 += ( data1[hap][locus] == '2' ) ? 1 : 0;
                    ancestralCount1 += ( data1[hap][locus] == '0' ) ? 1 : 0;
                    hetCount1 += ( data1[hap][locus] == '1' ) ? 1 : 0;
                }
                else{
                    derivedCount1 += ( data1[hap][locus] == '1' ) ? 1 : 0;
                    ancestralCount1 += ( data1[hap][locus] == '0' ) ? 1 : 0;
                }
            }
            //Pop2
            else
            {
                //snprintf(digit, "%d", data2[hap - nhaps1][locus]);
                haplotypeList2[hap - nhaps1] = data2[hap - nhaps1][locus];
                haplotypeListPooled[hap] = data2[hap - nhaps1][locus];
                if(unphased){
                    derivedCount2 += ( data2[hap - nhaps1][locus] == '2' ) ? 1 : 0;
                    ancestralCount2 += ( data2[hap - nhaps1][locus] == '0' ) ? 1 : 0;
                    hetCount2 += ( data2[hap - nhaps1][locus] == '1' ) ? 1 : 0;
                }
                else{
                    derivedCount2 += ( data2[hap - nhaps1][locus] == '1' ) ? 1 : 0;
                    ancestralCount2 += ( data2[hap - nhaps1][locus] == '0' ) ? 1 : 0;
                }
            }
        }

        derivedCountPooled = derivedCount1 + derivedCount2;
        ancestralCountPooled = ancestralCount1 + ancestralCount2;
        hetCountPooled = hetCount1 + hetCount2;

        //when calculating xp-ehh, ehh does not necessarily start at 1
        if (ALT)
        {
            double fD = double(derivedCount1) / double(nhaps1);
            double fA = double(ancestralCount1) / double(nhaps1);
            double fH = double(hetCount1) / double(nhaps1);

            current_pop1_ehh = fD * fD + fA * fA + fH * fH;
            previous_pop1_ehh = current_pop1_ehh;

            fD = double(derivedCount2) / double(nhaps2);
            fA = double(ancestralCount2) / double(nhaps2);
            fH = double(hetCount2) / double(nhaps2);
            current_pop2_ehh = fD * fD + fA * fA + fH * fH;
            previous_pop2_ehh = current_pop2_ehh;

            fD = double(derivedCountPooled) / double(nhaps1 + nhaps2);
            fA = double(ancestralCountPooled) / double(nhaps1 + nhaps2);
            fH = double(hetCountPooled) / double(nhaps1 + nhaps2);
            current_pooled_ehh = fD * fD + fA * fA + fH * fH;
            previous_pooled_ehh = current_pooled_ehh;
        }
        else
        {
            if (WAGH)
            {

                current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1,2) / (nCk(derivedCount1,2)+nCk(nhaps1-derivedCount1,2)) : 0;
                current_pop1_ehh += (nhaps1 - derivedCount1 > 1) ? nCk(nhaps1-derivedCount1,2) / (nCk(derivedCount1,2)+nCk(nhaps1-derivedCount1,2)) : 0;
                previous_pop1_ehh = current_pop1_ehh;

                current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / (nCk(derivedCount2,2)+nCk(nhaps2-derivedCount2,2)) : 0;
                current_pop2_ehh += (nhaps2 - derivedCount2 > 1) ? nCk(nhaps2 - derivedCount2, 2) / (nCk(derivedCount2,2)+nCk(nhaps2-derivedCount2,2)) : 0;
                previous_pop2_ehh = current_pop2_ehh;

            }
            else
            {
                current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1, 2) / nCk(nhaps1, 2) : 0;
                current_pop1_ehh += (ancestralCount1 > 1) ? nCk(ancestralCount1, 2) / nCk(nhaps1, 2) : 0;
                current_pop1_ehh += (hetCount1 > 1) ? nCk(hetCount1, 2) / nCk(nhaps1, 2) : 0;

                previous_pop1_ehh = current_pop1_ehh;

                current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / nCk(nhaps2, 2) : 0;
                current_pop2_ehh += (ancestralCount2 > 1) ? nCk(ancestralCount2, 2) / nCk(nhaps2, 2) : 0;
                current_pop2_ehh += (hetCount2 > 1) ? nCk(hetCount2, 2) / nCk(nhaps2, 2) : 0;
                previous_pop2_ehh = current_pop2_ehh;        

            }

            current_pooled_ehh = (derivedCountPooled > 1) ? nCk(derivedCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
            current_pooled_ehh += (ancestralCountPooled > 1) ? nCk(ancestralCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
            current_pooled_ehh += (hetCountPooled > 1) ? nCk(hetCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
            previous_pooled_ehh = current_pooled_ehh;    
        }

        while (current_pooled_ehh > EHH_CUTOFF)
        {
            if (nextLocus < 0)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[currentLocus] - physicalPos[nextLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[currentLocus] - physicalPos[nextLocus] << "bp > " << MAX_GAP 
                << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            //Check to see if the gap between the markers is huge, if so, scale it in an ad hoc way as in
            //Voight, et al. paper
            double scale = double(SCALE_PARAMETER) / double(physicalPos[currentLocus] - physicalPos[nextLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus--;

            map<string, int> hapCount1;
            map<string, int> hapCount2;
            map<string, int> hapCountPooled;

            //build haplotype strings
            for (int hap = 0; hap < nhaps1 + nhaps2; hap++)
            {
                //char digit[2];
                if (hap < nhaps1) //Pop1
                {
                    //snprintf(digit, "%d", data1[hap][currentLocus]);
                    haplotypeList1[hap] += data1[hap][currentLocus];
                    hapStr = haplotypeList1[hap];
                    if (hapCount1.count(hapStr) == 0) hapCount1[hapStr] = 1;
                    else hapCount1[hapStr]++;

                    //Pooled
                    haplotypeListPooled[hap] += data1[hap][currentLocus];
                    hapStr = haplotypeListPooled[hap];
                    if (hapCountPooled.count(hapStr) == 0) hapCountPooled[hapStr] = 1;
                    else hapCountPooled[hapStr]++;
                }
                else //Pop2
                {
                    //snprintf(digit, "%d", data2[hap - nhaps1][currentLocus]);
                    haplotypeList2[hap - nhaps1] += data2[hap - nhaps1][currentLocus];
                    hapStr = haplotypeList2[hap - nhaps1];
                    if (hapCount2.count(hapStr) == 0) hapCount2[hapStr] = 1;
                    else hapCount2[hapStr]++;

                    //Pooled
                    haplotypeListPooled[hap] += data2[hap - nhaps1][currentLocus];
                    hapStr = haplotypeListPooled[hap];
                    if (hapCountPooled.count(hapStr) == 0) hapCountPooled[hapStr] = 1;
                    else hapCountPooled[hapStr]++;
                }
            }

            if (WAGH)
            {
                current_pop1_ehh = calculateHomozygosity_Wagh(hapCount1,nhaps1,derivedCount1);
                current_pop2_ehh = calculateHomozygosity_Wagh(hapCount2,nhaps2,derivedCount2);

            }
            else
            {
                current_pop1_ehh = calculateHomozygosity(hapCount1, nhaps1, ALT);
                current_pop2_ehh = calculateHomozygosity(hapCount2, nhaps2, ALT);

            }

                       
            current_pooled_ehh = calculateHomozygosity(hapCountPooled, nhaps1 + nhaps2, ALT);

            //directly calculate ihh, iteratively
            //Trapezoid rule
            ihhPop1 += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_pop1_ehh + previous_pop1_ehh);
            previous_pop1_ehh = current_pop1_ehh;

            ihhPop2 += 0.5 * scale * (geneticPos[currentLocus + 1] - geneticPos[currentLocus]) * (current_pop2_ehh + previous_pop2_ehh);
            previous_pop2_ehh = current_pop2_ehh;

            previous_pooled_ehh = current_pooled_ehh;

            //check if currentLocus is beyond MAX_EXTEND
            if (!CALC_XPNSL && physicalPos[locus] - physicalPos[currentLocus] >= MAX_EXTEND) break;
            if (CALC_XPNSL && geneticPos[locus] - geneticPos[currentLocus] >= MAX_EXTEND) break;
        }

        delete [] haplotypeList1;
        delete [] haplotypeList2;
        delete [] haplotypeListPooled;

        if (skipLocus == 1)
        {
            ihh1[locus] = MISSING;
            ihh2[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        //calculate EHH to the right

        current_pooled_ehh = 1;
        previous_pooled_ehh = 1;
        derivedCountPooled = 0;

        current_pop1_ehh = 1;
        previous_pop1_ehh = 1;
        derivedCount1 = 0;

        current_pop2_ehh = 1;
        previous_pop2_ehh = 1;
        derivedCount2 = 0;

        //For unphased analyses
        ancestralCount1 = 0;
        ancestralCount2 = 0;
        ancestralCountPooled = 0;
        hetCount1 = 0;
        hetCount2 = 0;
        hetCountPooled = 0;

        currentLocus = locus;
        nextLocus = locus + 1;
        skipLocus = 0;

        //A list of all the haplotypes
        //Starts with just the focal snp and grows outward
        haplotypeList1 = new string[nhaps1];
        haplotypeList2 = new string[nhaps2];
        haplotypeListPooled = new string[nhaps1 + nhaps2];
for (int hap = 0; hap < nhaps1 + nhaps2; hap++)
        {
            //char digit[2];
            //Pop1
            if (hap < nhaps1)
            {
                //snprintf(digit, "%d", data1[hap][locus]);
                haplotypeList1[hap] = data1[hap][locus];
                haplotypeListPooled[hap] = data1[hap][locus];
                if(unphased){
                    derivedCount1 += ( data1[hap][locus] == '2' ) ? 1 : 0;
                    ancestralCount1 += ( data1[hap][locus] == '0' ) ? 1 : 0;
                    hetCount1 += ( data1[hap][locus] == '1' ) ? 1 : 0;
                }
                else{
                    derivedCount1 += ( data1[hap][locus] == '1' ) ? 1 : 0;
                    ancestralCount1 += ( data1[hap][locus] == '0' ) ? 1 : 0;
                }
            }
            //Pop2
            else
            {
                //snprintf(digit, "%d", data2[hap - nhaps1][locus]);
                haplotypeList2[hap - nhaps1] = data2[hap - nhaps1][locus];
                haplotypeListPooled[hap] = data2[hap - nhaps1][locus];
                if(unphased){
                    derivedCount2 += ( data2[hap - nhaps1][locus] == '2' ) ? 1 : 0;
                    ancestralCount2 += ( data2[hap - nhaps1][locus] == '0' ) ? 1 : 0;
                    hetCount2 += ( data2[hap - nhaps1][locus] == '1' ) ? 1 : 0;
                }
                else{
                    derivedCount2 += ( data2[hap - nhaps1][locus] == '1' ) ? 1 : 0;
                    ancestralCount2 += ( data2[hap - nhaps1][locus] == '0' ) ? 1 : 0;
                }
            }
        }

        derivedCountPooled = derivedCount1 + derivedCount2;
        ancestralCountPooled = ancestralCount1 + ancestralCount2;
        hetCountPooled = hetCount1 + hetCount2;

        //when calculating xp-ehh, ehh does not necessarily start at 1
        if (ALT)
        {
            double fD = double(derivedCount1) / double(nhaps1);
            double fA = double(ancestralCount1) / double(nhaps1);
            double fH = double(hetCount1) / double(nhaps1);

            current_pop1_ehh = fD * fD + fA * fA + fH * fH;
            previous_pop1_ehh = current_pop1_ehh;

            fD = double(derivedCount2) / double(nhaps2);
            fA = double(ancestralCount2) / double(nhaps2);
            fH = double(hetCount2) / double(nhaps2);
            current_pop2_ehh = fD * fD + fA * fA + fH * fH;
            previous_pop2_ehh = current_pop2_ehh;

            fD = double(derivedCountPooled) / double(nhaps1 + nhaps2);
            fA = double(ancestralCountPooled) / double(nhaps1 + nhaps2);
            fH = double(hetCountPooled) / double(nhaps1 + nhaps2);
            current_pooled_ehh = fD * fD + fA * fA + fH * fH;
            previous_pooled_ehh = current_pooled_ehh;
        }
        else
        {
            if (WAGH)
            {

                current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1,2) / (nCk(derivedCount1,2)+nCk(nhaps1-derivedCount1,2)) : 0;
                current_pop1_ehh += (nhaps1 - derivedCount1 > 1) ? nCk(nhaps1-derivedCount1,2) / (nCk(derivedCount1,2)+nCk(nhaps1-derivedCount1,2)) : 0;
                previous_pop1_ehh = current_pop1_ehh;

                current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / (nCk(derivedCount2,2)+nCk(nhaps2-derivedCount2,2)) : 0;
                current_pop2_ehh += (nhaps2 - derivedCount2 > 1) ? nCk(nhaps2 - derivedCount2, 2) / (nCk(derivedCount2,2)+nCk(nhaps2-derivedCount2,2)) : 0;
                previous_pop2_ehh = current_pop2_ehh;

            }
            else
            {
                current_pop1_ehh = (derivedCount1 > 1) ? nCk(derivedCount1, 2) / nCk(nhaps1, 2) : 0;
                current_pop1_ehh += (ancestralCount1 > 1) ? nCk(ancestralCount1, 2) / nCk(nhaps1, 2) : 0;
                current_pop1_ehh += (hetCount1 > 1) ? nCk(hetCount1, 2) / nCk(nhaps1, 2) : 0;

                previous_pop1_ehh = current_pop1_ehh;

                current_pop2_ehh = (derivedCount2 > 1) ? nCk(derivedCount2, 2) / nCk(nhaps2, 2) : 0;
                current_pop2_ehh += (ancestralCount2 > 1) ? nCk(ancestralCount2, 2) / nCk(nhaps2, 2) : 0;
                current_pop2_ehh += (hetCount2 > 1) ? nCk(hetCount2, 2) / nCk(nhaps2, 2) : 0;
                previous_pop2_ehh = current_pop2_ehh;        

            }

            current_pooled_ehh = (derivedCountPooled > 1) ? nCk(derivedCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
            current_pooled_ehh += (ancestralCountPooled > 1) ? nCk(ancestralCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
            current_pooled_ehh += (hetCountPooled > 1) ? nCk(hetCountPooled, 2) / nCk(nhaps1 + nhaps2, 2) : 0;
            previous_pooled_ehh = current_pooled_ehh;    
        }


        while (current_pooled_ehh > EHH_CUTOFF)
        {
            if (nextLocus > nloci - 1)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached chromosome edge before EHH decayed below " << EHH_CUTOFF
                        << ". ";
                if (!TRUNC){
                    skipLocus = 1;
                    (*flog) << "Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus];
                }
                (*flog) << "\n";
                pthread_mutex_unlock(&mutex_log);
                break;
            }
            else if (physicalPos[nextLocus] - physicalPos[currentLocus] > MAX_GAP)
            {
                pthread_mutex_lock(&mutex_log);
                (*flog) << "WARNING: Reached a gap of " << physicalPos[nextLocus] - physicalPos[currentLocus] << "bp > " << MAX_GAP << "bp. Skipping calculation at position " << physicalPos[locus] << " id: " << locusName[locus] << "\n";
                pthread_mutex_unlock(&mutex_log);
                skipLocus = 1;
                break;
            }

            double scale = double(SCALE_PARAMETER) / double(physicalPos[nextLocus] - physicalPos[currentLocus]);
            if (scale > 1) scale = 1;

            currentLocus = nextLocus;
            nextLocus++;

            map<string, int> hapCount1;
            map<string, int> hapCount2;
            map<string, int> hapCountPooled;

            //build haplotype strings
            for (int hap = 0; hap < nhaps1 + nhaps2; hap++)
            {
                //char digit[2];
                //string hapStr;

                //pop1
                if (hap < nhaps1)
                {
                    haplotypeList1[hap] += data1[hap][currentLocus];
                    hapStr = haplotypeList1[hap];
                    if (hapCount1.count(hapStr) == 0) hapCount1[hapStr] = 1;
                    else hapCount1[hapStr]++;

                    //Pooled
                    haplotypeListPooled[hap] += data1[hap][currentLocus];
                    hapStr = haplotypeListPooled[hap];
                    if (hapCountPooled.count(hapStr) == 0) hapCountPooled[hapStr] = 1;
                    else hapCountPooled[hapStr]++;
                }
                //Pop2
                else
                {
                    haplotypeList2[hap - nhaps1] += data2[hap - nhaps1][currentLocus];
                    hapStr = haplotypeList2[hap - nhaps1];
                    if (hapCount2.count(hapStr) == 0) hapCount2[hapStr] = 1;
                    else hapCount2[hapStr]++;

                    //Pooled
                    haplotypeListPooled[hap] += data2[hap - nhaps1][currentLocus];
                    hapStr = haplotypeListPooled[hap];
                    if (hapCountPooled.count(hapStr) == 0) hapCountPooled[hapStr] = 1;
                    else hapCountPooled[hapStr]++;
                }
            }

            if (WAGH)
            {
                current_pop1_ehh = calculateHomozygosity_Wagh(hapCount1,nhaps1,derivedCount1);
                current_pop2_ehh = calculateHomozygosity_Wagh(hapCount2,nhaps2,derivedCount2);

            }
            else
            {

                current_pop1_ehh = calculateHomozygosity(hapCount1, nhaps1, ALT);
                current_pop2_ehh = calculateHomozygosity(hapCount2, nhaps2, ALT);
            }


            current_pooled_ehh = calculateHomozygosity(hapCountPooled, nhaps1 + nhaps2, ALT);

            //directly calculate ihh1, iteratively
            //Trapezoid rule
            ihhPop1 += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_pop1_ehh + previous_pop1_ehh);
            previous_pop1_ehh = current_pop1_ehh;

            ihhPop2 += 0.5 * scale * (geneticPos[currentLocus] - geneticPos[currentLocus - 1]) * (current_pop2_ehh + previous_pop2_ehh);
            previous_pop2_ehh = current_pop2_ehh;

            previous_pooled_ehh = current_pooled_ehh;

            //check if currentLocus is beyond 1Mb
            if (!CALC_XPNSL && physicalPos[currentLocus] - physicalPos[locus] >= MAX_EXTEND) break;
            if (CALC_XPNSL && geneticPos[currentLocus] - geneticPos[locus] >= MAX_EXTEND) break;
        }

        delete [] haplotypeList1;
        delete [] haplotypeList2;
        delete [] haplotypeListPooled;

        if (skipLocus == 1)
        {
            ihh1[locus] = MISSING;
            ihh2[locus] = MISSING;
            skipLocus = 0;
            continue;
        }

        if (ihh1[locus] != MISSING)
        {
            ihh1[locus] = ihhPop1;
            //freq1[locus] = double(derivedCount1) / double(nhaps1);
        }

        if (ihh2[locus] != MISSING)
        {
            ihh2[locus] = ihhPop2;
            //freq2[locus] = double(derivedCount2) / double(nhaps2);
        }
    }
}

double calculateHomozygosity_Wagh(map<string, int> &count, int total, int derivedCount)
{
    double freq = 0;
    double homozygosity = 0;
    map<string, int>::iterator it;
    for (it = count.begin(); it != count.end(); it++)
    {
        homozygosity += (it->second > 1) ? nCk(it->second,2)/(nCk(derivedCount, 2) + nCk(total-derivedCount, 2)) : 0;
    }
        
    return homozygosity;
}


double calculateHomozygosity(map<string, int> &count, int total, bool ALT) // Called by XP-EHH
{
    double freq = 0;
    double homozygosity = 0;
    map<string, int>::iterator it;
    for (it = count.begin(); it != count.end(); it++)
    {
        if (ALT)
        {
            freq = double(it->second) / double(total);
            homozygosity += freq * freq;
        }

        else
        {
            homozygosity += (it->second > 1) ? nCk(it->second, 2) / nCk(total, 2) : 0;
        }
    }

    return homozygosity;
}

//Have to modify the function to handle the arbitrary declaration
//of EHH1K_VALUES
//If we pass, the MAX of that vector, plus the vector iteself
//We should be able to track up to the max in an array (up to what we do now for k = 2)
//return an array that is pre-allocated to the side of EKK1k_VALUES,
//and the index of that array corresponds to the indicies of EHH1K_VALUES
triplet_t calculateSoft(map<string, int> &count, int total)
{
    triplet_t res;
    res.h1 = 0;
    res.h2dh1 = 0;
    res.h12 = 0;

    double first = 0;
    double second = 0;
    //double freq = 0;
    double homozygosity = 0;
    map<string, int>::iterator it;
    for (it = count.begin(); it != count.end(); it++)
    {
        homozygosity += (it->second > 1) ? nCk(it->second, 2) / nCk(total, 2) : 0;
        //We can either do what is effectively a bubble sort here for the first K
        //sorted values
        //Or we can do a complete sort outside of the for loop
        //Will have to think about the scaling issues here
        if (it->second > first)
        {
            second = first;
            first = it->second;
        }
        else if (it->second > second)
        {
            second = it->second;
        }
    }

    //here we have a forloop over the EHH1k_VALUES elements
    double firstFreq = (first > 1) ? nCk(first, 2) / nCk(total, 2) : 0;
    double secondFreq = (second > 1) ? nCk(second, 2) / nCk(total, 2) : 0;
    double comboFreq = ((first + second) > 1) ? nCk((first + second), 2) / nCk(total, 2) : 0;

    res.h1 = homozygosity;
    res.h2dh1 = (homozygosity - firstFreq) / homozygosity;
    res.h12 = homozygosity - firstFreq - secondFreq + comboFreq;

    return res;
}

