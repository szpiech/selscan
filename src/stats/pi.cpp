

#include "pi.h"

void PI::main()
{
     double denominator = (hm->hapData->nhaps) * (hm->hapData->nhaps - 1) * 0.5;
    double pi = 0;

    int max_physpos = hm->mapData->mapEntries[hm->mapData->nloci - 1].physicalPos;
    int num_entries = ceil((max_physpos+1) / winsize);
    cout<<"DEBUG::: "<<"num entries"<<num_entries<<endl;

//    vector<double> index_by_start(num_entries+1);
    int site_counter;

    int locus = 0;
    for(int i =0; i<= num_entries; i++){
        (*fout)<<i*winsize + 1<<" " <<(i+1)*winsize  <<" ";

        int start = i*winsize + 1;
        int end = (i+1)*winsize;

        pi = 0;
        site_counter = 0;
        while(true){
            int pos = hm->mapData->mapEntries[locus].physicalPos;
            if(pos > end){
                break;
            }else{
                int n1 = hm->hapData->get_n_c1(locus);
                int n0 =  hm->hapData->nhaps - n1;
                pi += n1 * n0 ;
                site_counter++;

                locus++;
                if(locus == hm->mapData->nloci){
                    break;
                }
            }
        }

        (*fout) << pi / denominator   << endl;
    }

    
    exit(0);
    // for(int i = 1; i< mapData->nloci; i++){
    //     index_by_start.push_back(mapData->mapEntries[i].physicalPos);
    // }



    // int start = 1;
    // //int end = mapData->physicalPos[mapData->nloci - 1];
    // int startLocus = 0;
    // int endLocus = 0;
    // //Identify the start and end indicies in the window
    // int winStart = start;
    // int winEnd = winStart + winsize - 1;
    // int pos;
    // int length = 0;
    // double pi = 0;
    // double denominator = (hm->hapData->nhaps) * (hm->hapData->nhaps - 1) * 0.5;
    // int locus;


    // for (locus = 0; locus < hm->mapData->nloci; locus++)
    // {
    //     pos = hm->mapData->mapEntries[locus].physicalPos;
    //     while (pos > winEnd)
    //     {
    //         endLocus = locus - 1;
    //         length = endLocus - startLocus + 1;
    //         (*fout) << winStart << " " << winEnd << " ";
    //         //do stuff

    //         if (length == 0)
    //         {
    //             pi = 0;
    //         }
    //         else
    //         {
    //             int n1 = hm->hapData->get_n_c1(startLocus);
    //             int n0 =  hm->hapData->nhaps - n1;
    //             pi += n1 * n0 ;
    //             // for (int i = 0; i < hm->hapData->nhaps; i++)
    //             // {
    //             //     for (int j = i + 1; j < hm->hapData->nhaps; j++)
    //             //     {
    //             //         pi += hamming_dist_ptr(hapData->data[i] + startLocus, hapData->data[j] + startLocus, length);
    //             //     }
    //             // }
    //         }

            

    //         winStart += winsize;
    //         winEnd += winsize;

    //         startLocus = locus;
    //         //pi = 0;
    //         (*fout) << pi / denominator  / length  << endl;
    //     }
        
    // }


    // //final window
    // endLocus = locus - 1;
    // length = endLocus - startLocus + 1;
    // (*fout) << winStart << " " << winEnd << " ";
    // //do stuff

    // if (length == 0)
    // {
    //     pi = 0;
    // }
    // else
    // {
    //     int n1 = hm->hapData->get_n_c1(startLocus);
    //     int n0 =  hm->hapData->nhaps - n1;
    //     pi += n1 * n0;
    //     // for (int i = 0; i < hm->hapData->nhaps; i++)
    //     // {
    //     //     for (int j = i + 1; j < hm->hapData->nhaps; j++)
    //     //     {
    //     //         pi += hamming_dist_ptr(hapData->data[i] + startLocus, hapData->data[j] + startLocus, length);
    //     //     }
    //     // }
    // }

    // (*fout) << pi / denominator * length << endl;

    return;
}