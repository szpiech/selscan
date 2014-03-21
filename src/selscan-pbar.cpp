#include "selscan-pbar.h"

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
