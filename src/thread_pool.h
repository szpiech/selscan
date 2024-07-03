
#ifndef __THREAD_POOL_H__
#define __THREAD_POOL_H__

#include <iostream>
#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <unordered_map>
//#include "selscan-maintools.h"


/*
#include <atomic> 

class ThreadPool {
public:
    ThreadPool(size_t numThreads);
    ~ThreadPool();
    
    void enqueue(std::shared_ptr<MainTools> task);

private:
    std::vector<std::thread> workers;
    std::queue<std::shared_ptr<MainTools>> tasks;

    std::mutex queueMutex;
    std::condition_variable condition;
    std::atomic<bool> stop;

    void workerThread();
};
*/


class ThreadPool {
public:
    ThreadPool(size_t threads);
    ~ThreadPool();
    
    //template <class F>
    //void enqueue(F&& f);
    //void enqueue(const std::function<void(std::unordered_map< unsigned int, std::vector<unsigned int> >&)>& task, IHS& ehh_obj); ;
    void enqueue(const std::function<void(std::unordered_map< unsigned int, std::vector<unsigned int> >&  )>& task);
    //void enqueue(const std::function<void(std::map<int, std::string>&)>& task);


private:
    std::vector<std::thread> workers;
    std::vector<std::unordered_map< unsigned int, std::vector<unsigned int> > > worker_maps;

    std::queue<std::function<void(std::unordered_map< unsigned int, std::vector<unsigned int> >&)>> tasks;
    
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

// ThreadPool::ThreadPool(size_t threads) : stop(false) {
//     for(size_t i = 0; i < threads; ++i)
//         workers.emplace_back(
//             [this] {
//                 for(;;) {
//                     std::function<void()> task;
//                     {
//                         std::unique_lock<std::mutex> lock(this->queue_mutex);
//                         this->condition.wait(lock, [this]{ return this->stop || !this->tasks.empty(); });
//                         if(this->stop && this->tasks.empty())
//                             return;
//                         task = std::move(this->tasks.front());
//                         this->tasks.pop();
//                     }
//                     task();
//                 }
//             }
//         );
// }


#endif