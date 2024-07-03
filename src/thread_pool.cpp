#include "thread_pool.h"
//#include "stats/ihs.h"


#include <iostream>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <functional>
#include <memory>


ThreadPool::ThreadPool(size_t threads) : stop(false) {
    worker_maps.resize(threads);
    for (size_t i = 0; i < threads; ++i) {
        worker_maps[i] = {};
        workers.emplace_back(
            [this, i] {
                for (;;) {
                    std::function<void(std::unordered_map< unsigned int, std::vector<unsigned int> > &)> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task(this->worker_maps[i]);
                }
            }
        );
    }
}

ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(std::thread &worker: workers)
        worker.join();
}

// template <class F>
// void ThreadPool::enqueue(F&& f) {
//     {
//         std::unique_lock<std::mutex> lock(queue_mutex);
//         tasks.emplace(std::forward<F>(f));
//     }
//     condition.notify_one();
// }


void ThreadPool::enqueue(const std::function<void(std::unordered_map< unsigned int, std::vector<unsigned int> >&  )>& task) {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        tasks.emplace(task);
    }
    condition.notify_one();
}
