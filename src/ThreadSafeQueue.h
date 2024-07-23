#include <queue>
#include <mutex>
#include <condition_variable>

template <typename T>
class ThreadSafeQueue {
public:
    void push(const T& value) {
        std::lock_guard<std::mutex> lock(mtx);
        q.push(value);
        cv.notify_one();
    }

    bool try_pop(T& result) {
        std::lock_guard<std::mutex> lock(mtx);
        if (q.empty()) {
            return false;
        }
        result = q.front();
        q.pop();
        return true;
    }

private:
    std::queue<T> q;
    mutable std::mutex mtx;
    std::condition_variable cv;
};