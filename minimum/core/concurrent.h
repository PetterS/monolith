#pragma once

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <optional>
#include <queue>
#include <thread>

namespace minimum {
namespace core {

template <typename T>
class ConcurrentQueue {
   public:
	// Adds an element to the queue and notifies a waiting thread, if
	// available.
	template <class... S>
	void emplace(S&&... s) {
		{
			std::lock_guard<std::mutex> lock(m);
			if (closed) {
				return;
			}
			queue.emplace(std::forward<S>(s)...);
		}
		cond.notify_one();
	}

	// Closes the queue and no further processing will take place.
	void close() {
		{
			std::lock_guard<std::mutex> lock(m);
			closed = true;
		}
		cond.notify_all();
	}

	// Does not wait for an element to become available.
	// Returns empty iff no element is available or the queue
	// is closed.
	std::optional<T> possibly_get() {
		std::lock_guard<std::mutex> lock(m);
		if (queue.empty() || closed) {
			return std::nullopt;
		}
		auto elem = std::move(queue.front());
		queue.pop();
		return elem;
	}

	// Blocks until an element becomes available.
	// Returns empty iff the queue is closed in the meantime.
	std::optional<T> get() {
		std::unique_lock<std::mutex> lock(m);
		cond.wait(lock, [&]() { return !queue.empty() || closed; });
		if (closed) {
			return std::nullopt;
		}

		auto elem = std::move(queue.front());
		queue.pop();
		return elem;
	}

   private:
	std::queue<T> queue;
	bool closed = false;

	std::mutex m;
	std::condition_variable cond;
};

template <typename Input, typename Output>
class ConcurrentWorker {
   public:
	ConcurrentWorker(int num_threads = 1);
	// Will block if the worker is not already stopped.
	~ConcurrentWorker();

	// Closes the queues and waits for the threads to finish.
	void stop();

	template <class... S>
	void emplace(S&&... s) {
		input_queue.emplace(std::forward<S>(s)...);
	}

	// Does not wait for an element to become available.
	// Returns empty iff no element is available or the worker
	// is stopped.
	std::optional<Output> possibly_get() { return output_queue.possibly_get(); }

   protected:
	// Overload this to add zero or more outputs to the output
	// queue for every input.
	virtual void process(Input input) = 0;

	ConcurrentQueue<Output> output_queue;

   private:
	void thread_function();

	ConcurrentQueue<Input> input_queue;
	std::vector<std::thread> threads;
	std::atomic_flag stopped = ATOMIC_FLAG_INIT;

	std::mutex exception_mutex;
	std::exception_ptr exception;
};

template <typename Input, typename Output>
ConcurrentWorker<Input, Output>::ConcurrentWorker(int num_threads) {
	auto to_run = [this]() { thread_function(); };
	for (int i = 0; i < num_threads; ++i) {
		threads.emplace_back(to_run);
	}
}

template <typename Input, typename Output>
ConcurrentWorker<Input, Output>::~ConcurrentWorker() {
	stop();
}

template <typename Input, typename Output>
void ConcurrentWorker<Input, Output>::stop() {
	if (stopped.test_and_set()) {
		return;
	}
	input_queue.close();
	output_queue.close();
	for (auto& t : threads) {
		t.join();
	}
	if (exception) {
		std::rethrow_exception(exception);
	}
}

template <typename Input, typename Output>
void ConcurrentWorker<Input, Output>::thread_function() {
	try {
		while (true) {
			// Block waiting for data or a closed queue.
			auto data = input_queue.get();
			if (!data.has_value()) {
				return;
			}
			process(std::move(*data));
		}
	} catch (...) {
		std::lock_guard<std::mutex> lock(exception_mutex);
		exception = std::current_exception();
	}
}
}  // namespace core
}  // namespace minimum
