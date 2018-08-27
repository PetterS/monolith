#include <set>
#include <thread>

#include <catch.hpp>

#include <minimum/core/concurrent.h>
using namespace minimum::core;
using namespace std;

class TestWorker : public ConcurrentWorker<int, string> {
   public:
	using ConcurrentWorker::ConcurrentWorker;

   protected:
	void process(int i) override { output_queue.emplace(to_string(i)); }
};

TEST_CASE("ConcurrentWorker") {
	TestWorker worker(2);
	set<string> expected;
	for (int i = 0; i < 100; ++i) {
		worker.emplace(i);
		expected.emplace(to_string(i));
	}

	set<string> output;
	while (true) {
		auto s = worker.possibly_get();
		if (s.has_value()) {
			output.emplace(move(*s));
		} else {
			this_thread::yield();
		}
		if (output.size() == 100) {
			break;
		}
	}
	CHECK(output == expected);
	CHECK_FALSE(worker.possibly_get().has_value());
	worker.stop();
	worker.emplace(1);
	CHECK_FALSE(worker.possibly_get().has_value());
	worker.stop();
}

class FailingWorker : public ConcurrentWorker<int, int> {
   public:
	using ConcurrentWorker::ConcurrentWorker;
	static atomic<bool> has_started_process;

   protected:
	void process(int i) override {
		has_started_process.store(true);
		throw runtime_error("Error");
	}
};
atomic<bool> FailingWorker::has_started_process;

TEST_CASE("ConcurrentWorker_exception") {
	FailingWorker::has_started_process.store(false);
	FailingWorker worker(2);
	worker.emplace(1);
	worker.emplace(2);
	worker.emplace(3);
	while (!FailingWorker::has_started_process.load()) {
		this_thread::yield();
	}
	CHECK_FALSE(worker.possibly_get().has_value());
	CHECK_FALSE(worker.possibly_get().has_value());
	CHECK_THROWS_AS(worker.stop(), std::runtime_error);
}

TEST_CASE("ConcurrentQueue_basic") {
	ConcurrentQueue<vector<int>> queue;
	queue.emplace(5, 0);
	queue.emplace(3, 0);
	CHECK(queue.possibly_get()->size() == 5);
	CHECK(queue.possibly_get()->size() == 3);
	CHECK(!queue.possibly_get().has_value());
}

TEST_CASE("ConcurrentQueue_notify") {
	ConcurrentQueue<int> queue;
	queue.emplace(5);
	queue.emplace(3);
	CHECK(*queue.get() == 5);
	CHECK(*queue.get() == 3);
}

TEST_CASE("ConcurrentQueue_multi_thread") {
	ConcurrentQueue<int> queue;
	std::thread worker([&]() {
		this_thread::sleep_for(10ms);
		queue.emplace(1);
		this_thread::sleep_for(10ms);
		queue.emplace(2);
		this_thread::sleep_for(10ms);
		queue.emplace(3);
	});
	CHECK(*queue.get() == 1);
	CHECK(*queue.get() == 2);
	CHECK(*queue.get() == 3);
	worker.join();
	CHECK_FALSE(queue.possibly_get().has_value());
}

TEST_CASE("ConcurrentQueue_close") {
	ConcurrentQueue<int> queue;
	std::thread worker([&]() {
		this_thread::sleep_for(10ms);
		queue.emplace(1);
		this_thread::sleep_for(10ms);
		queue.close();
		this_thread::sleep_for(10ms);
		queue.emplace(3);
	});
	CHECK(*queue.get() == 1);
	CHECK_FALSE(queue.get().has_value());
	worker.join();
	CHECK_FALSE(queue.get().has_value());
}
