#include "flamegraph.h"

#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#include <minimum/core/flamegraph.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/core/time.h>

class FlamegraphGlobalState {
   public:
	FlamegraphGlobalState() {}

	void enable() {
		if (!enabled) {
			enabled = true;
			root.start();
		}
	}

	bool is_enabled() const { return enabled; }

	void write_to_file(const string& filename) {
		if (is_enabled()) {
			auto data = render_information();
			ofstream fout(filename);
			fout << data;
		}
	}

	void start(const char* name) {
		Node* node = current->get_child(name);
		node->start();
		current = node;
	}

	void stop() {
		if (current != &root) {
			current->stop();
			current = current->parent;
		}
	}

	string render_information() {
		if (!enabled) {
			return "";
		}
		root.stop();

		ostringstream result;
		vector<Node*> call_stack;
		call_stack.push_back(&root);
		dfs(&root, call_stack, [&](const vector<Node*>& stack) {
			bool first = true;
			double self_time = stack.back()->total_time;
			for (auto node : stack.back()->children) {
				self_time -= node->total_time;
			}
			for (auto node : stack) {
				if (!first) {
					result << ";";
				}
				first = false;
				result << node->name;
			}
			result << ' ' << (long long)(1e6 * self_time + 0.5) << '\n';
		});

		root.start();
		return result.str();
	}

	string render_pretty_string() {
		if (!enabled) {
			return "";
		}
		root.stop();

		ostringstream result;
		vector<Node*> call_stack;
		call_stack.push_back(&root);
		dfs(&root, call_stack, [&](const vector<Node*>& stack) {
			for (auto i : minimum::core::range(stack.size() - 1)) {
				result << "            ";
			}
			auto t = stack.back()->total_time;
			auto h = lround(floor(t / 3600));
			t -= h * 3600;
			auto m = lround(floor(t / 60));
			t -= m * 60;
			auto s = lround(floor(t));
			t -= s;
			auto frac = lround(floor(t * 100));
			result << setw(2) << setfill('0') << h << ':';
			result << setw(2) << setfill('0') << m << ':';
			result << setw(2) << setfill('0') << s << '.';
			result << setw(2) << setfill('0') << frac;
			result << " " << stack.back()->name << '\n';
		});

		root.start();
		return result.str();
	}

   private:
	struct Node {
		Node() : name("<global>"), parent(nullptr) {}

		~Node() {
			for (auto child : children) {
				delete child;
			}
		}

		void start() { start_time = minimum::core::wall_time(); }

		void stop() { total_time += minimum::core::wall_time() - start_time; }

		Node* get_child(const char* name) {
			Node* node = nullptr;
			for (auto child : children) {
				// Since we are deling with string constants known at
				// compile-time, we can get away with not using strcmp.
				if (name == child->name) {
					node = child;
					break;
				}
			}
			if (node == nullptr) {
				add_child(name);
				node = children.back();
			}
			return node;
		}

		double total_time = 0;
		const char* name;
		Node* parent;
		vector<Node*> children;

	   private:
		Node(const char* name_, Node* parent_) : name(name_), parent(parent_) {}

		void add_child(const char* name) { children.emplace_back(new Node(name, this)); }

		double start_time;
	};

	void dfs(Node* root,
	         vector<Node*>& call_stack,
	         const function<void(const vector<Node*>&)>& action) {
		action(call_stack);
		for (auto child : root->children) {
			call_stack.emplace_back(child);
			dfs(child, call_stack, action);
			call_stack.pop_back();
		}
	}

	bool enabled = false;
	Node root;
	Node* current = &root;
};

thread_local FlamegraphGlobalState flamegraph_state;

bool flamegraph::internal::is_enabled() { return flamegraph_state.is_enabled(); }

void flamegraph::internal::start(const char* name) { flamegraph_state.start(name); }

void flamegraph::internal::stop() { flamegraph_state.stop(); }

void flamegraph::enable() { flamegraph_state.enable(); }

void flamegraph::write_to_file(const std::string& filename) {
	flamegraph_state.write_to_file(filename);
}

void flamegraph::write_pretty(std::ostream& out) { out << flamegraph_state.render_pretty_string(); }

string flamegraph::internal::render_information() { return flamegraph_state.render_information(); }
