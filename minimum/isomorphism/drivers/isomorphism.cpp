
#include <iostream>
#include <random>

#include <minimum/core/time.h>

#include <minimum/isomorphism/graph-isomorphism.h>
using namespace minimum::isomorphism;

int main_function(const std::vector<std::string>& args) {
	using namespace std;

	string filter;
	if (args.size() > 0) {
		filter = args[0];
	}

	vector<string> file_names;
	ifstream all_files("database/all_files.txt");
	string file_name;
	while (all_files >> file_name) {
		file_names.push_back(move(file_name));
	}

	vector<string> gtr_file_names;
	ifstream gtr_files("database/all_gtr.txt");
	while (gtr_files >> file_name) {
		gtr_file_names.push_back(move(file_name));
	}

	map<string, int> ground_truth;
	for (auto gtr_file_name : gtr_file_names) {
		ifstream gtr_file("database/" + gtr_file_name);
		while (gtr_file >> file_name) {
			int num_isomorphisms = -1;
			gtr_file >> num_isomorphisms;
			ground_truth[file_name] = num_isomorphisms;
		}
	}

	// mt19937_64 engine(time(nullptr));
	mt19937_64 engine(101u);
	shuffle(file_names.begin(), file_names.end(), engine);

	ofstream result_file(filter + "results.txt");

	for (string file_name_1 : file_names) {
		if (file_name_1.find(filter) == string::npos) continue;
		// if (c++ > 100) break;

		auto file_name_2 = file_name_1;
		auto itr = file_name_2.find(".A");
		if (itr != string::npos) {
			itr++;
			file_name_2[itr] = 'B';
		}

		Eigen::SparseMatrix<size_t> graph1;
		read_graph("database/" + file_name_1, &graph1);
		Eigen::SparseMatrix<size_t> graph2;
		read_graph("database/" + file_name_2, &graph2);

		// if (graph1.rows() > 400) continue;

		cout << file_name_1 << " <--> " << file_name_2 << endl;
		cout << "-- n=" << graph1.rows()
		     << " avg. deg=" << double(graph1.nonZeros()) / graph1.rows() << endl;
		itr = file_name_1.find_last_of("\\/") + 1;
		auto base_file_name = file_name_1.substr(itr);
		cout << "-- Ground truth: " << ground_truth[base_file_name] << endl;

		try {
			double start_time = minimum::core::wall_time();
			int num_isomorphisms = is_isomorphic(graph1, graph2, nullptr);
			double total_time = minimum::core::wall_time() - start_time;
			if (num_isomorphisms >= 1) {
				// vector<int> isomorphism;
				// is_isomorphic(graph1, graph2, &isomorphism);
				// cout << "-- [";
				// for (auto e: isomorphism) {
				//	cout << e << " ";
				//}
				// cout << "]" << endl;
				cout << "-- Total number of isomorphisms: " << num_isomorphisms << " in "
				     << total_time << " seconds." << endl;
				if (num_isomorphisms != ground_truth[base_file_name]) {
					throw logic_error("Incorrect number of isomorphisms.");
				}
			} else {
				cout << "-- No isomorphisms found in " << total_time << " seconds." << endl;
				throw logic_error("At least one isomorphism is expected.");
			}

			result_file << base_file_name << " " << graph1.rows() << " " << total_time << endl;
		} catch (exception&) {
			cout << endl << endl;
			graph_to_unit_test(graph1, "graph1");
			graph_to_unit_test(graph2, "graph2");
			cout << endl;
			throw;
		}
	}

	return 0;
}

int main(int arg_count, char* arrayargs[]) {
	using namespace std;
	try {
		std::vector<string> args(arrayargs + 1, arrayargs + arg_count);
		return main_function(args);
	} catch (exception& err) {
		cerr << endl << "ERROR: " << err.what() << endl;
	}
	return 1;
}
