#pragma once

namespace minimum {
namespace algorithms {
class DisjointSetElement {
   public:
	DisjointSetElement() : parent(this) {}
	DisjointSetElement(const DisjointSetElement&) = delete;
	DisjointSetElement(DisjointSetElement&&) = delete;

	DisjointSetElement* find() {
		// https://en.wikipedia.org/wiki/Disjoint-set_data_structure
		if (this->parent != this) {
			this->parent = this->parent->find();
		}
		return this->parent;
	}

	bool join(DisjointSetElement* that) {
		// https://en.wikipedia.org/wiki/Disjoint-set_data_structure
		auto this_root = find();
		auto that_root = that->find();
		// if x and y are already in the same set (i.e., have the same root or representative)
		if (this_root == that_root) {
			return false;
		}
		// x and y are not in same set, so we merge them
		if (this_root->rank < that_root->rank) {
			this_root->parent = that_root;
		} else if (this_root->rank > that_root->rank) {
			that_root->parent = this_root;
		} else {
			that_root->parent = this_root;
			this_root->rank += 1;
		}
		return true;
	}

   private:
	DisjointSetElement* parent;
	int rank = 0;
};
}  // namespace algorithms
}  // namespace minimum
