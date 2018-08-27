// Petter Strandmark 2013.
#ifndef SPII_SYM_ILDL_CONVERSIONS_H
#define SPII_SYM_ILDL_CONVERSIONS_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <lilc_matrix.h>

namespace spii
{

void eigen_to_lilc(const Eigen::MatrixXd& A, lilc_matrix<double>* Alilc_input)
{
	minimum_core_assert(Alilc_input);
	lilc_matrix<double>& Alilc = *Alilc_input;

	auto m = A.rows();
	auto n = A.cols();
	minimum_core_assert(m == n);

	int count = 0;
	Alilc.resize(n, n);
	fill(Alilc.first.begin(), Alilc.first.end(), 0);
	for (int i = 0; i < n; i++) {

		for (int j = 0; j < n; j++) {

			if (j < i) {
				continue;
			}

			Alilc.m_idx[i].push_back(j);
			Alilc.m_x[i].push_back(A(i, j));
			if (i != j) 
				Alilc.list[j].push_back(i);
			count++;
		}
	}
	Alilc.nnz_count = count;
}

Eigen::MatrixXd lilc_to_eigen(const lilc_matrix<double>& Alilc, bool symmetric=false)
{
	Eigen::MatrixXd A{Alilc.n_rows(), Alilc.n_cols()};
	A.fill(0.0);

	for (int i = 0; i < Alilc.n_rows(); i++) {
		for (std::size_t ind = 0; ind < Alilc.m_idx.at(i).size(); ind++) {
			auto j = Alilc.m_idx.at(i)[ind];
			auto value = Alilc.m_x.at(i).at(ind);
			if (symmetric) {
				A(i, j) = value;
			}
			A(j, i) = value;
		}
	}
	return std::move(A);
}

Eigen::DiagonalMatrix<double, Eigen::Dynamic> diag_to_eigen(const block_diag_matrix<double>& Ablock)
{
	minimum_core_assert(Ablock.n_rows() == Ablock.n_cols());
	minimum_core_assert(Ablock.off_diag.empty());

	Eigen::DiagonalMatrix<double, Eigen::Dynamic> A{Ablock.n_rows()};
	A.setZero();

	for (std::size_t i = 0; i < Ablock.main_diag.size(); ++i) {
		A.diagonal()[i] = Ablock.main_diag[i];
	}

	return std::move(A);
}

Eigen::MatrixXd block_diag_to_eigen(const block_diag_matrix<double>& Ablock)
{
	Eigen::MatrixXd A{Ablock.n_rows(), Ablock.n_cols()};
	A.fill(0.0);

	for (std::size_t i = 0; i < Ablock.main_diag.size(); ++i) {
		A(i, i) = Ablock.main_diag[i];
	}

	for (const auto& i_and_value: Ablock.off_diag) {
		auto i = i_and_value.first;
		auto value = i_and_value.second;
		A(i+1, i) = value;
		A(i, i+1) = value;
	}

	return std::move(A);
}

class MyPermutation
	: public Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>
{
public:
	MyPermutation(const std::vector<int>& perm)
		: Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>(perm.size())
	{
		for (std::size_t i = 0; i < perm.size(); ++i) {
			m_indices(i) = perm[i];
		}
	}
};

// Solve A*x = b, where
//
// A = S.inverse() * (P * L * D * L.transpose() * P.transpose()) * S.inverse()
//
void solve_system_ildl_dense(const Eigen::MatrixXd& D,
                             const Eigen::MatrixXd& L,
                             const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& S,
                             const MyPermutation& P,
                             const Eigen::VectorXd& lhs,
                             Eigen::VectorXd* x_output)
{
	minimum_core_assert(x_output);
	auto& x = *x_output;
	auto n = lhs.rows();
	minimum_core_assert(D.rows() == n);
	minimum_core_assert(D.cols() == n);
	minimum_core_assert(L.rows() == n);
	minimum_core_assert(L.cols() == n);
	minimum_core_assert(S.rows() == n);
	minimum_core_assert(S.cols() == n);
	
	// TODO: make more efficient if needed.

	x = S * lhs;
	x = P.inverse() * x;
	x = L.lu().solve(x);
	x = D.lu().solve(x);
	x = L.transpose().lu().solve(x);
	x = P.inverse() * x;
	x = S * x;
}

}

#endif
