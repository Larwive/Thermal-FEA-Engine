#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>

#include "thermal_fea/linalg/Matrix.hpp"
#include "thermal_fea/linalg/Vector.hpp"
#include "thermal_fea/mesh/Element.hpp"
#include "thermal_fea/mesh/Mesh.hpp"
#include "thermal_fea/mesh/Node.hpp"
#include "thermal_fea/solver/assembly.hpp"
#include "thermal_fea/solver/conjugate_gradient.hpp"

namespace {

std::vector<thermal_fea::mesh::Node>
generate_nodes_perturbed_grid(std::size_t nx, std::size_t ny,
															double perturb) {
	std::vector<thermal_fea::mesh::Node> nodes;
	nodes.reserve(nx * ny);

	std::mt19937 rng(42);
	std::uniform_real_distribution<double> dist(-perturb, perturb);

	std::size_t id = 0;
	for (std::size_t i = 0; i < nx; ++i) {
		for (std::size_t j = 0; j < ny; ++j) {
			double x = static_cast<double>(i) / static_cast<double>(nx - 1) +
				     dist(rng);
			double y = static_cast<double>(j) / static_cast<double>(ny - 1) +
				     dist(rng);

			x = std::min(std::max(x, 0.0), 1.0);
			y = std::min(std::max(y, 0.0), 1.0);

			nodes.push_back(thermal_fea::mesh::Node{id, x, y});
			++id;
		}
	}

	return nodes;
}

std::vector<thermal_fea::mesh::Element>
generate_elements_grid(std::size_t nx, std::size_t ny) {
	std::vector<thermal_fea::mesh::Element> elements;
	elements.reserve((nx - 1) * (ny - 1) * 2);

	auto node_id = [ny](std::size_t i, std::size_t j) {
		return i * ny + j;
	};

	for (std::size_t i = 0; i + 1 < nx; ++i) {
		for (std::size_t j = 0; j + 1 < ny; ++j) {
			std::size_t n00 = node_id(i, j);
			std::size_t n10 = node_id(i + 1, j);
			std::size_t n01 = node_id(i, j + 1);
			std::size_t n11 = node_id(i + 1, j + 1);

			elements.push_back(thermal_fea::mesh::Element{{n00, n10, n11}});
			elements.push_back(thermal_fea::mesh::Element{{n00, n11, n01}});
		}
	}

	return elements;
}

} // namespace

int main() {
	constexpr double perturb = 0.01;
	constexpr double k = 1.0;
	constexpr double diag_shift = 1e-3;

	{
		constexpr std::size_t nx = 12;
		constexpr std::size_t ny = 12;

		auto nodes = generate_nodes_perturbed_grid(nx, ny, perturb);
		auto elements = generate_elements_grid(nx, ny);

		thermal_fea::mesh::Mesh mesh(nodes, elements);
		thermal_fea::linalg::Matrix A = thermal_fea::solver::K(mesh, k);
		for (std::size_t i = 0; i < A.rows(); ++i) {
			A(i, i) += diag_shift;
		}

		thermal_fea::linalg::Vector b(mesh.nb_nodes());
		for (std::size_t i = 0; i < b.size(); ++i) {
			b[i] = 1.0 + 0.05 * static_cast<double>(i % 19);
		}

		thermal_fea::linalg::Vector x_ref =
				thermal_fea::solver::conjugate_gradient(A, b);
		thermal_fea::linalg::Vector x_opt =
				thermal_fea::solver::conjugate_gradient_optimized(A, b);

		thermal_fea::linalg::Vector r_ref = b - A * x_ref;
		thermal_fea::linalg::Vector r_opt = b - A * x_opt;
		thermal_fea::linalg::Vector diff = x_opt - x_ref;

		std::cout << "CG validity check (small grid)\n";
		std::cout << "Nodes: " << mesh.nb_nodes() << "\n";
		std::cout << "Elements: " << mesh.nb_elements() << "\n";
		std::cout << "Solution diff: " << diff.norm() << "\n\n";
	}

	{
		constexpr std::size_t nx = 50;
		constexpr std::size_t ny = 50;

		auto nodes = generate_nodes_perturbed_grid(nx, ny, perturb);
		auto elements = generate_elements_grid(nx, ny);

		thermal_fea::mesh::Mesh mesh(nodes, elements);
		thermal_fea::linalg::Matrix A = thermal_fea::solver::K(mesh, k);
		for (std::size_t i = 0; i < A.rows(); ++i) {
			A(i, i) += diag_shift;
		}

		thermal_fea::linalg::Vector b(mesh.nb_nodes());
		for (std::size_t i = 0; i < b.size(); ++i) {
			b[i] = 1.0 + 0.05 * static_cast<double>(i % 19);
		}

		auto start_opt = std::chrono::steady_clock::now();
		thermal_fea::linalg::Vector x_opt =
				thermal_fea::solver::conjugate_gradient_optimized(A, b);
		auto end_opt = std::chrono::steady_clock::now();

		std::chrono::duration<double> elapsed_opt = end_opt - start_opt;

		thermal_fea::linalg::Vector r_opt = b - A * x_opt;

		std::cout << "CG performance (large grid)\n";
		std::cout << "Nodes: " << mesh.nb_nodes() << "\n";
		std::cout << "Elements: " << mesh.nb_elements() << "\n";
		std::cout << "Elapsed opt (s): " << elapsed_opt.count() << "\n";
	}

	return 0;
}
