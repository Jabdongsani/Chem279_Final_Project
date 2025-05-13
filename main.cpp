#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <complex>

#include "periodic_solid.hpp"


int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }
    std::cout << "Input Path: " << argv[1] << std::endl;

    // Read geometry and lattice
    auto [atoms, lattice] = read_input(argv[1]);
    std::cout << "Read " << atoms.size() << " atoms, lattice dimension " << lattice.dim << std::endl;
    for (size_t i = 0; i < atoms.size(); ++i) {
        std::cout << i << ": " << atoms[i].element << " ("
                  << atoms[i].x << ", " << atoms[i].y << ", " << atoms[i].z << ")" << std::endl;
    }

    // Build basis set
    auto basis_functions = build_basis_set(atoms);
    std::cout << "Total basis functions: " << basis_functions.size() << std::endl;
    print_basis_info(basis_functions);

    // Compute normalization constants
    compute_normalization_constants(basis_functions);

    // Set EHT parameters
    //double K_eht = (atoms[0].element == "C") ? 2.8 : 2.3; // 2.8 for graphene, 2.3 for silicon
    double K_eht = 2.8;
    double cutoff_radius = 9.0 * angstrom_to_bohr; // Angstrom

    // Compute band structure for periodic systems
    if (lattice.dim > 0) {
        auto [k_path, k_labels] = get_graphene_k_path();
        compute_band_structure(basis_functions, lattice, k_path, k_labels, K_eht);
    } else {
        // For molecules, compute as before
        arma::mat S;
        arma::cx_mat H;
        build_matrices(basis_functions, S, H, K_eht, cutoff_radius);
        arma::vec eigvals;
        arma::cx_mat eigvecs;
        arma::mat X_mat = solve_eigenvalue(H, S, eigvals, eigvecs);
        // Compute total energy (as in original code)
        int total_electrons = 0;
        for (const auto &atom : atoms) {
            total_electrons += (atom.element == "H") ? 1 : (atom.element == "C") ? 4 : (atom.element == "Si") ? 4 : 0;
        }
        if (total_electrons % 2 != 0) {
            throw std::runtime_error("Invalid: Total number of valence electrons must be even.");
        }
        int occupied_orbitals = total_electrons / 2;
        double E_total = 0.0;
        for (int i = 0; i < occupied_orbitals; ++i) {
            E_total += 2 * eigvals(i);
        }
        std::cout << "Molecular energy: " << E_total << " eV" << std::endl;
    }

    return 0;
}