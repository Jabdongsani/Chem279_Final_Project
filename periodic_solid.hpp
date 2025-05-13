#ifndef PERIODIC_SOLID_HPP
#define PERIODIC_SOLID_HPP

#include <string>
#include <vector>
#include <array>
#include <armadillo>


const double angstrom_to_bohr = 1.8897259886;
const double cutoff_radius_bohr = 9.0 * angstrom_to_bohr;
// const double k_eht = 5.0;

struct Atom {
    std::string element;
    double x, y, z;
};

struct Lattice {
    std::array<double, 3> a1, a2, a3; // Lattice vectors
    int dim; // 1D, 2D, or 3D
};

struct val_exp_coeff {
    arma::vec exps;
    arma::vec coeffs;
    double onsite_energy; // E_mu_mu for EHT
};

struct basis_function {
    double x, y, z; // Atom coordinates
    int l, m, n; // Angular momentum exponents
    std::string orbital_type; // e.g., "2s", "2p", "3d"
    arma::vec exps; // Gaussian exponents
    arma::vec coeffs; // Contraction coefficients
    arma::vec norms; // Normalization constants
    double onsite_energy; // EHT diagonal element
    std::string element; // Associated element
};



std::pair<std::vector<Atom>, Lattice> read_input(std::string file_path);

val_exp_coeff get_valence_data(std::string element, std::string orbital_type);

std::vector<std::array<int, 3>> cartesian_functions(std::string orbital_type);

std::vector<basis_function> build_basis_set(std::vector<Atom> atoms);

void print_basis_info(const std::vector<basis_function> &basis_functions);

int calc_double_factorial(int n);

int binomial_coeff(int n, int k);

double compute_primitive_normalization(int l, int m, int n, double zeta);

void compute_normalization_constants(std::vector<basis_function> &basis_functions);

std::array<double, 3> compute_product_center(const std::array<double, 3> &R_A, double alpha,
                                             const std::array<double, 3> &R_B, double beta);

std::pair<double, double> compute_exp_prefactor(double alpha, double beta, double R_A, double R_B);

double compute_double_summation(int l_A, int l_B, double R_P, double R_A, double R_B, double alpha, double beta); 

double primitive_overlap_3D(const basis_function &bfA, int exp_a, const basis_function &bfB, int exp_b); 

double contracted_overlap(const basis_function &bfA, const basis_function &bfB, const std::array<double, 3> &R_shift = {0.0, 0.0, 0.0}, double cutoff_radius = cutoff_radius_bohr);

void build_matrices(const std::vector<basis_function> &basis,
                    arma::mat &S, arma::cx_mat &H,
                    double K_eht = 2.8, double cutoff_radius = cutoff_radius_bohr);

void apply_bloch_phase(const std::vector<basis_function> &basis,
                       const Lattice &lattice,
                       const std::array<double, 3> &k_point,
                       arma::cx_mat &S_k, arma::cx_mat &H_k,
                       double K_eht = 5.0, double cutoff_radius = cutoff_radius_bohr);

arma::vec solve_k_eigenvalue(const arma::cx_mat &H_k, const arma::cx_mat &S_k);

void compute_band_structure(const std::vector<basis_function> &basis,
                            const Lattice &lattice,
                            const std::vector<std::array<double, 3>> &k_path,
                            const std::vector<std::string> &k_labels,
                            double K_eht = 5.0);

std::pair<std::vector<std::array<double, 3>>, std::vector<std::string>>
get_graphene_k_path(double a_cc, std::vector<double>& k_dist, std::vector<double>& k_ticks);

std::pair<std::vector<std::array<double, 3>>, std::vector<std::string>> get_graphene_k_path(double a_cc = 1.44 * angstrom_to_bohr);

arma::mat solve_eigenvalue(arma::cx_mat &H,arma::mat &S, arma::vec &eigvals, arma::cx_mat &eigvecs);


#endif