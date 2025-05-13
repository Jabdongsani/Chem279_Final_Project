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



// Read input file with lattice vectors and atomic coordinates
std::pair<std::vector<Atom>, Lattice> read_input(std::string file_path) {
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        throw std::runtime_error("File path does not exist!");
    }

    int atom_count, dim;
    infile >> atom_count >> dim;

    Lattice lattice;
    lattice.dim = dim;

    // Read and convert lattice vectors to bohr
    infile >> lattice.a1[0] >> lattice.a1[1] >> lattice.a1[2];
    lattice.a1[0] *= angstrom_to_bohr;
    lattice.a1[1] *= angstrom_to_bohr;
    lattice.a1[2] *= angstrom_to_bohr;

    if (dim >= 2) {
        infile >> lattice.a2[0] >> lattice.a2[1] >> lattice.a2[2];
        lattice.a2[0] *= angstrom_to_bohr;
        lattice.a2[1] *= angstrom_to_bohr;
        lattice.a2[2] *= angstrom_to_bohr;
    } else {
        lattice.a2 = {0.0, 0.0, 0.0};
    }

    if (dim == 3) {
        infile >> lattice.a3[0] >> lattice.a3[1] >> lattice.a3[2];
        lattice.a3[0] *= angstrom_to_bohr;
        lattice.a3[1] *= angstrom_to_bohr;
        lattice.a3[2] *= angstrom_to_bohr;
    } else {
        lattice.a3 = {0.0, 0.0, 0.0};
    }

    std::vector<Atom> atoms;
    for (int i = 0; i < atom_count; ++i) {
        int element_int;
        Atom atom;
        infile >> element_int >> atom.x >> atom.y >> atom.z;

        // Convert position from angstrom to bohr
        atom.x *= angstrom_to_bohr;
        atom.y *= angstrom_to_bohr;
        atom.z *= angstrom_to_bohr;

        atom.element = (element_int == 1) ? "H" :
                       (element_int == 6) ? "C" :
                       (element_int == 14) ? "Si" : "Unknown";

        atoms.push_back(atom);
    }

    return {atoms, lattice};
}


// Get STO basis parameters for an element and orbital
val_exp_coeff get_valence_data(std::string element, std::string orbital_type) {
    val_exp_coeff data;
    // divide by 3.751
    if (element == "H") {
        if (orbital_type == "1s") {
            data.exps = {3.42525091, 0.62391373, 0.16885540};
            data.coeffs = {0.15432897, 0.53532814, 0.44463454};
            data.onsite_energy = -13.6; // eV
        }
    } else if (element == "C") {
        if (orbital_type == "2s") {
            data.exps = {2.037}; // C-sp from Table I
            data.coeffs = {0.741};
            data.onsite_energy = -20.316;
        } else if (orbital_type == "2p") {
            data.exps = {0.88, 1.3};
            data.coeffs = {0.640, 0.412};
            data.onsite_energy = -13.670;
            /* 
            data.exps = {1.777, 3.249};
            data.coeffs = {0.640, 0.412};
            data.onsite_energy = -13.670;
            */
        }
    } else if (element == "Si") {
        if (orbital_type == "3s") {
            data.exps = {1.634};
            data.coeffs = {0.750};
            data.onsite_energy = -17.3;
        } else if (orbital_type == "3p") {
            data.exps = {1.428, 2.500};
            data.coeffs = {0.680, 0.320};
            data.onsite_energy = -9.2;
        } else if (orbital_type == "3d") {
            data.exps = {1.100};
            data.coeffs = {0.500};
            data.onsite_energy = -1.0;
        }
    }
    return data;
}

// Define Cartesian functions for orbitals
std::vector<std::array<int, 3>> cartesian_functions(std::string orbital_type) {
    std::vector<std::array<int, 3>> shell_functions;
    if (orbital_type == "1s" || orbital_type == "2s" || orbital_type == "3s") {
        shell_functions.push_back({0, 0, 0});
    } else if (orbital_type == "2p" || orbital_type == "3p") {
        shell_functions.push_back({1, 0, 0}); // px
        shell_functions.push_back({0, 1, 0}); // py
        shell_functions.push_back({0, 0, 1}); // pz
    } else if (orbital_type == "3d") {
        shell_functions.push_back({2, 0, 0}); // dxx
        shell_functions.push_back({0, 2, 0}); // dyy
        shell_functions.push_back({0, 0, 2}); // dzz
        shell_functions.push_back({1, 1, 0}); // dxy
        shell_functions.push_back({1, 0, 1}); // dxz
    }
    return shell_functions;
}

// Build basis set for the molecule
std::vector<basis_function> build_basis_set(std::vector<Atom> atoms) {
    std::vector<basis_function> basis_functions;
    for (const auto &atom : atoms) {
        std::vector<std::string> orbital_types;
        if (atom.element == "H") {
            orbital_types = {"1s"};
        } else if (atom.element == "C") {
            orbital_types = {"2s", "2p"};
        } else if (atom.element == "Si") {
            orbital_types = {"3s", "3p", "3d"};
        }

        for (const auto &orb : orbital_types) {
            val_exp_coeff data = get_valence_data(atom.element, orb);
            auto shell_functions = cartesian_functions(orb);
            for (const auto &sf : shell_functions) {
                basis_function bf;

                bf.x = atom.x;
                bf.y = atom.y;
                bf.z = atom.z;

                bf.l = sf[0];
                bf.m = sf[1];
                bf.n = sf[2];
                bf.orbital_type = orb;
                bf.exps = data.exps;
                bf.coeffs = data.coeffs;
                bf.onsite_energy = data.onsite_energy;
                bf.element = atom.element;
                basis_functions.push_back(bf);
            }
        }
    }
    return basis_functions;
}

// Print basis set information
void print_basis_info(const std::vector<basis_function> &basis_functions) {
    std::cout << "Basis set info:" << std::endl;
    for (size_t i = 0; i < basis_functions.size(); ++i) {
        const auto &bf = basis_functions[i];
        std::cout << "Basis function " << i << ": "
                  << bf.element << ", " << bf.orbital_type << ", "
                  << "center=(" << bf.x << "," << bf.y << "," << bf.z << "), "
                  << "(l,m,n)=(" << bf.l << "," << bf.m << "," << bf.n << ")" << std::endl;
        std::cout << "   exps = " << bf.exps.t();
        std::cout << "   coeffs = " << bf.coeffs.t();
        std::cout << "   onsite_energy = " << bf.onsite_energy << " eV" << std::endl;
    }
}

// Compute double factorial
int calc_double_factorial(int n) {
    if (n <= 0) return 1;
    int result = 1;
    for (int i = n; i > 0; i -= 2) result *= i;
    return result;
}

// Compute binomial coefficient
int binomial_coeff(int n, int k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    int result = 1;
    for (int i = 0; i < k; ++i) {
        result *= (n - i);
        result /= (i + 1);
    }
    return result;
}

// Compute normalization constant for STO (simplified for single exponent)
double compute_primitive_normalization(int l, int m, int n, double zeta) {
    // For STO: N = (2*zeta)^(n+1/2) / sqrt((2n)!), where n = l + m + n
    int n_total = l + m + n;
    double norm = std::pow(2.0 * zeta, n_total + 0.5) / std::sqrt(calc_double_factorial(2 * n_total - 1));
    return norm;
}

// Compute normalization constants for all basis functions
void compute_normalization_constants(std::vector<basis_function> &basis_functions) {
    for (auto &bf : basis_functions) {
        bf.norms.set_size(bf.exps.size());
        for (size_t k = 0; k < bf.exps.size(); ++k) {
            bf.norms[k] = compute_primitive_normalization(bf.l, bf.m, bf.n, bf.exps[k]);
        }
    }
}

// Compute product Gaussian center
std::array<double, 3> compute_product_center(const std::array<double, 3> &R_A, double alpha,
                                             const std::array<double, 3> &R_B, double beta) {
    std::array<double, 3> R_P;
    double denom = alpha + beta;
    for (int i = 0; i < 3; ++i) {
        R_P[i] = (alpha * R_A[i] + beta * R_B[i]) / denom;
    }
    return R_P;
}

// Compute exponential prefactor
std::pair<double, double> compute_exp_prefactor(double alpha, double beta, double R_A, double R_B) {
    double denom = alpha + beta;
    double exp_pref = std::exp(-alpha * beta * std::pow(R_A - R_B, 2) / denom);
    double sqrt_term = std::sqrt(M_PI / denom);
    return {exp_pref, sqrt_term};
}

// Compute double summation for overlap
double compute_double_summation(int l_A, int l_B, double R_P, double R_A, double R_B, double alpha, double beta) {
    auto prefactor_sqrt = compute_exp_prefactor(alpha, beta, R_A, R_B);
    double exp_pref = prefactor_sqrt.first;
    double sqrt_term = prefactor_sqrt.second;
    double sum = 0.0;
    for (int i = 0; i <= l_A; ++i) {
        for (int j = 0; j <= l_B; ++j) {
            if ((i + j) % 2 == 0) {
                double double_fact = calc_double_factorial(i + j - 1);
                double num = double_fact * std::pow(R_P - R_A, l_A - i) * std::pow(R_P - R_B, l_B - j);
                double denom = std::pow(2.0 * (alpha + beta), (i + j) / 2.0);
                int binom_A = binomial_coeff(l_A, i);
                int binom_B = binomial_coeff(l_B, j);
                sum += binom_A * binom_B * (num / denom);
            }
        }
    }
    return exp_pref * sqrt_term * sum;
}

// Compute 3D primitive overlap
double primitive_overlap_3D(const basis_function &bfA, int exp_a, const basis_function &bfB, int exp_b) {
    double alpha = bfA.exps[exp_a];
    double beta = bfB.exps[exp_b];
    std::array<double, 3> RA = {bfA.x, bfA.y, bfA.z};
    std::array<double, 3> RB = {bfB.x, bfB.y, bfB.z};
    auto RP = compute_product_center(RA, alpha, RB, beta);
    double S_x = compute_double_summation(bfA.l, bfB.l, RP[0], RA[0], RB[0], alpha, beta);
    double S_y = compute_double_summation(bfA.m, bfB.m, RP[1], RA[1], RB[1], alpha, beta);
    double S_z = compute_double_summation(bfA.n, bfB.n, RP[2], RA[2], RB[2], alpha, beta);
    return S_x * S_y * S_z;
}

// Compute contracted overlap
double contracted_overlap(const basis_function &bfA, const basis_function &bfB, const std::array<double, 3> &R_shift, double cutoff_radius) {
    // Check distance between centers
    double dx = bfA.x - (bfB.x + R_shift[0]);
    double dy = bfA.y - (bfB.y + R_shift[1]);
    double dz = bfA.z - (bfB.z + R_shift[2]);
    double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (distance > cutoff_radius) return 0.0; // Apply cutoff

    double S = 0.0;
    for (size_t a = 0; a < bfA.exps.size(); ++a) {
        for (size_t b = 0; b < bfB.exps.size(); ++b) {
            std::array<double, 3> RA = {bfA.x, bfA.y, bfA.z};
            std::array<double, 3> RB = {bfB.x + R_shift[0], bfB.y + R_shift[1], bfB.z + R_shift[2]};
            auto RP = compute_product_center(RA, bfA.exps[a], RB, bfB.exps[b]);
            double Sx = compute_double_summation(bfA.l, bfB.l, RP[0], RA[0], RB[0], bfA.exps[a], bfB.exps[b]);
            double Sy = compute_double_summation(bfA.m, bfB.m, RP[1], RA[1], RB[1], bfA.exps[a], bfB.exps[b]);
            double Sz = compute_double_summation(bfA.n, bfB.n, RP[2], RA[2], RB[2], bfA.exps[a], bfB.exps[b]);
            double overlap = Sx * Sy * Sz;
            S += bfA.coeffs[a] * bfB.coeffs[b] * bfA.norms[a] * bfB.norms[b] * overlap;
        }
    }
    return S;
}

// Build real-space overlap and Hamiltonian matrices
void build_matrices(const std::vector<basis_function> &basis,
                    arma::mat &S, arma::cx_mat &H,
                    double K_eht, double cutoff_radius) {
    int N = basis.size();
    S.zeros(N, N);
    H.zeros(N, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j <= i; ++j) {
            double s_ij = contracted_overlap(basis[i], basis[j], {0.0, 0.0, 0.0}, cutoff_radius);
            S(i, j) = S(j, i) = s_ij;
            if (i == j) {
                H(i, i) = basis[i].onsite_energy; // diagonal은 emprical onsite energy
            } else {
                double h_ij = 0.5 * K_eht * (basis[i].onsite_energy + basis[j].onsite_energy) * s_ij;
                H(i, j) = H(j, i) = h_ij;
            }
        }
    }
}

// Apply Bloch phase for periodic systems
void apply_bloch_phase(const std::vector<basis_function> &basis,
                       const Lattice &lattice,
                       const std::array<double, 3> &k_point,
                       arma::cx_mat &S_k, arma::cx_mat &H_k,
                       double K_eht, double cutoff_radius) {
    int N = basis.size(); // 우리 경우에는 8
    // overlap과 hamiltonian matrix 0으로 reset
    S_k.zeros(N, N);
    H_k.zeros(N, N);

    int n_cells = 1; // number of neighbor cells in each direction

    for (int n1 = -n_cells; n1 <= n_cells; ++n1) { // 격자의 주기성 고려, 나부터 양쪽으로 n_cells 만큼
        // y dimension 있으면 또 똑같이 n_cells 만큼 
        for (int n2 = (lattice.dim >= 2 ? -n_cells : 0); n2 <= (lattice.dim >= 2 ? n_cells : 0); ++n2) {
            // z dimension 있으면 또 똑같이 n_cells 만큼
            for (int n3 = (lattice.dim == 3 ? -n_cells : 0); n3 <= (lattice.dim == 3 ? n_cells : 0); ++n3) {

                std::array<double, 3> R_m = {
                    n1 * lattice.a1[0] + n2 * lattice.a2[0] + n3 * lattice.a3[0],
                    n1 * lattice.a1[1] + n2 * lattice.a2[1] + n3 * lattice.a3[1],
                    n1 * lattice.a1[2] + n2 * lattice.a2[2] + n3 * lattice.a3[2]
                };

                for (int i = 0; i < N; ++i) {
                    std::array<double, 3> R_i = {basis[i].x, basis[i].y, basis[i].z};

                    for (int j = 0; j < N; ++j) {
                        std::array<double, 3> R_j = {
                            basis[j].x + R_m[0],
                            basis[j].y + R_m[1],
                            basis[j].z + R_m[2]
                        };

                        // phase = k · (R_i - R_j)
                        double phase = k_point[0] * (R_i[0] - R_j[0]) +
                                       k_point[1] * (R_i[1] - R_j[1]) +
                                       k_point[2] * (R_i[2] - R_j[2]);
                        if (i == 0 && j == 0) {
                        // std::cout << "k = (" << k_point[0] << ", " << k_point[1] << ", " << k_point[2] << "), ";
                        // std::cout << "phase = " << phase << std::endl;
}

                        // std::exp(std::complex<double>(0.0, phase)) 
                        // → std::complex<double>(cos(phase), sin(phase)) 
                        std::complex<double> bloch_factor = std::exp(std::complex<double>(0.0, phase)); // 여기가 i 부분

                        double s_ij = contracted_overlap(basis[j], basis[i], R_m, cutoff_radius); // overlap을 구함
                        // bloch phase를 곱해줌
                        S_k(j, i) += bloch_factor * s_ij;
                        // if (std::abs(s_ij) > 1e-4 && i != j && n1 == 0 && n2 == 0 && n3 == 0) {
                        // std::cout << "[s_ij] i=" << i << ", j=" << j
                        // << ", s_ij=" << s_ij
                        // << ", h_ij=" << 0.5 * K_eht * (basis[i].onsite_energy + basis[j].onsite_energy) * s_ij
                        // << std::endl;

                        if (i == j && n1 == 0 && n2 == 0 && n3 == 0) {
                            H_k(i, i) = basis[i].onsite_energy;
                        } else {
                            double h_ij = 0.5 * K_eht * (basis[i].onsite_energy + basis[j].onsite_energy) * s_ij;
                            H_k(j, i) += bloch_factor * h_ij;
                        }
                        // if (i != j && std::abs(H_k(j, i)) > 1e-3 && n1 == 0 && n2 == 0 && n3 == 0) {
                        // std::cout << "[H_k] j=" << j << ", i=" << i
                        // << ", H_k=" << H_k(j, i) << std::endl;

                    }
                }
            }
        }
    }
}



// Solve generalized eigenvalue problem for k-point
arma::vec solve_k_eigenvalue(const arma::cx_mat &H_k, const arma::cx_mat &S_k) {
    arma::vec eigvals;
    arma::cx_mat eigvecs;

    arma::mat S_real = arma::real(S_k);

    arma::cx_mat H_copy = H_k;

    solve_eigenvalue(const_cast<arma::cx_mat&>(H_k), S_real, eigvals, eigvecs);

    return eigvals;
}

// Compute band structure along a k-path
void compute_band_structure(const std::vector<basis_function> &basis,
                            const Lattice &lattice,
                            const std::vector<std::array<double, 3>> &k_path,
                            const std::vector<std::string> &k_labels,
                            double K_eht) {
    std::ofstream outfile("band_structure.dat");
    int N = basis.size();
    // double fermi_shift = -13.7;  // fermi level shift (eV)
    arma::cx_mat S_k(N, N), H_k(N, N);
    for (size_t i = 0; i < k_path.size(); ++i) {
        apply_bloch_phase(basis, lattice, k_path[i], S_k, H_k, K_eht);
        arma::vec energies = solve_k_eigenvalue(H_k, S_k);

        outfile << i << " ";
        for (const auto &e : energies) {
            outfile << e << " ";
        }
        outfile << std::endl;
    }
    outfile.close();
    std::cout << "Band structure written to band_structure.dat" << std::endl;
}

// Define k-path for graphene (M -> Gamma -> K -> M)
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <utility>

std::pair<std::vector<std::array<double, 3>>, std::vector<std::string>>
get_graphene_k_path(double a_cc, std::vector<double>& k_dist, std::vector<double>& k_ticks) {
    // 1. Real-space lattice vectors for graphene
    std::array<double, 3> a1 = {a_cc * 1.5,  a_cc * std::sqrt(3) / 2, 0.0};
    std::array<double, 3> a2 = {a_cc * 1.5, -a_cc * std::sqrt(3) / 2, 0.0};

    // 2. Reciprocal lattice vectors using cross product
    double volume = a1[0] * a2[1] - a1[1] * a2[0];  // 2D cross product z-component
    std::array<double, 3> b1 = {
        2 * M_PI *  a2[1] / volume,
        -2 * M_PI * a2[0] / volume,
        0.0
    };
    std::array<double, 3> b2 = {
        -2 * M_PI * a1[1] / volume,
        2 * M_PI *  a1[0] / volume,
        0.0
    };

    // 3. High-symmetry points
    std::array<double, 3> Gamma = {0.0, 0.0, 0.0};
    std::array<double, 3> K = {
        (b1[0] + b2[0]) / 3,
        (b1[1] + b2[1]) / 3,
        0.0
    };
    std::array<double, 3> M = {
        (b1[0]) / 2,
        (b1[1]) / 2,
        0.0
    };

    // 4. Interpolation helper
    auto interpolate = [](const std::array<double, 3>& start,
                          const std::array<double, 3>& end,
                          int n_steps,
                          std::vector<std::array<double, 3>>& k_path,
                          std::vector<std::string>& labels,
                          const std::string& label_start,
                          std::vector<double>& k_dist,
                          std::vector<double>& k_ticks,
                          bool include_last = true) {
        double accum_dist = k_dist.empty() ? 0.0 : k_dist.back();

        for (int i = 0; i <= n_steps; ++i) {
            if (i == n_steps && !include_last) break;
            double t = static_cast<double>(i) / n_steps;
            std::array<double, 3> k = {
                (1 - t) * start[0] + t * end[0],
                (1 - t) * start[1] + t * end[1],
                0.0
            };
            k_path.push_back(k);

            // 누적 거리 계산
            if (!k_dist.empty()) {
                const auto& prev = k_path[k_path.size() - 2];
                double dk = std::sqrt(
                    std::pow(k[0] - prev[0], 2) +
                    std::pow(k[1] - prev[1], 2));
                accum_dist += dk;
            }
            k_dist.push_back(accum_dist);

            if (i == 0 && !label_start.empty()) {
                labels.push_back(label_start);
                k_ticks.push_back(accum_dist);
            } else {
                labels.push_back("");
            }
        }
    };

    // 5. Build path M → Γ → K → M
    std::vector<std::array<double, 3>> k_path;
    std::vector<std::string> k_labels;
    k_dist.clear();
    k_ticks.clear();

    int n_points = 50;
    interpolate(M, Gamma, n_points, k_path, k_labels, "M", k_dist, k_ticks);
    interpolate(Gamma, K, n_points, k_path, k_labels, "Γ", k_dist, k_ticks);
    interpolate(K, M, n_points, k_path, k_labels, "K", k_dist, k_ticks);

    return {k_path, k_labels};
}


// SOLVE EIGENVALUE PROBLEM
arma::mat solve_eigenvalue(arma::cx_mat &H,arma::mat &S, arma::vec &eigvals, arma::cx_mat &eigvecs)
{
    // diagonalize S to get its eigenvalues and eigenvectors
    arma::vec s_eigvals;
    arma::mat s_eigvecs;

    arma::eig_sym(s_eigvals, s_eigvecs, S); 

    // Regularize small eigenvalues to avoid instability
    /*double eps = 1e-6;
    for (size_t i = 0; i < s_eigvals.n_elem; ++i) {
        if (s_eigvals(i) < eps) s_eigvals(i) = eps;
    }*/

    // build X = S^(-1/2) = U diag(1/sqrt(sigma)) U^T
    // columns of X will be used to transform the basis so that the overlap matrix becomes the identity
    arma::vec inv_sqrt_vals = 1.0 / arma::sqrt(s_eigvals);
    arma::mat diag_inv_sqrt = arma::diagmat(inv_sqrt_vals);
    arma::mat X = s_eigvecs * diag_inv_sqrt * s_eigvecs.t();

    // enforce that each column's diagonal element is positive
    for (size_t j = 0; j < X.n_cols; j++) {
        if (X(j, j) < 0) {
            X.col(j) *= -1;
        }
    }

    // transform Hamiltonian: H' = X^T H X
    // we have orthonormal basis where overlap is the identity
    arma::cx_mat H_prime = X.t() * H * X;

    H_prime = 0.5 * (H_prime + H_prime.t());  // make it symmetric (Hermitian)

    // solve standard eigenvalue problem for H'
    arma::vec tmp_eigvals;
    arma::cx_mat tmp_eigvecs;
    arma::eig_sym(tmp_eigvals, tmp_eigvecs, H_prime);

    // transform eigenvectors back to original basis so that H C = S C e, following the generalized problem
    // C = X * V
    eigvals = tmp_eigvals;  
    eigvecs = X * tmp_eigvecs;

    // return X for printing
    return X;
}
