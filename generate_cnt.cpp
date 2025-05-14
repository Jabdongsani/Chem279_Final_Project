#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <array>
#include <string>


// we can test (6, 0) semimetallic CNT, 5, 0 -> metallic 9,0  5.5 metalic 

const double PI = 3.14159265358979323846;
const double angstrom_to_bohr = 1.8897259886;

struct Atom {
    std::string element;
    double x, y, z;
};

// Graphene lattice vectors in 2D
std::array<double, 3> a1(double a_cc) {
    return {a_cc * 1.5,  a_cc * std::sqrt(3) / 2.0, 0.0};
}
std::array<double, 3> a2(double a_cc) {
    return {a_cc * 1.5, -a_cc * std::sqrt(3) / 2.0, 0.0};
}

// Calculate chiral vector C_h = n a1 + m a2
std::array<double, 3> chiral_vector(int n, int m, double a_cc) {
    auto a1_v = a1(a_cc);
    auto a2_v = a2(a_cc);
    return {
        n * a1_v[0] + m * a2_v[0],
        n * a1_v[1] + m * a2_v[1],
        0.0
    };
}

// Length of chiral vector
double chiral_length(int n, int m, double a_cc) {
    auto ch = chiral_vector(n, m, a_cc);
    return std::sqrt(ch[0]*ch[0] + ch[1]*ch[1]);
}

// Map (x, y) on graphene to CNT surface (x, y, z)
Atom wrap_atom_to_cnt(double x, double y, double z_graphene, double radius, double chiral_len) {
    double theta = 2.0 * PI * x / chiral_len;
    double new_x = radius * std::cos(theta);
    double new_y = radius * std::sin(theta);
    return {"C", new_x, new_y, z_graphene};
}

// Generate 2-atom graphene basis in rectangular unit cell
std::vector<Atom> graphene_unit_cell(double a_cc) {
    std::vector<Atom> atoms;
    atoms.push_back({"C", 0.0, 0.0, 0.0});
    atoms.push_back({"C", a_cc, 0.0, 0.0});
    return atoms;
}

// Main CNT generation routine
void generate_cnt_structure(int n, int m, double a_cc, int num_z_repeat, const std::string& output_file) {
    double ch_len = chiral_length(n, m, a_cc);
    double radius = ch_len / (2.0 * PI);
    double unit_height = a_cc * std::sqrt(3); // approximate unit cell height

    std::vector<Atom> atoms;
    auto base_atoms = graphene_unit_cell(a_cc);

    for (int iz = 0; iz < num_z_repeat; ++iz) {
        for (const auto& atom : base_atoms) {
            double z = iz * unit_height;
            atoms.push_back(wrap_atom_to_cnt(atom.x, atom.y, z, radius, ch_len));
        }
    }

    std::ofstream fout(output_file);
    fout << atoms.size() << " 1\n";
    fout << ch_len * angstrom_to_bohr << " 0.0 0.0\n";
    fout << "0.0 0.0 0.0\n0.0 0.0 0.0\n";

    for (const auto& atom : atoms) {
        fout << "6 "
             << atom.x * angstrom_to_bohr << " "
             << atom.y * angstrom_to_bohr << " "
             << atom.z * angstrom_to_bohr << "\n";
    }
    fout.close();
    std::cout << "CNT structure written to " << output_file << std::endl;
}

int main() {
    int n = 5, m = 5;              // chirality (n,m)
    double a_cc = 1.42;            // C–C bond length in angstroms
    int num_repeat = 4;           // number of unit cells along z
    std::string output = "cnt_input_5_5.dat";

    generate_cnt_structure(n, m, a_cc, num_repeat, output);
    return 0;
}
