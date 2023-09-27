#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>
#include <eigen3/Eigen/Dense>
#include <iomanip>
using namespace std;
typedef std::complex<double> dcomp;

void read_object(int steps, int num_orbitals, std::vector<std::vector<Eigen::MatrixXcd>> &object, std::string filename, std::vector<double> &energy) {
	dcomp j1 = -1;
	j1 = sqrt(j1);
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
            fstream my_file;
            std::ostringstream oss;
            oss << filename << "_" << i << "_" << j << ".dat";
            std::string var = oss.str(), line;
            std::cout << "reading file " << var << std::endl;
            //double delta_energy;
            //double num_old = 0;
            int line_count = -1;
            my_file.open(var, ios::in);
            if (my_file.is_open()) {
                while (std::getline(my_file, line)) {
                    line_count++;
                    std::istringstream iss(line);
                    double num1, num2, num3;

                    if (iss >> num1 >> num2 >> num3) {
                        object.at(line_count).at(0)(i, j) = num2 + num3 * j1;
                        object.at(line_count).at(1)(i, j) = num2 + num3 * j1;
                        energy.at(line_count) = num1;
                    }
                }
                my_file.close();
            } else {
                std::cout << var << " File cannot be opened for reading." << std::endl;
            }
        }
    }
}

void get_hybridisation(int steps, int num_orbitals, std::vector<std::vector<Eigen::MatrixXcd>> &gf_retarded,
    std::vector<std::vector<Eigen::MatrixXcd>> &gf_lesser, std::vector<std::vector<Eigen::MatrixXcd>> &se_lesser, 
    std::vector<std::vector<Eigen::MatrixXcd>> &hybridisation_lesser) {
    for (int r = 0; r < steps; r ++) {
        for (int spin = 0; spin < 2; spin++) {
            Eigen::MatrixXcd gf_advanced(num_orbitals, num_orbitals);
            Eigen::MatrixXcd gf_advanced_inverse(num_orbitals, num_orbitals);
            Eigen::MatrixXcd gf_retarded_inverse(num_orbitals, num_orbitals);
            gf_advanced = gf_retarded.at(r).at(spin).adjoint();
            gf_advanced_inverse = gf_advanced.inverse();
            gf_retarded_inverse = gf_retarded.at(r).at(spin).inverse();
            hybridisation_lesser.at(r).at(spin) = gf_retarded_inverse * gf_lesser.at(r).at(spin) * gf_advanced_inverse - se_lesser.at(r).at(spin);

       }
    }
}

void print_to_file(int steps, int num_orbitals, std::vector<std::vector<Eigen::MatrixXcd>> &object, std::string filename, std::vector<double> &energy) {
    for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
			std::ostringstream ossgf;
			ossgf << filename << "_" << i << "_" << j << ".dat";
			std::string var = ossgf.str();
			std::ofstream gf_local_file;
			gf_local_file.open(var);
			for (int r = 0; r < steps; r++) {
				gf_local_file << std::setprecision(15) << energy.at(r) << "  " << object.at(r).at(0)(i, j).real() << "   " << object.at(r).at(0)(i, j).imag() << "\n";
			}
			gf_local_file.close();	
        }
    }   
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " number_of_energy_points number_of_orbitals " <<   std::endl;
        return 1;
    }
    int steps = stoi(argv[1]);
    int num_orbitals = stoi(argv[2]);

    std::vector<double> energy(steps, 0);

    //in the format gfmat(neg, spin, nc, nc)
    std::vector<std::vector<Eigen::MatrixXcd>> gf_retarded(steps, std::vector<Eigen::MatrixXcd>(2, Eigen::MatrixXcd::Zero(num_orbitals, num_orbitals)));
    std::vector<std::vector<Eigen::MatrixXcd>> gf_lesser(steps, std::vector<Eigen::MatrixXcd>(2, Eigen::MatrixXcd::Zero(num_orbitals, num_orbitals)));
    std::vector<std::vector<Eigen::MatrixXcd>> se_lesser(steps, std::vector<Eigen::MatrixXcd>(2, Eigen::MatrixXcd::Zero(num_orbitals, num_orbitals)));
    std::vector<std::vector<Eigen::MatrixXcd>> hybridisation_lesser(steps, std::vector<Eigen::MatrixXcd>(2, Eigen::MatrixXcd::Zero(num_orbitals, num_orbitals)));
    
    read_object(steps, num_orbitals, gf_retarded, "gf_retarded", energy);
    read_object(steps, num_orbitals, gf_lesser, "gf_lesser", energy);
    read_object(steps, num_orbitals, se_lesser, "se_lesser", energy);
    get_hybridisation(steps, num_orbitals, gf_retarded, gf_lesser, se_lesser, hybridisation_lesser);
    print_to_file(steps, num_orbitals, hybridisation_lesser, "hybridisation_lesser", energy);
    return 1; 
}
