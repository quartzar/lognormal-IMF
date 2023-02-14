#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

#include "gnuplot-iostream.h"

int main() {
    int N = 1000;
    double mu = 0.1;        // Mean in base_10 log space       [Chabrier 2002]
    double sigma = 0.627;   // Std. Deviations in log 10 space [Chabrier 2002]
    
    double ln_mu = mu * std::log(10); // [base_e conversion]
    double ln_sigma = sigma * std::log(10); // [base_e conversion]
    
    // Generate the lognormal distribution random masses in natural log space
    std::vector<double> ln_mass(N);
    std::random_device rd;
    std::mt19937_64 rng(rd()); // Mersenne Twister 19937 generator with 64-bit output
    std::lognormal_distribution<double> dist(ln_mu, ln_sigma);
    for (int i = 0; i < N; i++) {
        ln_mass[i] = dist(rng);
    }
    
    // Convert the natural logarithm values back to base-10 logarithm space
    std::vector<double> log_mass(N);
    for (int i = 0; i < N; i++) {
        log_mass[i] = std::log10(std::exp(ln_mass[i]));
    }
    
    Gnuplot gp;
    gp << "set terminal pngcairo font 'DejaVu Sans Mono,12' enhanced size 800,600\n";
    // gp << "set terminal x11 font 'DejaVu Sans Mono.12'\n";
    gp << "set output 'logmass_histogram.png'\n";
    gp << "set boxwidth 0.9 relative\n";
    gp << "set style fill solid 0.5\n";
    gp << "set xlabel 'Log mass (solar masses)'\n";
    gp << "set ylabel 'Number of stars'\n";
    gp << "binwidth=0.1\n";
    gp << "bin(x,width)=width*floor(x/width)+binwidth/2.0\n";
    gp << "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes title 'Lognormal initial mass function'\n";
    for (int i = 0; i < N; i++) {
        gp << std::log10(log_mass[i]) << "\n";
    }
    gp << "e\n";
    
    return 0;
}