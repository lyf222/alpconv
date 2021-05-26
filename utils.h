#include <vector>
#include <string>
#include <iostream>
#include <fstream>

void printv(const std::vector<double>& v);

std::vector<double> linspace(double start, double stop, int num);

std::vector<double> arange(double start, double stop, double step);

void split(const std::string& s, std::vector<std::string>& tokens, const std::string& delimiters = " ");

double quad(std::vector<double> &xx, std::vector<double> &yy);

double interp1d(std::vector<double> &xvec, std::vector<double> &yvec, double xx);

std::vector<std::vector<double>> readtxt(const char* filename, int ncols, int nSkip);
