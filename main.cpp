#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <complex>
#include <thread>
#include "bitmap/bitmap_image.hpp"
using namespace std;


typedef vector<vector<complex<double>>> matrix;
typedef vector<complex<double>> ComplVector;

const double pi = 3.141592653589793238463;

complex<double> imagePart(0, 1);

int WIDTH = 0;
int HEIGHT = 0;


complex<double> adamarProduct(complex<double> a, complex<double> b) {
	complex<double> numerator = a * conj(b);
	complex<double> denumerator = complex<double>(abs(numerator));
	complex<double> res = complex<double>(numerator/denumerator);
	return res;
}

matrix readImage(string filename) {
	ifstream fin;
	fin.open(filename);
	fin >> HEIGHT >> WIDTH;
	matrix image;
	for (int x = 0; x < ::HEIGHT; x++) {
		image.emplace_back();
		for (int y = 0; y < ::WIDTH; y++) {
			int val;
			fin >> val;
			image[x].emplace_back(val, 0);
		}
	}
	fin.close();
	return image;
}

void writeImage(string filename, matrix image) {
	bitmap_image res_image(image[0].size(), image.size());
	for (int x = 0; x < image.size(); x++) {
		for (int y = 0; y < image[0].size(); y++) {
			rgb_t c;
			c.blue = char(image[x][y].real());
			c.green = char(image[x][y].real());
			c.red = char(image[x][y].real());
			res_image.set_pixel(y, x, c);
		}
	}
	res_image.save_image(filename);
}

matrix product(matrix a, matrix b) {
	matrix c(a.size());
	for (int i = 0; i < a.size(); i++) {
		c[i].resize(b.size());
		for (int k = 0; k < a[0].size(); k++) {
			for (int j = 0; j < b[0].size(); j++) {
				c[i][j] += complex<double>(a[i][k] * b[k][j]);
			}
			// cout << i << ' ' << k << "\n";
		}
	}
	return c;
}

matrix productCOMP(matrix a, matrix b) {
	matrix c(a.size(), ComplVector(b[0].size(), complex<double>(0, 0)));
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < a[0].size(); j++) {
			complex<double> val(0,0);
			for (int k = 0; k < b[0].size(); k++) {
				val += complex<double>(a[i][k] * b[k][j]);
			}
			c[i][j] = val;
		}
	}
	return c;
}

matrix createF(int N) {
	double PI = pi*(2);
	matrix F(N, ComplVector(N, complex<double>(0, 0)));
	for(int k = 0; k < N; k++) {
		for(int n = 0; n < N; n++) {
			F[k][n] = exp(complex<double>(0.0 , PI * k * n / N));
		}
	}
	return F;
}

matrix createHCF(int N) {
	double PI = pi*(2);
	matrix FH(N, ComplVector(N, complex<double>(0, 0)));
	for (int k = 0; k < N; k++) {
		for (int n = 0; n < N; n++) {
			FH[k][n] = conj(exp(complex<double>(0.0, PI * k * n / N)));
		}
	}
	return FH;
}

matrix ermitTransp(matrix m) {
	matrix res(m[0].size());
	for(int x = 0; x < m.size(); x++) {
		for(int y = 0; y < m[0].size(); y++) {
			res[y].emplace_back(conj(m[x][y]));
		}
	}
	return res;
}

//Digital Furie Transform
matrix DFT(matrix image) {
	matrix spectr(image.size(), ComplVector(image[0].size(), complex<double>(0, 0)));
	
	int count = 0;
	spectr = product(image, createHCF(image[0].size()));
	spectr = ermitTransp(spectr);
	spectr = productCOMP(spectr, createHCF(image.size()));
	return spectr;
}

matrix crossPowerSpectr(matrix spectr1, matrix spectr2) {
	matrix crossPowSpectr;

	int count = 0;
	cout << spectr1.size() << ' ' << spectr1[0].size() << '\n';
	for(int x = 0; x < spectr1.size(); x++) {
		crossPowSpectr.emplace_back();
		for(int y = 0; y < spectr1[0].size(); y++) {
			crossPowSpectr[x].emplace_back(adamarProduct(spectr1[x][y], spectr2[x][y]));
			count++;
		}
	}
	return crossPowSpectr;
}


//Inverse Furie Transform
matrix IFT(matrix CrossPowSpectr) {
	matrix image(CrossPowSpectr.size(), ComplVector(CrossPowSpectr[0].size(), complex<double>(0, 0)));
	int count = 0;
	image = product(CrossPowSpectr, createF(CrossPowSpectr[0].size()));
	image = ermitTransp(image);
	image = productCOMP(image, createF(CrossPowSpectr.size()));
	return image;
}

matrix peakByCorrelation(matrix correlation) {

	double max = 0;
	double min = correlation[0][0].real();
	int max_X = 0;
	int max_Y = 0;

	matrix res;

	for(int x = 0; x < correlation.size(); x++) {
		for(int y = 0; y < correlation[0].size(); y++) {
			if(correlation[x][y].real() > max) {
				max = correlation[x][y].real();
				max_X = x;
				max_Y = y;
			}
			if(correlation[x][y].real() < min) {
				min = correlation[x][y].real();
			}
		}
	}
	for(int x = 0; x < correlation.size(); x++) {
		res.emplace_back();
		for(int y = 0; y < correlation[0].size(); y++) {
			res[x].emplace_back();
			res[x][y].real((correlation[x][y].real() - min) / (max-min)*255);
		}
	}
	cout << "x: " << max_X << " | y: " << max_Y << '\n';
	return res;
}

int main() {
	ofstream fout;
	fout << fixed << setprecision(22);

	system("python3 ./readImag.py tank_first.bmp image1.txt");
	system("python3 ./readImag.py tank_second.bmp image2.txt");

	matrix image1 = readImage("image1.txt");
	matrix image2 = readImage("image2.txt");

	cout << WIDTH << " * " << HEIGHT << " | " << WIDTH * HEIGHT << '\n';

	cout << "DFT1" << '\n';
	matrix spectr1 = DFT(image1);
	cout << "DFT2" << '\n';
	matrix spectr2 = DFT(image2);
	
	cout << "POW" << '\n';
	matrix crossPowerSpectrum = crossPowerSpectr(spectr1, spectr2);

	cout << "IFT" << '\n';
	matrix phaseCorrelation = IFT(crossPowerSpectrum);

	cout << "COR" << '\n';
	matrix resCorrelation = peakByCorrelation(phaseCorrelation);

	writeImage("tank_res.bmp", resCorrelation);
	system("xdg-open tank_res.bmp");
}