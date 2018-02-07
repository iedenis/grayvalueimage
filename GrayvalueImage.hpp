/* The program was created by Denis Ievlev
 * ID 321582918
 * Was checked with valgrind with flag --leak-check=full
 * Was compiled successfull by makefile
 */
#pragma once

#include "dice.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>
#include <cassert>
#include <functional>
#include <algorithm>
#include <limits>
#include <iterator>
#include "gtest/gtest.h" // for FRIEND_TEST, for testing private methods
template<typename F = unsigned char> class GrayvalueImage {

private:

	// A synomon for the type of the grayvalues
	typedef unsigned int grayvalue_t;
	static_assert(std::numeric_limits<grayvalue_t>::max()<=
			std::numeric_limits<size_t>::max(),
			"grayvalue_t maximum should be smaller than size_t maximum");

	// Number of rows
	size_t _R;
	// Number of columns
	size_t _C;
	// Maximum grayvalue
	grayvalue_t _MAX_G;
	// Pixels' grayvalues
	std::vector<grayvalue_t> _pixels;

	//Using for new indexes of the image
	std::vector<grayvalue_t> _cols;
	std::vector<grayvalue_t> _rows;

	// Reports an error that has occurred when reading
	// from the file filename and closes the stream and
	// exits the program.
	static void errorReport(const std::string& filename, std::fstream& fin) {
		std::cerr << "Image " << filename << " has a problem in its format"
				<< std::endl;
		fin.close();
		exit(1);
	}

	// Reads i from the stream and report if an error has occurred
	template<typename T>
	static void readi(const std::string& filename, std::fstream& fin, T& i) {
		fin >> i;
		if (!fin)
			errorReport(filename, fin);
	}

	void ctorsAdditionalWork() {

		assert(_R > 0 && _C > 0 && _MAX_G > 0);
		for (auto p : _pixels)
			assert(p <= _MAX_G);
		assert(_pixels.size() == _R * _C);

	}

public:

	/// Creates an image with R rows and C columns
	/// where the grayvalues are given in pixels and are
	/// between 0 and MAX_G
	GrayvalueImage(const size_t& R, const size_t& C, const grayvalue_t& MAX_G,
			const std::vector<grayvalue_t>& pixels) :
			_R(R), _C(C), _MAX_G(MAX_G), _pixels(pixels) {
		ctorsAdditionalWork();

	}

	/// Like previous constructor, only pixel values are uniform random number
	/// between 0 and MAX_G
	GrayvalueImage(const size_t& R, const size_t& C, const grayvalue_t& MAX_G) :
			GrayvalueImage(R, C, MAX_G, std::vector<grayvalue_t>(R * C)) {

		std::generate(std::begin(_pixels), std::end(_pixels),
				dice<grayvalue_t>(0, MAX_G));

	}

	/// Creates an image from a file
	/// Ignores mistakes in input
	GrayvalueImage(const std::string& filename) {

		std::fstream fin(filename.c_str());

		readi(filename, fin, _R);
		readi(filename, fin, _C);
		readi(filename, fin, _MAX_G);
		_rows.resize(_R);
		_cols.resize(_C);

		//vectors for indexes
		vector_initialization(_rows);
		vector_initialization(_cols);

		grayvalue_t tmp;
		while (fin >> tmp) {
			_pixels.push_back(tmp);
		}

		ctorsAdditionalWork();

	}

	/// Writes an image to a file
	void writeToFile(const std::string& filename) const {

		std::fstream fout(filename.c_str(),
				std::fstream::out | std::fstream::trunc);

		fout << _R << " " << _C << " " << _MAX_G << std::endl;

		for (size_t r = 0; r < R(); ++r) {
			for (size_t c = 0; c < C(); ++c) {
				fout << getPixel(r, c) << " ";
			}
			fout << getPixel(r, C() - 1) << std::endl;
		}

	}

	/// Returns the number of rows
	const size_t R() const {
		return _R;
	}

	/// Returns the number of columns
	const size_t C() const {
		return _C;
	}

	/// Retruns the maximum grayvalue
	const grayvalue_t MAX_G() const {
		return _MAX_G;
	}

	/// Returns the grayvalue of a pixel in row r and column c
	const grayvalue_t getPixel(const size_t& r, const size_t& c) const {
		assert(r < R() && c < C());
		//return _pixels[r * C() + c];
		return _pixels[_rows[r] * C() + _cols[c]];
	}

	/// Adds "pepper" noise.
	/// Each pixel has (percentage/100) probability to change to 0
	void addPepper(const size_t& percentage) {
		assert(percentage >= 0 && percentage <= 100);

		dice<size_t> d(1, 100);

		for (auto& p : _pixels) {
			if (d() < percentage)
				p = 0;
		}
	}

	
	void medianFilter(const size_t& radius, const unsigned int& method) {
		if (radius == 0)
			return;
		switch (method) {
		case 0:
			this->sorting(radius);
			break;
		case 1:
			this->nth_element(radius);
			break;
		case 2:
			this->histograms(radius);
			break;
		default:
			break;
		}
	}

	/*
	 * Using std::sort to find median in the vector by taking all
	 * neighborhood pixels and sorting them inside vector
	 */

	void sorting(const size_t& radius) {
		//Creating temp picture - newIm
		std::vector<std::vector<grayvalue_t> > newIm(this->_R,
				std::vector<grayvalue_t>(this->_C));

		for (int i = 0; i < this->_R; ++i) {
			for (int j = 0; j < this->_C; ++j) {
				int tmp_size = num_of_neighborhood_pixels(i, j, radius);
				std::vector<grayvalue_t> tmp(tmp_size);
				get_neighborhood_pixels(i, j, radius, tmp);
				std::sort(tmp.begin(), tmp.end());
				newIm[i][j] = tmp[(tmp.size() - 1) / 2];

			}
		}
		write_from_newIm(newIm);
	}

	/*
	 * Using std::nth_element to find median in the vector by taking all
	 * neighborhood pixels inside vector
	 */
	void nth_element(const size_t& radius) {
		std::vector<std::vector<grayvalue_t> > newIm(this->_R,
				std::vector<grayvalue_t>(this->_C));
		for (size_t i = 0; i < this->_R; ++i) {
			for (size_t j = 0; j < this->_C; ++j) {
				size_t tmp_size = num_of_neighborhood_pixels(i, j, radius);
				std::vector<grayvalue_t> tmp(tmp_size);
				get_neighborhood_pixels(i, j, radius, tmp);
				std::nth_element(tmp.begin(),
						tmp.begin() + (tmp.size() - 1) / 2, tmp.end());
				newIm[i][j] = tmp[(tmp.size() - 1) / 2];
			}
		}

		write_from_newIm(newIm);
	}

	void histograms(const size_t& rad) {
		std::vector<std::vector<grayvalue_t> > resultIm(this->_R,
				std::vector<grayvalue_t>(this->_C));

		size_t radius = rad;

		/* The second vector of histograms of all columns of the picture
		 * each column is (radius*2+1)
		 */
		std::vector<std::vector<grayvalue_t> > tmp_row(this->_C,
				std::vector<grayvalue_t>(this->_MAX_G + 1));
		//int tmp_row[_C][_MAX_G + 1]; Not in use because of error after running make
		// init
		for (int v = 0; v < _C; ++v) {
			for (int c = 0; c < _MAX_G + 1; ++c) {
				tmp_row[v][c] = 0;
			}
		}
		for (int v = 0; v < radius + 1; ++v) {
			for (int b = 0; b < _C; ++b) {
				int tmp = this->getPixel(v, b);
				tmp_row[b][tmp] += 1;
			}
		}

		// calculates histogram
		for (int i = 0; i < _R; ++i) {
			// Using arr 2 tmp_row
			for (int j = 0; j < _C; ++j) {
				if (isInImage(i - radius, j)) {
					tmp_row[j][this->getPixel(i - radius, j)] -= 1;
				}
			}
			for (int j = 0; j < _C; ++j) {
				if (isInImage(i + radius, j)) {
					tmp_row[j][this->getPixel(i + radius, j)] += 1;
				}
			}

			// tmp_cell - the first vector of histogram of all neighborhood pixels for the specific pixel

			/*grayvalue_t tmp_cell[_MAX_G + 1]; Not in use because of an error by make
			 * const size of array
			 */
			std::vector<grayvalue_t> tmp_cell(_MAX_G + 1);
			/*for (size_t k = 0; k < _MAX_G + 1; ++k) {
			 tmp_cell[k] = 0;
			 }*/
			for (size_t j = 0; j < radius + 1; ++j) {
				for (size_t k = 0; k < _MAX_G + 1; ++k) {
					tmp_cell[k] += tmp_row[j][k];
				}
			}
			for (int j = 0; j < _C; ++j) {
				// Uses tmp arr 1 - tmp_cell
				if (isInImage(0, j - radius)) {
					for (int k = 0; k < _MAX_G + 1; ++k) {
						tmp_cell[k] -= tmp_row[j - radius][k];
					}
				}
				if (isInImage(0, j + radius)) {
					for (int k = 0; k < _MAX_G + 1; ++k) {
						tmp_cell[k] += tmp_row[j + radius][k];
					}
				}
				//calculates a median from tmp_cell of histograms for pixel
				int acc = 0;
				int k;
				for (k = 0; acc < num_of_neighborhood_pixels(i, j, radius) / 2;
						++k) {
					acc += tmp_cell[k];
				}
				//((_MAX_G * (_MAX_G + 1)) / 2);

				resultIm[i][j] = k;
			}
		}

		write_from_newIm(resultIm);
	}

/// Return true if and only if
/// the dimensions of the images are equal,
/// the maximum grayvalue of the images are equal,
/// all the pixels values are equal
/// Time complexity: O(R()*C())
	const bool operator==(const GrayvalueImage& b) const {
		if (this->_pixels.size() != b._pixels.size()
				|| this->_MAX_G != b._MAX_G)
			return false;
		for (std::vector<grayvalue_t>::const_iterator i = this->_pixels.begin();
				i != this->_pixels.end(); i++)
			for (size_t j = 0; j != this->_C; j++) {
				if (this->_pixels != b._pixels)
					return false;
			}
		return true;
	}

/// Returns the opposite of operator==
/// Time complexity: O(R()*C())
	const bool operator!=(const GrayvalueImage& b) const {
		return !(*this == b);
	}

/// swaps between rows r1 and r2
/// Time complexity: O(1)
//swaps between rows using iter_swap from <algorithm>
	void swap_rows(const size_t& r1, const size_t& r2) {
		std::iter_swap(_rows.begin() + r1, _rows.begin() + r2);
	}

/// swaps between columns c1 and c2
/// Time complexity: O(1)
//swaps between columns using iter_swap from <algorithm>
	void swap_cols(const size_t& c1, const size_t& c2) {
		std::iter_swap(_cols.begin() + c1, _cols.begin() + c2);
	}

private:

	/*Returning number of neighborhood pixels around the specific pixel with radius
	 * without the specific pixel
	 */
	inline size_t num_of_neighborhood_pixels(const size_t row, const size_t col,
			int radius) {
		int counter = 0;
		for (int i = -radius; i <= radius; ++i) {
			for (int j = -radius; j <= radius; ++j) {
				counter += this->isInImage(row + i, col + j) ? 1 : 0;
			}
		}
		return counter - 1;
	}

	//Taking all neighborhood pixel around the specific pixel with radius
	void get_neighborhood_pixels(const size_t row, const size_t col, int radius,
			std::vector<grayvalue_t> & target) {
		int counter = 0;
		for (int i = -radius; i <= radius; ++i) {
			for (int j = -radius; j <= radius; ++j) {
				int cur_row = row + i;
				int cur_col = col + j;
				if (this->isInImage(cur_row, cur_col)
						&& !(cur_row == row && cur_col == col)) {
					target[counter] = this->getPixel(cur_row, cur_col);
					counter++;
				}
			}
		}

	}

	//Checking if the specific pixel if inside the picture or outside
	const inline bool isInImage(int x, int y) {
		return x >= 0 && y >= 0 && y < this->_C && x < this->_R;
	}

	//Printing vector
	void print_vector(std::vector<grayvalue_t>& vec) const {
		for (std::vector<grayvalue_t>::const_iterator i = vec.begin();
				i != vec.end(); ++i) {
			std::cout << *i << " ";
		}
	}

	//Changing the original picture taken from temporary newIm picture
	void write_from_newIm(std::vector<std::vector<grayvalue_t> > newIm) {
		for (int i = 0; i < this->_R; ++i) {
			for (int j = 0; j < this->_C; ++j) {
				this->_pixels[i * this->_C + j] = newIm[i][j];
			}
		}
	}

	//Initializing vector of indexes (rows and columns) using std::generate and lambda function inside
	void vector_initialization(std::vector<grayvalue_t>& v) {
		size_t i = 0;
		std::generate(v.begin(), v.end(), [&] {return i++;});
	}

};

