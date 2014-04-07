#ifndef MATRIX_H
#define MATRIX_H

#include <limits>
#include <sstream>
#include <cassert>
#include <iomanip>

#include <stdint.h>

/**
 * A matrix of bits.
 *
 * Bits are stored efficiently in blocks of an integer type (currently uint64_t).
 *
 * Example usage:
 *
 * Matrix M(4, 4);
 * M(1, 1) = true;
 * std::cout << M << std::endl;
 *
 * Output:
 *
 * 0 0 0 0
 * 0 1 0 0
 * 0 0 0 0
 * 0 0 0 0
 */
class Matrix {
public:
    // Integer type used for blocks of bits.
    typedef uint64_t Block;

    // Number of bits per block.
    const static uint32_t BitsPerBlock = std::numeric_limits<Block>::digits;

    // Helper class to reference a single bit.
    //
    // Used as return type for operator() (int, int) to allow e.g.
    // M(i, j) = true and if (!M(i, j)) { ... }.
    class Value {
    public:
        friend class Matrix;

        // Assign a bool to a bit.
        inline Value& operator=(bool rhs) {
            m_block = rhs ? m_block | m_mask : m_block & ~m_mask;
            return *this;
        }

        // Implicit conversion of bit to bool.
        inline operator bool() const {
            return (m_block & m_mask) != 0;
        }

        // Flip a bit.
        inline void flip() {
            m_block ^= m_mask;
        }

    private:
        Value();
        Value(const Matrix &matrix, uint32_t row, uint32_t col) :
            m_block(matrix.m_matrix[row][col / BitsPerBlock]),
            m_mask(Block(1) << col % BitsPerBlock) { }

        Block& m_block;  // Block this bit occurs in.
        Block m_mask;    // Mask to get bit in the block.
    };

    // Construct a new matrix with the given dimensions initialized with zeroes.
    Matrix(uint32_t rows, uint32_t cols) :
        m_blocksPerRow(cols / BitsPerBlock + 1),
        m_rows(rows),
        m_cols(cols)
    {
        m_matrix = new Block*[m_rows];
        for (uint32_t i = 0; i < m_rows; ++i) {
            m_matrix[i] = new Block[m_blocksPerRow];
            std::fill(m_matrix[i], m_matrix[i] + m_blocksPerRow, 0);
        }
    }

    // Construct a new matrix from another matrix.
    Matrix(const Matrix& other) {
        m_rows = other.rows();
        m_cols = other.cols();
        m_blocksPerRow = m_cols / BitsPerBlock + 1;

        m_matrix = new Block*[m_rows];
        for (uint32_t i = 0; i < m_rows; ++i) {
            m_matrix[i] = new Block[m_blocksPerRow];
            for (uint32_t j = 0; j < m_cols; ++j)
                (*this)(i, j) = (bool)other(i, j);
        }
    }

    // Load matrix data from a string.
    void load(const std::string &in) {
        std::stringstream ss(in);
        for (uint32_t i = 0; i < m_rows; ++i) {
            for (uint32_t j = 0; j < m_cols; ++j) {
                bool bit; ss >> bit;
                (*this)(i, j) = bit;
            }
        }
    }

    // Destructor.
    ~Matrix() {
        for (uint32_t i = 0; i < m_rows; ++i)
            delete[] m_matrix[i];
        delete[] m_matrix;
    }

    // Number of rows / cols.
    inline uint32_t rows() const { return m_rows; }
    inline uint32_t cols() const { return m_cols; }

    // Returns a reference to a single bit (const).
    inline Value operator() (uint32_t row, uint32_t col) const {
        return Value(*this, row, col);
    }

    // Returns a reference to a single bit (non-const).
    inline Value operator() (uint32_t row, uint32_t col) {
        return Value(*this, row, col);
    }

    // Adds row i to row j (mod 2), storing the result in row j.
    inline void addRows(uint32_t i, uint32_t j) {
        for (uint32_t k = 0; k < m_blocksPerRow; ++k)
            m_matrix[j][k] ^= m_matrix[i][k];
    }

    // Swaps row i with row j.
    inline void swapRows(uint32_t i, uint32_t j) {
        std::swap(m_matrix[i], m_matrix[j]);
    }

    // Clears row i, setting all elements to 0.
    inline void clearRow(uint32_t i) {
        std::fill(m_matrix[i], m_matrix[i] + m_blocksPerRow, 0);
    }

    // Reduces the matrix to row echelon form using Gaussian elimination.
    inline void reduce() {
        uint32_t i = 0, j = 0;
        while (i < rows() && j < cols()) {
            uint32_t maxi = i;
            // Find pivot element.
            for (uint32_t k = i + 1; k < rows(); ++k) {
                if ((*this)(k, j)) {
                    maxi = k;
                    break;
                }
            }
            if ((*this)(maxi, j)) {
                // Pivot.
                swapRows(i, maxi);
                for (uint32_t l = i + 1; l < rows(); ++l) {
                    if ((*this)(l, j)) {
                        addRows(i, l);
                    }
                }
                ++i;
            }
            ++j;
        }
    }

    /*
     * Perform back-substitution on a copy of the matrix and return a solution vector
     * x to Ax = b. Matrix is assumed to have been reduced and the system to be
     * underdetermined.
     */
    std::vector<uint32_t> solve() const {
        Matrix M(*this); // Work on a copy.

        std::vector<uint32_t> x(cols() - 1, 0);
        int32_t i = rows() - 1;
        while (i >= 0) {
            // Count the 1:s in the current row.
            int32_t count = 0;
            int32_t current = -1;
            for (uint32_t j = 0; j < cols() - 1; ++j) {
                count += M(i, j);
                current = M(i, j) ? j : current;
            }
            if (count == 0) {
                --i;
                continue; // Row was empty, proceed upwards.
            }

            // Introduce some randomness to avoid the trivial solution.
            uint32_t x_current = count > 1 ? rand() % 2 : M(i, cols() - 1);
            x[current] = x_current;

            for (int32_t k = 0; k <= i; ++k) {
                if (M(k, current)) {
                    if (x_current == 1)
                        M(k, cols() - 1).flip(); // Add to RHS.
                    M(k, current) = false;       // Remove from LHS.
                }
            }
            if (count == 1)
                --i; // Done with this row, proceed upwards.
        }

        return x;
    }

private:
    Matrix();     // Default constructor private.

    Block **m_matrix;        // Matrix blocks (row major order).
    uint32_t m_blocksPerRow; // Number of blocks per row.
    uint32_t m_rows;         // Number of rows (in bits).
    uint32_t m_cols;         // Number of columns (in bits).
};

// Output stream operator.
std::ostream& operator<<(std::ostream& os, const Matrix &matrix)
{
    for (uint32_t i = 0; i < matrix.rows(); ++i) {
        for (uint32_t j = 0; j < matrix.cols(); ++j) {
            os << std::left << std::setw(2) << matrix(i, j);
        }
        if (i < matrix.rows() - 1)
            os << std::endl;
    }
    return os;
}

#endif // MATRIX_H
