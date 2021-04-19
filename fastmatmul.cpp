#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
using std::ostream;
using std::cout;
using std::endl;
#include <vector>
using std::vector;

using Matrix = vector<vector<int>>;

class Benchmark {
public:
    Benchmark() {}
    void timestamp() {
        t1 = t2;
        t2 = std::chrono::high_resolution_clock::now();
    }
    double count() {
        std::chrono::duration<double, std::milli> n = t2 - t1;
        return n.count();
    }
private:
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;
};

int randint(int a, int b) {
    // random integer in range [a, b]
    return rand() % (b - a + 1) + a;
}

size_t nextPow2(size_t n) {
    return (size_t)pow(2.0, (ceil(log2((int)n))));
}

Matrix zeros(size_t N) {
    return Matrix(N, vector<int>(N));
}

// mxn matrix of random ints
Matrix rand_matrix(size_t N) {
    Matrix M = zeros(N);
    for (size_t r = 0; r < N; ++r) {
        for (size_t c = 0; c < N; ++c) {
            M[r][c] = randint(1,100);
        }
    }
    return M;
}

Matrix add(const Matrix& A, const Matrix& B) {
    size_t N = A.size();
    assert(N == A[0].size());
    assert(N == B.size());
    assert(N == B[0].size());
    Matrix C = zeros(N);
    for (size_t r = 0; r < N; ++r) {
        for (size_t c = 0; c < N; ++c) {
            C[r][c] = A[r][c] + B[r][c];
        }
    }
    return C;
}

Matrix add(const Matrix& A, const Matrix& B,
           size_t N,
           size_t A_r_begin, size_t A_c_begin,
           size_t B_r_begin, size_t B_c_begin) {
    Matrix C = zeros(N);
    for (size_t r = 0; r < N; ++r) {
        for (size_t c = 0; c < N; ++c) {
            C[r][c] = A[A_r_begin + r][A_c_begin + c] + B[B_r_begin + r][B_c_begin + c];
        }
    }
    return C;
}

Matrix sub(const Matrix& A, const Matrix& B) {
    size_t N = A.size();
    assert(N == A[0].size());
    assert(N == B.size());
    assert(N == B[0].size());
    Matrix C = zeros(N);
    for (size_t r = 0; r < N; ++r) {
        for (size_t c = 0; c < N; ++c) {
            C[r][c] = A[r][c] - B[r][c];
        }
    }
    return C;
}

Matrix sub(const Matrix& A, const Matrix& B,
           size_t N,
           size_t A_r_begin, size_t A_c_begin,
           size_t B_r_begin, size_t B_c_begin) {
    Matrix C = zeros(N);
    for (size_t r = 0; r < N; ++r) {
        for (size_t c = 0; c < N; ++c) {
            C[r][c] = A[A_r_begin + r][A_c_begin + c] - B[B_r_begin + r][B_c_begin + c];
        }
    }
    return C;
}

// returns submatrix of size sz starting at r, c
Matrix submatrix(const Matrix& M, size_t sz, size_t r_begin, size_t c_begin) {
    Matrix subm = zeros(sz);
    for (size_t dr = 0; dr < sz; ++dr) {
        for (size_t dc = 0; dc < sz; ++dc) {
            subm[dr][dc] = M[r_begin + dr][c_begin + dc];
        }
    }
    return subm;
}

// copies submatrix onto region of bigger matrix
void copy_submatrix(Matrix& M, const Matrix& subm, size_t r_begin, size_t c_begin) {
    size_t sz = subm.size();
    for (size_t dr = 0; dr < sz; ++dr) {
        for (size_t dc = 0; dc < sz; ++dc) {
            M[r_begin + dr][c_begin + dc] = subm[dr][dc];
        }
    }
}

ostream& operator<<(ostream& os, const Matrix& M) {
    size_t N = M.size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            os << M[i][j] << ' ';
        }
        os << endl;
    }
    return os;
}

bool operator==(const Matrix& A, const Matrix& B) {
    size_t N = A.size();
    assert(N == A[0].size());
    assert(N == B.size());
    assert(N == B[0].size());
    for (size_t r = 0; r < N; ++r) {
        for (size_t c = 0; c < N; ++c) {
            if (A[r][c] != B[r][c]) return false;
        }
    }
    return true;
}

// A*B = C
// row-major
// O(N^3)
Matrix matmul_naive(const Matrix& A, const Matrix& B) {
    size_t N = A.size();
    assert(N == A[0].size());
    assert(N == B.size());
    assert(N == B[0].size());
    Matrix C = zeros(N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            int sum = 0;
            for (size_t k = 0; k < N; ++k) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
    return C;
}

// only allows square matrices with dimensions as powers of 2
// amend threshold to have a naive-Strassen hybrid (use naive below N=threshold)
// O(N^2.8074)
Matrix matmul_strassen(const Matrix& A, const Matrix& B, size_t threshold=1) {
    size_t N = A.size();
    assert(N == A[0].size());
    assert(N == B.size());
    assert(N == B[0].size());

    Matrix C = zeros(N);

    if (N <= threshold) { // base case
        return matmul_naive(A, B);
    }

    Matrix A11 = submatrix(A, N/2, 0, 0),
           A12 = submatrix(A, N/2, 0, N/2),
           A21 = submatrix(A, N/2, N/2, 0),
           A22 = submatrix(A, N/2, N/2, N/2),
           B11 = submatrix(B, N/2, 0, 0),
           B12 = submatrix(B, N/2, 0, N/2),
           B21 = submatrix(B, N/2, N/2, 0),
           B22 = submatrix(B, N/2, N/2, N/2);

    Matrix M1 = matmul_strassen(add(A11, A22), add(B11, B22), threshold),
           M2 = matmul_strassen(add(A21, A22), B11, threshold),
           M3 = matmul_strassen(A11, sub(B12, B22), threshold),
           M4 = matmul_strassen(A22, sub(B21, B11), threshold),
           M5 = matmul_strassen(add(A11, A12), B22, threshold),
           M6 = matmul_strassen(sub(A21, A11), add(B11, B12), threshold),
           M7 = matmul_strassen(sub(A12, A22), add(B21, B22), threshold);

    Matrix C11 = add(sub(add(M1, M4), M5), M7),
           C12 = add(M3, M5),
           C21 = add(M2, M4),
           C22 = add(sub(add(M1, M3), M2), M6);

    copy_submatrix(C, C11, 0, 0);
    copy_submatrix(C, C12, 0, N/2);
    copy_submatrix(C, C21, N/2, 0);
    copy_submatrix(C, C22, N/2, N/2);

    return C;
}

double vector_avg(const vector<double>& v) {
    double sum = 0;
    for (double n : v) {
        sum += n;
    }
    return sum / (double)v.size();
}

double vector_stdev(const vector<double>& v) {
    double avg = vector_avg(v);
    double sum_var = 0;
    for (double n : v) {
        sum_var += (n-avg)*(n-avg);
    }
    return sqrt(sum_var / (double)v.size());
}

// constructs random matrices of size N, and compares matmul performances
void benchmark_comp(size_t N, int numTrials=1) {
    Benchmark bench;
    vector<double> naive_perf;
    vector<double> strassen_perf;
    vector<double> hybrid_perf;
    for (int i = 0; i < numTrials; ++i) {
        Matrix A = rand_matrix(N);
        Matrix B = rand_matrix(N);

        bench.timestamp();
        Matrix C_naive = matmul_naive(A, B);
        bench.timestamp();
        cout << "naive: " << bench.count() << endl;
        naive_perf.push_back(bench.count());

        bench.timestamp();
        Matrix C_strassen = matmul_strassen(A, B, 16);
        bench.timestamp();
        cout << "strassen: " << bench.count() << endl;
        strassen_perf.push_back(bench.count());

        bench.timestamp();
        Matrix C_hybrid = matmul_strassen(A, B, 512);
        bench.timestamp();
        cout << "hybrid: " << bench.count() << endl;
        hybrid_perf.push_back(bench.count());

        assert(C_naive == C_strassen);
        assert(C_naive == C_hybrid);
    }
    cout << "Benchmarking for N=" << N << endl;
    cout << "Num trials: " << numTrials << endl;
    cout << "naive_perf: avg=" << vector_avg(naive_perf) << ", stdev=" << vector_stdev(naive_perf) << endl;
    cout << "strassen_perf: avg=" << vector_avg(strassen_perf) << ", stdev=" << vector_stdev(strassen_perf) << endl;
    cout << "hybrid_perf: avg=" << vector_avg(hybrid_perf) << ", stdev=" << vector_stdev(hybrid_perf) << endl;
    cout << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: fastmatmul N [numTrials]" << endl;
        cout << "Conducts benchmarking for square matrix multiplications\n" <<
                "of size N over numTrials trials" << endl;
        return 0;
    }
    int N_int = atoi(argv[1]);
    size_t N = (size_t)N_int;
    if (nextPow2(N) != N) {
        cout << "N must be a power of 2" << endl;
        return 0;
    }
    if (argc == 2) {
        benchmark_comp(N);
    } else {
        int numTrials = atoi(argv[2]);
        benchmark_comp(N, numTrials);
    }
    return 0;
}
