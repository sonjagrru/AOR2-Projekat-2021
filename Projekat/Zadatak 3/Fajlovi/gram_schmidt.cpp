/**
 * @file
 * @brief [Gram Schmidt Orthogonalisation
 * Process](https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process)
 *
 * @details
 * Takes the input of Linearly Independent Vectors,
 * returns vectors orthogonal to each other.
 *
 * ### Algorithm
 * Take the first vector of given LI vectors as first vector of Orthogonal
 * vectors. Take projection of second input vector on the first vector of
 * Orthogonal vector and subtract it from the 2nd LI vector. Take projection of
 * third vector on the second vector of Othogonal vectors and subtract it from
 * the 3rd LI vector. Keep repeating the above process until all the vectors in
 * the given input array are exhausted.
 *
 * For Example:
 * In R2,
 * Input LI Vectors={(3,1),(2,2)}
 * then Orthogonal Vectors= {(3, 1),(-0.4, 1.2)}
 *
 *  Have defined maximum dimension of vectors to be 10 and number of vectors
 *  taken is 20.
 *  Please do not give linearly dependent vectors
 *
 *
 * @author [Akanksha Gupta](https://github.com/Akanksha-Gupta920)
 */

#include <array>     /// for std::array
#include <cassert>   /// for assert
#include <cmath>     /// for fabs
#include <iostream>  /// for io operations
#include <immintrin.h>

#include "sys/timeb.h"
#include "emmintrin.h"
#include "stdio.h"
#include "math.h"
#include "_Timer.h"

using namespace std;

 /**
  * @namespace linear_algebra
  * @brief Linear Algebra algorithms
  */
namespace linear_algebra {
    /**
     * @namespace gram_schmidt
     * @brief Functions for [Gram Schmidt Orthogonalisation
     * Process](https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process)
     */
    namespace gram_schmidt {
        /**
         * Dot product function.
         * Takes 2 vectors along with their dimension as input and returns the dot
         * product.
         * @param x vector 1
         * @param y vector 2
         * @param c dimension of the vectors
         *
         * @returns sum
         */
        double dot_product(const std::array<double, 30>& x,
            const std::array<double, 30>& y, const int& c) {
            /*ORIGINAL*/
            double sum = 0;
            for (int i = 0; i < c; i++) {
                sum += x[i] * y[i];
            }
            return sum;
        }

        double dot_productO(const std::array<double, 30>& x,
            const std::array<double, 30>& y, const int& c) {
            /*OPTIMIZOVANO*/
           // if (c > 9) {
                double sum = 0;
                register int i = 0;
               // for (; i < 9; i++)
                //{
                    _mm_prefetch((char*)(&x[i]), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i]+1), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i]+2), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i]+3), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i]+4), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 5), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 6), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 7), _MM_HINT_T1);
                    _mm_prefetch((char*)(&y[i]), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i]+1), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i]+2), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i]+3), _MM_HINT_T1);
                    _mm_prefetch((char*)(&y[i]+4), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 5), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 6), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 7), _MM_HINT_T1);
                //}

                for (i = 1; i < c; i++) {
                    _mm_prefetch((char*)(&x[i]), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 1), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 2), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 3), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 4), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 5), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 6), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 7), _MM_HINT_T1);
                    _mm_prefetch((char*)(&y[i]), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 1), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 2), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 3), _MM_HINT_T1);
                    _mm_prefetch((char*)(&y[i] + 4), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 5), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 6), _MM_HINT_T1);
                    _mm_prefetch((char*)(&x[i] + 7), _MM_HINT_T1);
                    sum += x[i - 1] * y[i - 1];
                }
                sum += x[c - 1] * y[c - 1];
               /* for ( i = c-9; i < c; i++)
                {
                    sum += x[i] * y[i];
                }*/
                return sum;
            //}
            //else return dot_product(x, y, c);
        }

        /*OPTIMIZOVANO VEKTORSKI*/
        double dot_productO1(const std::array<double, 30>& x,
            const std::array<double, 30>& y, const int& c) {
            /*ORIGINAL*/
            double sum = 0;
            register int i = 0;
            int ostatak = c % 4;
            int kolicnik = c / 4;
            if (kolicnik != 0) {
                __m256d v1;
                __m256d v2;
                __m256d suma = _mm256_set1_pd(0);
                for (; i < kolicnik; i++) {
                    v1 = _mm256_load_pd(&x[4*i]);
                    v2 = _mm256_load_pd(&y[4*i]);
                    suma = _mm256_fmadd_pd(v1, v2, suma);
                }
                for (i = 0; i < 4; i++)
                {
                    sum += suma.m256d_f64[i];
                }
                v1 = _mm256_load_pd(&x[4 * kolicnik]);
                v2 = _mm256_load_pd(&y[4 * kolicnik]);
                suma = _mm256_mul_pd(v1, v2);
                for (i = 0;i < ostatak; i++)
                {
                    sum += suma.m256d_f64[i];
                }
                return sum;
            }
            else return dot_product(x, y, c);
        }

        /**
         * Projection Function
         * Takes input of 2 vectors along with their dimension and evaluates their
         * projection in temp
         *
         * @param x Vector 1
         * @param y Vector 2
         * @param c dimension of each vector
         *
         * @returns factor
         */
        double projection(const std::array<double, 30>& x,
            const std::array<double, 30>& y, const int& c) {
            double dot =
                dot_product(x, y, c);  /// The dot product of two vectors is taken
            double anorm =
                dot_product(y, y, c);  /// The norm of the second vector is taken.
            double factor =
                dot /
                anorm;  /// multiply that factor with every element in a 3rd vector,
                        /// whose initial values are same as the 2nd vector.
            return factor;
        }

        /*OPTIMIZOVANO*/
        double projectionO(const std::array<double, 30>& x,
            const std::array<double, 30>& y, const int& c) {
            double dot =
                dot_productO(x, y, c);  /// The dot product of two vectors is taken
            double anorm =
                dot_productO(y, y, c);  /// The norm of the second vector is taken.
            double factor =
                dot /
                anorm;  /// multiply that factor with every element in a 3rd vector,
                        /// whose initial values are same as the 2nd vector.
            return factor;
        }

        /*OPTIMIZOVANO VEKTORSKI*/
        double projectionO1(const std::array<double, 30>& x,
            const std::array<double, 30>& y, const int& c) {
            double dot =
                dot_productO1(x, y, c);  /// The dot product of two vectors is taken
            double anorm =
                dot_productO1(y, y, c);  /// The norm of the second vector is taken.
            double factor =
                dot /
                anorm;  /// multiply that factor with every element in a 3rd vector,
                        /// whose initial values are same as the 2nd vector.
            return factor;
        }

        /**
         * Function to print the orthogonalised vector
         *
         * @param r number of vectors
         * @param c dimenaion of vectors
         * @param B stores orthogonalised vectors
         *
         * @returns void
         */
        void display(const int& r, const int& c,
            const std::array<std::array<double, 30>, 30>& B) {
            /*ORIGINAL*/
            for (int i = 0; i < r; ++i) {
                std::cout << "Vector " << i + 1 << ": ";
                for (int j = 0; j < c; ++j) {
                    std::cout << B[i][j] << " ";
                }
                std::cout << '\n';
            }
        }

        void displayO(const int& r, const int& c,
            const std::array<std::array<double, 30>, 30>& B) {
            /*OPTIMIZOVANO*/
           // if (c > 8) {
                for (int i = 0; i < r; ++i) {
                    std::cout << "Vector " << i + 1 << ": ";

                    _mm_prefetch((char*)(&B[i][0]), _MM_HINT_T1);
                    _mm_prefetch((char*)(&B[i][0]+1), _MM_HINT_T1);
                    _mm_prefetch((char*)(&B[i][0]+2), _MM_HINT_T1);
                    _mm_prefetch((char*)(&B[i][0]+3), _MM_HINT_T1);
                    _mm_prefetch((char*)(&B[i][0]+4), _MM_HINT_T1);
                    _mm_prefetch((char*)(&B[i][0]+5), _MM_HINT_T1);
                    _mm_prefetch((char*)(&B[i][0]+6), _MM_HINT_T1);
                    _mm_prefetch((char*)(&B[i][0]+7), _MM_HINT_T1);

                    for (int j = 1; j < c; ++j) {
                        _mm_prefetch((char*)(&B[i][0]), _MM_HINT_T1);
                        _mm_prefetch((char*)(&B[i][0] + 1), _MM_HINT_T1);
                        _mm_prefetch((char*)(&B[i][0] + 2), _MM_HINT_T1);
                        _mm_prefetch((char*)(&B[i][0] + 3), _MM_HINT_T1);
                        _mm_prefetch((char*)(&B[i][0] + 4), _MM_HINT_T1);
                        _mm_prefetch((char*)(&B[i][0] + 5), _MM_HINT_T1);
                        _mm_prefetch((char*)(&B[i][0] + 6), _MM_HINT_T1);
                        _mm_prefetch((char*)(&B[i][0] + 7), _MM_HINT_T1);
                        std::cout << B[i][j - 1] << " ";
                    }
                   /* std::cout << B[i][c - 8] << " ";
                    std::cout << B[i][c - 7] << " ";
                    std::cout << B[i][c - 6] << " ";
                    std::cout << B[i][c - 5] << " ";
                    std::cout << B[i][c - 4] << " ";
                    std::cout << B[i][c - 3] << " ";
                    std::cout << B[i][c - 2] << " ";*/
                    std::cout << B[i][c - 1] << " ";
                    std::cout << '\n';
                }
           // }
            //else display(r, c, B);
        }
        
        

        /**
         * Function for the process of Gram Schimdt Process
         * @param r number of vectors
         * @param c dimension of vectors
         * @param A stores input of given LI vectors
         * @param B stores orthogonalised vectors
         *
         * @returns void
         */
        void gram_schmidt(int r, const int& c,
            const std::array<std::array<double, 30>, 30>& A,
            std::array<std::array<double, 30>, 30> B) {
            if (c < r) {  /// we check whether appropriate dimensions are given or not.
                std::cout << "Dimension of vector is less than number of vector, hence "
                    "\n first "
                    << c << " vectors are orthogonalised\n";
                r = c;
            }

            int k = 1;

            while (k <= r) {
                if (k == 1) {
                    for (int j = 0; j < c; j++)
                        B[0][j] = A[0][j];  /// First vector is copied as it is.
                }

                else {
                    std::array<double, 30>
                        all_projection{};  /// array to store projections
                    for (int i = 0; i < c; ++i) {
                        all_projection[i] = 0;  /// First initialised to zero
                    }

                    int l = 1;
                    while (l < k) {
                        std::array<double, 30>
                            temp{};           /// to store previous projected array
                        double factor = NAN;  /// to store the factor by which the
                                              /// previous array will change
                        factor = projection(A[k - 1], B[l - 1], c);

                        /*ORIGINAL*/
                        for (int i = 0; i < c; ++i) {
                            temp[i] = B[l - 1][i] * factor;  /// projected array created
                        }
                        for (int j = 0; j < c; ++j) {
                            all_projection[j] =
                                all_projection[j] +
                                temp[j];  /// we take the projection with all the
                                          /// previous vector and add them.
                        }
                        l++;
                    }
                    for (int i = 0; i < c; ++i) {
                        B[k - 1][i] =
                            A[k - 1][i] -
                            all_projection[i];  /// subtract total projection vector
                                                /// from the input vector
                    }
                }
                k++;
            }
            display(r, c, B);  // for displaying orthogoanlised vectors
        }

        void gram_schmidtO1(int r, const int& c,
            const std::array<std::array<double, 30>, 30>& A,
            std::array<std::array<double, 30>, 30> B) {
            if (c < r) {  /// we check whether appropriate dimensions are given or not.
                std::cout << "Dimension of vector is less than number of vector, hence "
                    "\n first "
                    << c << " vectors are orthogonalised\n";
                r = c;
            }

            int k = 1;

            while (k <= r) {
                if (k == 1) {
                    for (int j = 0; j < c; j++)
                        B[0][j] = A[0][j];  /// First vector is copied as it is.
                }

                else {
                    std::array<double, 30>
                        all_projection{};  /// array to store projections
                    for (int i = 0; i < c; ++i) {
                        all_projection[i] = 0;  /// First initialised to zero
                    }

                    int l = 1;
                    while (l < k) {
                        std::array<double, 30>
                            temp{};           /// to store previous projected array
                        double factor = NAN;  /// to store the factor by which the
                                              /// previous array will change
                        factor = projectionO1(A[k - 1], B[l - 1], c);
                        
                        register int i = 0;
                        int ostatak = c % 4;
                        int kolicnik = c / 4;
                        __m256 v2 = _mm256_set1_pd(factor);
                        if (kolicnik != 0) {
                            __m256d v1,suma;
                            for (; i < kolicnik; ++i) {
                                v1 = _mm256_load_pd(&B[l-1][i*4]);
                                suma = _mm256_load_pd(&all_projection[i*4]);
                                suma = _mm256_fmadd_pd(v1, v2, suma);
                                _mm256_store_pd(&all_projection[i*4], suma);
                            }
                        }

                        for(i = 4*kolicnik; i < c; ++i)
                        {
                            all_projection[i] += B[l - 1][i] * factor;
                        }

                        l++;
                    }
                    for (int i = 0; i < c; ++i) {
                        B[k - 1][i] =
                            A[k - 1][i] -
                            all_projection[i];  /// subtract total projection vector
                                                /// from the input vector
                    }
                }
                k++;
            }
            display(r, c, B);  // for displaying orthogoanlised vectors
        }

        void gram_schmidtO(int r, const int& c,
            const std::array<std::array<double, 30>, 30>& A,
            std::array<std::array<double, 30>, 30> B) {
            if (c < r) {  /// we check whether appropriate dimensions are given or not.
                std::cout << "Dimension of vector is less than number of vector, hence "
                    "\n first "
                    << c << " vectors are orthogonalised\n";
                r = c;
            }

            int k = 1;

            while (k <= r) {
                if (k == 1) {
                    for (int j = 0; j < c; j++)
                        B[0][j] = A[0][j];  /// First vector is copied as it is.
                }

                else {
                    std::array<double, 30>
                        all_projection{};  /// array to store projections
                    for (int i = 0; i < c; ++i) {
                        all_projection[i] = 0;  /// First initialised to zero
                    }

                    int l = 1;
                    while (l < k) {
                        std::array<double, 30>
                            temp{};           /// to store previous projected array
                        double factor = NAN;  /// to store the factor by which the
                                              /// previous array will change
                        factor = projectionO(A[k - 1], B[l - 1], c);
                        /*OPTIMIZOVANO*/
                       
                            register int i = 0;
                            
                                _mm_prefetch((char*)(&B[l - 1][i] + 0), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 1), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 2), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 3), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 4), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 5), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 6), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 7), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 0), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 1), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 2), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 3), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 4), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 5), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 6), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 7), _MM_HINT_T1);
                            
                            for (i = 1; i < c; ++i) {
                                _mm_prefetch((char*)(&B[l - 1][i] + 0), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 1), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 2), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 3), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 4), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 5), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 6), _MM_HINT_T1);
                                _mm_prefetch((char*)(&B[l - 1][i] + 7), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 0), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 1), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 2), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 3), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 4), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 5), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 6), _MM_HINT_T1);
                                _mm_prefetch((char*)(&all_projection[i] + 7), _MM_HINT_T1);
                                all_projection[i - 1] += B[l - 1][i - 1] * factor;
                            }
                            all_projection[c-1] += B[l - 1][c-1] * factor;
                        
                        l++;
                    }
                    for (int i = 0; i < c; ++i) {
                        B[k - 1][i] =
                            A[k - 1][i] -
                            all_projection[i];  /// subtract total projection vector
                                                /// from the input vector
                    }
                }
                k++;
            }
            displayO(r, c, B);  // for displaying orthogoanlised vectors
        }
    }  // namespace gram_schmidt
}  // namespace linear_algebra
/**
 * Test Function. Process has been tested for 3 Sample Inputs
 * @returns void
 */
static void test() {
    std::array<std::array<double, 30>, 30> a1 = {
        {{1, 0, 1, 0}, {1, 1, 1, 1}, {0, 1, 2, 1}} };
    std::array<std::array<double, 30>, 30> b1 = { {0} };
    double dot1 = 0;
    int flag = 1;
    StartTimer(ORIGINAL1)
    linear_algebra::gram_schmidt::gram_schmidt(3, 4, a1, b1);
    for (int i = 0; i < 2; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            dot1 = fabs(
                linear_algebra::gram_schmidt::dot_product(b1[i], b1[j], 4));
            if (dot1 > 0.1) {
                flag = 0;
                break;
            }
        }
    }
    if (flag == 0)
        std::cout << "Vectors are linearly dependent\n";
    assert(flag == 1);
    EndTimer

    a1 = { {{1, 0, 1, 0}, {1, 1, 1, 1}, {0, 1, 2, 1}} };
    b1 = { {0} };

   StartTimer(OPTIMIZOVANO1)
   linear_algebra::gram_schmidt::gram_schmidtO(3, 4, a1, b1);
    for (int i = 0; i < 2; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            dot1 = fabs(
                linear_algebra::gram_schmidt::dot_productO(b1[i], b1[j], 4));
            if (dot1 > 0.1) {
                flag = 0;
                break;
            }
        }
    }
    if (flag == 0)
        std::cout << "Vectors are linearly dependent\n";
    assert(flag == 1);
    EndTimer

    a1 = { {{1, 0, 1, 0}, {1, 1, 1, 1}, {0, 1, 2, 1}} };
    b1 = { {0} };

    StartTimer(OPTIMIZOVANOVEKT1)
        linear_algebra::gram_schmidt::gram_schmidtO1(3, 4, a1, b1);
    for (int i = 0; i < 2; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            dot1 = fabs(
                linear_algebra::gram_schmidt::dot_productO1(b1[i], b1[j], 4));
            if (dot1 > 0.1) {
                flag = 0;
                break;
            }
        }
    }
    if (flag == 0)
        std::cout << "Vectors are linearly dependent\n";
    assert(flag == 1);
    EndTimer
    
    std::cout << "Passed Test Case 1\n ";

    std::array<std::array<double, 30>, 30> a2 = { {{3, 1}, {2, 2}} };
    std::array<std::array<double, 30>, 30> b2 = { {0} };
    double dot2 = 0;
    StartTimer(ORIGINAL2)
    linear_algebra::gram_schmidt::gram_schmidt(2, 2, a2, b2);
    flag = 1;
    for (int i = 0; i < 1; ++i) {
        for (int j = i + 1; j < 2; ++j) {
            dot2 = fabs(
                linear_algebra::gram_schmidt::dot_product(b2[i], b2[j], 2));
            if (dot2 > 0.1) {
                flag = 0;
                break;
            }
        }
    }
    if (flag == 0)
        std::cout << "Vectors are linearly dependent\n";
    assert(flag == 1);
    EndTimer
    a2 = { {{3, 1}, {2, 2}} };
    b2 = { {0} };
    StartTimer(OPTIMIZOVANO2)
    linear_algebra::gram_schmidt::gram_schmidtO(2, 2, a2, b2);
    flag = 1;
    for (int i = 0; i < 1; ++i) {
        for (int j = i + 1; j < 2; ++j) {
            dot2 = fabs(
                linear_algebra::gram_schmidt::dot_productO(b2[i], b2[j], 2));
            if (dot2 > 0.1) {
                flag = 0;
                break;
            }
        }
    }
    if (flag == 0)
        std::cout << "Vectors are linearly dependent\n";
    assert(flag == 1);
    EndTimer
    a2 = { {{3, 1}, {2, 2}} };
    b2 = { {0} };
    StartTimer(OPTIMIZOVANOVEKT2)
        linear_algebra::gram_schmidt::gram_schmidtO1(2, 2, a2, b2);
    flag = 1;
    for (int i = 0; i < 1; ++i) {
        for (int j = i + 1; j < 2; ++j) {
            dot2 = fabs(
                linear_algebra::gram_schmidt::dot_productO1(b2[i], b2[j], 2));
            if (dot2 > 0.1) {
                flag = 0;
                break;
            }
        }
    }
    if (flag == 0)
        std::cout << "Vectors are linearly dependent\n";
    assert(flag == 1);
    EndTimer
    std::cout << "Passed Test Case 2\n";

    std::array<std::array<double, 30>, 30> a3;// = { {{1, 2, 2,4,5,6}, {-4, 3, 2,-1,-5,-6},{6,8,7,9,8,4},
                                              //{-1,2,-3,4,-5,6},{7,-3,-5,-1,0,2},{1,1,1,1,1,1} }};
    std::array<std::array<double, 30>, 30> b3;// = { {0} };
    for (int i = 0; i < 30; i++)
    {
        for (int j = 0; j < 30; j++)
        {
            a3[i][j] = rand()*10;
            b3[i][j] = 0;
        }
    }
    double dot3 = 0;

    StartTimer(ORIGINAL3)
    linear_algebra::gram_schmidt::gram_schmidt(30, 30, a3, b3);
    flag = 1;
    for (int i = 0; i < 29; ++i) {
        for (int j = i + 1; j < 29; ++j) {
            dot3 = fabs(
                linear_algebra::gram_schmidt::dot_product(b3[i], b3[j], 30));
            if (dot3 > 0.1) {
                flag = 0;
                break;
            }
        }
    }
    if (flag == 0)
        std::cout << "Vectors are linearly dependent\n";
    assert(flag == 1);
    EndTimer
        for (int i = 0; i < 30; i++)
        {
            for (int j = 0; j < 30; j++)
            {
                a3[i][j] = rand() * 10;
                b3[i][j] = 0;
            }
        }
    StartTimer(OPTIMIZOVANO3)
    linear_algebra::gram_schmidt::gram_schmidtO(30, 30, a3, b3);
    flag = 1;
    for (int i = 0; i < 29; ++i) {
        for (int j = i + 1; j < 29; ++j) {
            dot3 = fabs(
                linear_algebra::gram_schmidt::dot_productO(b3[i], b3[j], 30));
            if (dot3 > 0.1) {
                flag = 0;
                break;
            }
        }
    }
    if (flag == 0)
        std::cout << "Vectors are linearly dependent\n";
    assert(flag == 1);
    EndTimer
        for (int i = 0; i < 30; i++)
        {
            for (int j = 0; j < 30; j++)
            {
                a3[i][j] = rand() * 10;
                b3[i][j] = 0;
            }
        }
    StartTimer(OPTIMIZOVANOVEKT3)
        linear_algebra::gram_schmidt::gram_schmidtO1(30, 30, a3, b3);
    flag = 1;
    for (int i = 0; i < 29; ++i) {
        for (int j = i + 1; j < 29; ++j) {
            dot3 = fabs(
                linear_algebra::gram_schmidt::dot_productO1(b3[i], b3[j], 30));
            if (dot3 > 0.1) {
                flag = 0;
                break;
            }
        }
    }
    if (flag == 0)
        std::cout << "Vectors are linearly dependent\n";
    assert(flag == 1);
    EndTimer
    std::cout << "Passed Test Case 3\n";
}

/**
 * @brief Main Function
 * @return 0 on exit
 */
int main() {
    int r = 0, c = 0;
    test();  // perform self tests
    std::cout << "Enter the dimension of your vectors\n";
    std::cin >> c;
    std::cout << "Enter the number of vectors you will enter\n";
    std::cin >> r;

    std::array<std::array<double, 30>, 30>
        A{};  /// a 2-D array for storing all vectors
    std::array<std::array<double, 30>, 30> B = {
        {0} };  /// a 2-D array for storing orthogonalised vectors
    /// storing vectors in array A
    for (int i = 0; i < r; ++i) {
        std::cout << "Enter vector " << i + 1
            << '\n';  /// Input of vectors is taken
        for (int j = 0; j < c; ++j) {
            std::cout << "Value " << j + 1 << "th of vector: ";
            std::cin >> A[i][j];
        }
        std::cout << '\n';
    }

    StartTimer(ORIGINAL)

    linear_algebra::gram_schmidt::gram_schmidt(r, c, A, B);

    double dot = 0;
    int flag = 1;  /// To check whether vectors are orthogonal or  not

    
        for (int i = 0; i < r - 1; ++i) {
             for (int j = i + 1; j < r; ++j) {
                dot =
                 fabs(linear_algebra::gram_schmidt::dot_product(B[i], B[j], c));
                if (dot > 0.1)  /// take make the process numerically stable, upper
                            /// bound for the dot product take 0.1
                {
                    flag = 0;
                    break;
                }
             }
        }

        if (flag == 0)
            std::cout << "Vectors are linearly dependent\n";

        EndTimer

    StartTimer(OPTIMIZOVANO)

            linear_algebra::gram_schmidt::gram_schmidt(r, c, A, B);

       double dot = 0;
        int flag = 1;  /// To check whether vectors are orthogonal or  not


        for (int i = 0; i < r - 1; ++i) {
            for (int j = i + 1; j < r; ++j) {
                dot =
                    fabs(linear_algebra::gram_schmidt::dot_product(B[i], B[j], c));
                if (dot > 0.1)  /// take make the process numerically stable, upper
                            /// bound for the dot product take 0.1
                {
                    flag = 0;
                    break;
                }
            }
        }

        if (flag == 0)
            std::cout << "Vectors are linearly dependent\n";

        EndTimer

    return 0;
}
