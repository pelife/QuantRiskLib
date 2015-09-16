using System;
using System.Collections.Generic;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///
    ///References:
    ///Miller, Michael B. 2012. Mathematics and Statistics for Financial Risk Management. New York: John Wiley & Sons.
    ///Chapters 6 and 7

    public class Matrix
    {
        public readonly int NRows;              //the number of rows in the matrix
        public readonly int NColumns;           //the number of columns in the matrix
        private readonly double[][] _elements;    //the elements of the matrix. Even though all rows have the same lenght, [][] is faster than [,] for many computations.

        public new string ToString
        {
            get 
            { 
                string s = "";
                for(int r = 0; r < NRows; r++)
                {
                    s = s + "|" + _elements[r][0].ToString("0.00");
                    for(int c = 1; c < NColumns; c++)
                        s = s + ", " + _elements[r][c].ToString("0.00");
                    
                }
                s = s + "|";
                return s;
            }
        }

        #region Constructors
        /// <summary>
        /// Creates a new matrix with elements corresponding to the elements of the array.
        /// </summary>
        /// <param name="array">The elements of this array will be the initial elements of the matrix.</param>
        public Matrix(double[,] array)
        {
            NRows = array.GetLength(0);
            NColumns = array.GetLength(1);
            if(NRows == 0 || NColumns == 0)
                throw new ArgumentException("The elements of the array must have at least one row and one column");
            _elements = new double[NRows][];
            for (int r = 0; r < NRows; r++)
            {
                _elements[r] = new double[NColumns];
                for (int j = 0; j < NColumns; j++)
                    _elements[r][j] = array[r, j];
            }
        }

        /// <summary>
        /// Creates a vector, a matrix with one column.
        /// </summary>
        /// <param name="array">The elements of this array will be the initial elements of the matrix.</param>
        public Matrix(double[] array)
        {
            NRows = array.Length;
            if (NRows == 0)
                throw new ArgumentException("The array must have at least one element");
            NColumns = 1;
            _elements = new double[NRows][];
            for(int r = 0; r < NRows; r++)
            {
                _elements[r] = new double[1];
                _elements[r][0] = array[r];
            }
        }

        /// <summary>
        /// Creates a new matrix, nRows x nColumns, with all elements equal to zero.
        /// </summary>
        /// <param name="nRows">Number of rows.</param>
        /// <param name="nColumns">Number of columns.</param>
        public Matrix(int nRows, int nColumns)
        {
            NRows = nRows;
            NColumns = nColumns;
            _elements = new double[NRows][];
            for (int r = 0; r < NRows; r++)
                _elements[r] = new double[NColumns];
        }
        #endregion

        /// <summary>
        /// Creates an nxn identity matrix.
        /// </summary>
        /// <param name="n">Size of idenity matrix = number of rows = number of columns.</param>
        public static Matrix IdentityMatrix(int n)
        {
            Matrix M = new Matrix(n, n);
            for (int i = 0; i < n; i++)
                M[i, i] = 1.0;
            return M;
        }

        /// <summary>
        /// For an nxm matrix, M, the nxn centering matrix, C, has the property that each of the columns in CM is mean zero.
        /// </summary>
        public static Matrix CenteringMatrix(int n)
        {
            Matrix U = new Matrix(n, 1);
            for(int i = 0; i < n; i++)
                U[i, 0] = 1;

            Matrix V = U.Multiply(U.Transpose()).Multiply(1.0/n);
            Matrix C = IdentityMatrix(n).Subtract(V);
            return C;
        }

        /// <summary>
        /// Centers a matrix so that each of the columns is mean zero.
        /// </summary>
        public Matrix Center()
        {
            Matrix C = CenteringMatrix(this.NRows);
            return C.Multiply(this);
        }


        /// <summary>
        /// Returns element [r,c] from the matrix.
        /// </summary>
        /// <param name="r">row index</param>
        /// <param name="c">column index</param>
        public double this[int r, int c]
        {
            get { return _elements[r][c]; }
            set { _elements[r][c] = value; }
        }

        public double[,] ToArray()
        {
            double[,] array = new double[NRows,NColumns];
            for (int r = 0; r < NRows; r++)
                for(int c = 0; c < NColumns; c++)
                    array[r, c] = _elements[r][c];
            return array;
        }

        public double[] ToVectorArray()
        {
            if(NColumns != 1)
                throw new Exception("This method should only be used with vectors, matrices with only one column.");
            double[] array = new double[NRows];
            for (int r = 0; r < NRows; r++)
                array[r] = _elements[r][0];
            return array;
        }

        public Matrix Copy()
        {
            Matrix M = new Matrix(NRows, NColumns);
            for (int r = 0; r < NRows; r++)
                for (int c = 0; c < NColumns; c++)
                    M[r, c] = _elements[r][c];
            return M;
        }

        /// <summary>
        /// Returns the sum of matrix M and the existing matrix.
        /// </summary>
        /// <param name="M">Must be the same size as the existing matrix.</param>
        public Matrix Add(Matrix M)
        {
            if (NRows != M.NRows || NColumns != M.NColumns)
                throw new ArgumentException("Cannot add. The matrices must be the same size.");
            Matrix sum = new Matrix(NRows, NColumns);
            for (int r = 0; r < NRows; r++)
                for (int c = 0; c < NColumns; c++)
                    sum._elements[r][c] = _elements[r][c] + M._elements[r][c];
            return sum;
        }

        /// <summary>
        /// Matrix subtraction. X - M.
        /// </summary>
        /// <param name="M">The matrix to be subtracted.</param>
        public Matrix Subtract(Matrix M)
        {
            if (NRows != M.NRows || NColumns != M.NColumns)
                throw new ArgumentException("Cannot subtract. The matrices must be the same size.");
            Matrix S = new Matrix(NRows, NColumns);
            for (int r = 0; r < NRows; r++)
                for (int c = 0; c < NColumns; c++)
                    S._elements[r][c] = _elements[r][c] - M._elements[r][c];
            return S;
        }

        /// <summary>
        /// Returns the transpose of the existing matrix.
        /// </summary>
        public Matrix Transpose()
        {
            Matrix T = new Matrix(NColumns, NRows);
            for (int r = 0; r < NRows; r++)
                for (int c = 0; c < NColumns; c++)
                    T._elements[c][r] = _elements[r][c];
            return T;
        }

        /// <summary>
        /// Matrix multiplication. Post multiplies the existing matrix by another matrix, M.
        /// </summary>
        /// <param name="M">The number of rows in m should be equal to the number of columns in the existing matrix.</param>
        public Matrix Multiply(Matrix M)
        {
            if (NColumns != M.NRows)
                throw new ArgumentException("Cannot multiply. Number of columns in first matrix must equal the number of rows in the second.");

            Matrix product = new Matrix(NRows, M.NColumns);
            for (int r = 0; r < product.NRows; r++)
                for (int c = 0; c < product.NColumns; c++)
                {
                    double d = 0;
                    for (int i = 0; i < NColumns; i++)
                        d += _elements[r][i] * M._elements[i][c];   //NOTE: explicitly accessing the elements, M._elements[i, c], is faster than using get, M[i, c]. 
                    product._elements[r][c] = d;
                }
            return product;
        }

        /// <summary>
        /// Raises a matrix to the n-th power.
        /// </summary>
        public Matrix Power(int n)
        {
            Matrix M = Copy();
            for (int i = 0; i < n - 1; i++)
                M = M.Multiply(this);
            return M;
        }

        /// <summary>
        /// Scalar multiplication. Multiplies the existing matrix by d.
        /// </summary>
        /// <param name="d">The scalar by which the matrix is being multiplied.</param>
        /// <returns></returns>
        public Matrix Multiply(double d)
        {
            Matrix product = new Matrix(NRows, NColumns);
            for (int r = 0; r < product.NRows; r++)
                for (int c = 0; c < product.NColumns; c++)
                    product[r, c] = d * _elements[r][c];
            return product;
        }

        /// <summary>
        /// Returns true if matrix is square and X[r, c] = X[c, r] for all r and c.
        /// </summary>
        /// <returns></returns>
        public bool IsSymmetric()
        {
            if (NRows != NColumns) return false;
            for (int r = 0; r < NRows; r++)
                for (int c = 0; c < r; c++)
                    if (Math.Abs(_elements[r][c] - _elements[c][r]) > double.Epsilon)
                        return false;
            return true;
        }

        /// <summary>
        /// Kronecker product of two matrices, A and B.
        /// If A is m x n and B is p x q, then the return matrix will be mp x nq.
        /// </summary>
        public static Matrix KroneckerProduct(Matrix A, Matrix B)
        {
            Matrix M = new Matrix(A.NRows * B.NRows, A.NColumns * B.NColumns);
            for (int rowA = 0; rowA < A.NRows; rowA++)
            {
                for (int rowB = 0; rowB < B.NRows; rowB++)
                {
                    for (int colA = 0; colA < A.NColumns; colA++)
                    {
                        for (int colB = 0; colB < B.NColumns; colB++)
                        {
                            M[rowA * B.NRows + rowB, colA * B.NColumns + colB] = A[rowA, colA] * B[rowB, colB];
                        }
                    }
                }
            }
            return M;
        }

        /// <summary>
        /// Calculates the mean of all of the values in a matrix and then subtracts that mean from each entry.
        /// </summary>
        public Matrix DeviationsFromMean()
        {
            double mean = 0;
            for (int r = 0; r < NRows; r++)
                for (int c = 0; c < NColumns; c++)
                    mean += this[r, c];
            mean /= (NRows * NColumns);

            Matrix dMatrix = new Matrix(NRows, NColumns);
            for (int r = 0; r < NRows; r++)
                for (int c = 0; c < NColumns; c++)
                    dMatrix[r, c] = this[r, c] - mean;
            return dMatrix;
        }

        #region Vector math
        /// <summary>
        /// For a matrix M, often noted vec(M).
        /// If M is a matrix of column vectors, M = [c1 c2 c3], then vec(M) would be a vector with c1 stacked on top of c2 on top of c3.
        /// </summary>
        public Matrix Vectorization()
        {
            Matrix M = new Matrix(NRows * NColumns, 1);
            for (int c = 0; c < NColumns; c++)
                for (int r = 0; r < NRows; r++)
                    M[c * NRows + r, 0] = this[r, c];
            return M;
        }

        /// <summary>
        /// Inverse function of matrix vectorization.
        /// </summary>
        /// <param name="nColumns">The number of columns of the resulting matrix.</param>
        /// <returns></returns>
        public Matrix VectorizationInverse(int nColumns)
        {
            if(NColumns != 1)
                throw new ArgumentException("To perform VectorizationInverse matrix must be a vector");
            int rem;
            int nRows = Math.DivRem(NRows, nColumns, out rem);
            if (rem != 0)
                throw new ArgumentException("Cannot evenly divide this vector into this many columns");
            Matrix M = new Matrix(nRows, nColumns);
            for (int c = 0; c < nColumns; c++)
                for (int r = 0; r < nRows; r++)
                    M[r, c] = this[c * nRows + r, 0];
            return M;
        }

        /// <summary>
        /// Converts a matrix into an array of column matrices.
        /// </summary>
        public Matrix[] MatrixToColumMatrices()
        {
            Matrix[] cMatrices = new Matrix[NColumns];
            for (int i = 0; i < NColumns; i++)
            {
                cMatrices[i] = new Matrix(NRows, 1);
                for (int j = 0; j < NRows; j++)
                    cMatrices[i][j, 0] = _elements[j][i];
            }
            return cMatrices;
        }

        /// <summary>
        /// A.k.a. the dot product.
        /// a∙b = a[1] x b[1] + a[2] x b[2] + .... + a[n] x b[n] = a'b.
        /// a∙b = a'b.
        /// </summary>
        /// <param name="a">A vector, that is a matrix with only one column.</param>
        /// <param name="b">Another vector.</param>
        public static double InnerProduct(Matrix a, Matrix b)
        {
            if (a.NColumns != 1 || b.NColumns != 1)
                throw new ArgumentException("a and b should both be vectors (i.e. matrices with only one column).");
            if (a.NRows != b.NRows)
                throw new ArgumentException("a and b should be the same length.");

            double ip = 0;
            for (int r = 0; r < a.NRows; r++)
                ip += a._elements[r][0] * b._elements[r][0];
            return ip;
        }

        /// <summary>
        /// Project of b onto a = (a∙b/a∙a)a
        /// </summary>
        /// <param name="a">A vector.</param>
        /// <param name="b">A vector.</param>
        public static Matrix VectorProjection(Matrix a, Matrix b)
        {
            if (a.NColumns != 1 || b.NColumns != 1)
                throw new ArgumentException("a and b should both be vectors (i.e. matrices with only one column).");

            double d1 = InnerProduct(a, b);
            double d2 = InnerProduct(a, a);

            return a.Multiply(d1 / d2);
        }

        /// <summary>
        /// The magnitude (a.k.a. length, norm, or Euclidean norm). Often denoted ||v||.
        /// </summary>
        public double VectorMagnitude()
        {
            if (NColumns != 1)
                throw new ArgumentException("v should be a vector (i.e. a matrix with only one column)");

            double d = InnerProduct(this, this);
            return Math.Sqrt(d);
        }

        /// <summary>
        /// Transforms a vector into a square matrix with the vector elemetns on the diagonal and zero everywhere else.
        /// </summary>
        public Matrix VectorToDiagonalMatrix()
        {
            if (NColumns != 1)
                throw new ArgumentException("The input matrix for this function should be a vector, a matrix with only one column.");

            Matrix M = new Matrix(NRows, NRows);
            for (int i = 0; i < NRows; i++)
                M[i, i] = _elements[i][0];
            return M;
        }
        #endregion

        #region Matrix Decomposition
        /// <summary>
        /// Decomposes a symmetric square matrix into upper and lower triangular matrices, L and U, where L = U'.
        /// </summary>
        /// <param name="U">Upper triangular matrix</param>
        /// <param name="L">Lower triangular matrix</param>
        public void CholeskyDecomposition(out Matrix L, out Matrix U)
        {
            if (NRows != NColumns)
                throw new ArgumentException("Can only decompose square matrix. Number of rows and columns must be equal.");
            if (!IsSymmetric())
                throw new ArgumentException("Can only use Cholesky with a symmetric matrix.");
            int n = NRows;
            L = new Matrix(n, n);

            for(int r = 0; r < n; r++)
                for (int c = 0; c <= r; c++)
                {
                    if (r == c)
                    {
                        double d = this[r, c];
                        for (int i = 0; i < r; i++)
                            d -= L[r, i] * L[r, i];
                        L[r, c] = Math.Sqrt(d);
                    }
                    else
                    {
                        double d = this[r, c];
                        for (int i = 0; i < c; i++)
                            d -= L[r, i] * L[c, i];
                        d /= L[c, c];
                        L[r, c] = d;
                    }
                }

            U = L.Transpose();
        }

        /// <summary>
        /// Decomposes a square matrix into upper and lower triangular matrices, L and U, such that X = LU. Crout method.
        /// </summary>
        /// <param name="U">Upper triangular matrix</param>
        /// <param name="L">Lower triangular matrix</param>
        /// <param name="checkDecomposition">If true, will throw exception if matrix is not *exactly* equal to LU.</param>
        public void LuDecomposition(out Matrix L, out Matrix U, bool checkDecomposition = false)
        {
            if (NRows != NColumns)
                throw new ArgumentException("Can only decompose square matrix. Number of rows and columns must be equal.");

            int n = NRows;
            L = new Matrix(n, n);
            U = new Matrix(n, n);

            for (int i = 0; i < n; i++)
                U._elements[i][i] = 1;

            for (int j = 0; j < n; j++)
            {
                for (int i = j; i < n; i++)
                {
                    double sum = 0;
                    for (int k = 0; k < j; k++)
                        sum = sum + L._elements[i][k]*U._elements[k][j];

                    L._elements[i][j] = _elements[i][j] - sum;
                }

                for (int i = j; i < n; i++)
                {
                    double sum = 0;
                    for (int k = 0; k < j; k++)
                        sum = sum + L._elements[j][k]*U._elements[k][i];

                    if (L._elements[j][j] == 0)
                        throw new Exception("LuDecomposition: cannot divide by 0");
                    
                    U._elements[j][i] = (_elements[j][i] - sum) / L._elements[j][j];
                }
            }

            //check decomposition
            if (checkDecomposition)
            {
                Matrix C = this.Subtract(L.Multiply(U));
                double sumAbsDiff = C.SumAbsoluteValueOfElements();
                if (sumAbsDiff != 0)
                    throw new Exception("LuDecomposition check failed.");
            }
        }

        /// <summary>
        /// The basic Gram–Schmidt process for QR decomposition of a matrix.
        /// M is decomposed into matrices Q and R, such that M = QR, Q is orthogonal and R is upper triangular. 
        /// Warning: this method is subject to rounding error and is not efficient for large matrices.
        /// </summary>
        /// <param name="Q">An orthogonal matrix.</param>
        /// <param name="R">An upper triangular matrix.</param>
        public void GramSchmidtQrDecomposition(out Matrix Q, out Matrix R)
        {
            if (NColumns != NRows)
                throw new ArgumentException("This algorithm currently only works for square matrices.");

            int n = NColumns;
            Matrix[] cMatrices = MatrixToColumMatrices();
            Matrix[] uMatrices = new Matrix[n];
            Matrix[] eMatrices = new Matrix[n];

            for (int i = 0; i < n; i++)
            {
                uMatrices[i] = cMatrices[i].Copy();
                for (int j = 0; j < i; j++)
                {
                    uMatrices[i] = uMatrices[i].Subtract(VectorProjection(eMatrices[j], cMatrices[i]));
                }
                eMatrices[i] = uMatrices[i].Multiply(1 / uMatrices[i].VectorMagnitude());
            }

            Q = new Matrix(n, n);
            for (int c = 0; c < n; c++)
                for (int r = 0; r < n; r++)
                    Q[r, c] = eMatrices[c][r, 0];

            R = new Matrix(n, n);
            for (int c = 0; c < n; c++)
            {
                for (int r = 0; r < n; r++)
                {
                    if (c < r)
                        R[r, c] = 0;
                    else
                        R[r, c] = InnerProduct(eMatrices[r], cMatrices[c]);
                }
            }
        }

        /// <summary>
        /// The Modified Gram–Schmidt process for QR decomposition of a matrix.
        /// M is decomposed into matrices Q and R, such that M = QR, Q is orthogonal and R is upper triangular. 
        /// Warning: this method should be more accurate than the plain GramSchmidt decomposition but is still subject to rounding error.
        /// </summary>
        /// <param name="Q">An orthogonal matrix.</param>
        /// <param name="R">An upper triangular matrix.</param>
        public void ModifiedGramSchmidtQrDecomposition(out Matrix Q, out Matrix R)
        {
            if (NColumns != NRows)
                throw new ArgumentException("This algorithm currently only works for square matrices.");

            int n = NColumns;
            Matrix[] cMatrices = MatrixToColumMatrices();
            Matrix[] uMatrices = new Matrix[n];
            Matrix[] eMatrices = new Matrix[n];

            for (int i = 0; i < n; i++)
            {
                if(i == 0)
                    uMatrices[i] = cMatrices[i].Copy();
                else
                    uMatrices[i] = cMatrices[i].Subtract(VectorProjection(uMatrices[0], cMatrices[i]));

                for (int j = 1; j < i; j++)
                {
                    uMatrices[i] = uMatrices[i].Subtract(VectorProjection(uMatrices[j], uMatrices[i]));
                }
                eMatrices[i] = uMatrices[i].Multiply(1 / uMatrices[i].VectorMagnitude());
            }

            Q = new Matrix(n, n);
            for (int c = 0; c < n; c++)
                for (int r = 0; r < n; r++)
                    Q[r, c] = eMatrices[c][r, 0];

            R = new Matrix(n, n);
            for (int c = 0; c < n; c++)
            {
                for (int r = 0; r < n; r++)
                {
                    if (c < r)
                        R[r, c] = 0;
                    else
                        R[r, c] = InnerProduct(eMatrices[r], cMatrices[c]);
                }
            }
        }

        /// <summary>
        /// Produces a tridiagonal matrix using the Householder transformation, H = QAQ'.
        /// </summary>
        /// <param name="Q">H = QAQ', where A is input matrix and B is output matrix.</param>
        public Matrix HouseholderTridiagonalization(out Matrix Q)
        {
            Matrix H = HouseholderTridiagonalizationStep(1, out Q);
            int n = this.NRows - 1;
            for (int i = 2; i < n; i++)
            {
                Matrix P;
                H = H.HouseholderTridiagonalizationStep(i, out P);
                Q = P.Multiply(Q);
            }

            //The following loop is not necessary in theory, but it ensures the resulting matrix is symmetric (could have precision issues).
            for (int r = 1; r < NRows; r++)
            {
                double d = 0.5 * (H._elements[r - 1][r] + H._elements[r][r - 1]);
                H._elements[r - 1][r] = d;
                H._elements[r][r - 1] = d;
            }

            return H;
        }

        /// <param name="k">k = 1, 2, ..., n</param>
        /// <param name="P">B = PAP', where A is input matrix and B is output matrix.</param>
        private Matrix HouseholderTridiagonalizationStep(int k, out Matrix P)
        {
            double a = 0;
            for (int r = k; r < NRows; r++)
                a += _elements[r][k - 1] * _elements[r][k - 1];
            a = -Math.Sign(_elements[k][k - 1]) * Math.Sqrt(a);
            double c = 1.0 / (a * a - _elements[k][k - 1] * a);

            Matrix v = new Matrix(NRows, 1);
            for (int r = 0; r < k - 1; r++)
                v._elements[r][0] = 0;
            v._elements[k][0] = (_elements[k][k - 1] - a);
            for (int r = k + 1; r < NRows; r++)
                v._elements[r][0] = _elements[r][k - 1];
            
            Matrix N = v.Multiply(v.Transpose()).Multiply(c);
            P = IdentityMatrix(NRows).Subtract(N);
            Matrix H = P.Multiply(this).Multiply(P);

            //In theory this is not necessary, but precision errors can cause these values to differ slightly from zero.
            for (int r = k + 1; r < NRows; r++)
            {
                H._elements[r][k - 1] = 0;
                H._elements[k - 1][r] = 0;
            }

            return H;
        }
        #endregion

        #region Matrix Inversion
        /// <summary>
        /// Used for regression analysis.
        /// Calculates the inverse of the existing matrix. Matrix must be square.
        /// If the matrix contains a diagonal element equal to zero, returns the inverse matrix of the remaining rows and columns, with zero rows and columns inserted.
        /// A diagonal element equal to zero happens when one of the columns of the X matrix is all zeros.
        /// </summary>
        public Matrix PseudoInverse()
        {
            if (NRows != NColumns)
                throw new ArgumentException("Can only invert square matrix. Number of rows and columns must be equal.");
            int nBig = NRows;

            List<int> zeroIndexes = new List<int>();
            for (int i = 0; i < nBig; i++)
                if(_elements[i][i] == 0)
                    zeroIndexes.Add(i);

            if (zeroIndexes.Count == 0)
                return Inverse();

            int nSmall = nBig - zeroIndexes.Count;
            Matrix M = new Matrix(nSmall, nSmall);
            int newRowIndex, newColIndex = 0;
            for (int c = 0; c < NColumns; c++)
            {
                if(zeroIndexes.Contains(c)) continue;
                newRowIndex = 0;
                for (int r = 0; r < NRows; r++)
                {
                    if (zeroIndexes.Contains(r)) continue;
                    M._elements[newRowIndex][newColIndex] = _elements[r][c];
                    newRowIndex++;
                }
                newColIndex++;
            }

            Matrix MInv = M.Inverse();
            Matrix Inv = new Matrix(NRows, NColumns);
            newRowIndex = 0;
            for (int r = 0; r < NRows; r++)
            {
                if (zeroIndexes.Contains(r)) continue;  //leave row all zeros.
                newColIndex = 0;
                for (int c = 0; c < NColumns; c++)
                {
                    if (zeroIndexes.Contains(c))
                        Inv._elements[r][c] = 0;
                    else
                    {
                        Inv._elements[r][c] = MInv._elements[newRowIndex][newColIndex];
                        newColIndex++;
                    }
                }
                newRowIndex++;
            }
            return Inv;
        }

        /// <summary>
        /// Calculates the inverse of the existing matrix. Matrix must be square.
        /// </summary>
        public Matrix Inverse()
        {
            //This method first decomposes the existing matrix into a lower and upper tringular matrix. 
            //These can be more easily inverted.
            //There may be more precise or efficient algorithms, but this is very easy to follow.

            if (NRows != NColumns)
                throw new ArgumentException("Can only invert square matrix. Number of rows and columns must be equal.");

            Matrix L, U;
            LuDecomposition(out L, out U);

            Matrix Linv = L.LowerTriangularInverse();
            Matrix Uinv = U.Transpose().LowerTriangularInverse().Transpose();

            Matrix inverse = Uinv.Multiply(Linv);

            return inverse;
        }

        /// <summary>
        /// Checks that the inverse of a matrix multiplied by the matrix returns the idenity matrix. 
        /// If M and N are inverses, calculates C = MN - I, and the average absolute value of all the elements in C. 
        /// In theory, returned value should be zero. In practice, floating-point precision errors can cause result to be non-zero.
        /// </summary>
        public double InverseMeanError(Matrix Inverse)
        {
            if (NRows != NColumns)
                throw new ArgumentException("Can only invert square matrix. Number of rows and columns must be equal.");
            Matrix E1 = this.Multiply(Inverse);
            Matrix E2 = Inverse.Multiply(this);
            for (int r = 0; r < NRows; r++)
            {
                E1._elements[r][r] -= 1.0;
                E2._elements[r][r] -= 1.0;
            }
            return (E1.SumAbsoluteValueOfElements() + E2.SumAbsoluteValueOfElements())/(2.0*NRows*NRows);
        }

        /// <summary>
        /// The sum of the absolute values of all of the elements in a matrix.
        /// </summary>
        public double SumAbsoluteValueOfElements()
        {
            double sum = 0;
            for(int r = 0; r<NRows; r++)
                for (int c = 0; c < NColumns; c++)
                {
                    sum += Math.Abs(_elements[r][c]);
                }
            return sum;
        }

        /// <summary>
        /// Returns true if the matrix is lower triangular.
        /// Matrix must be square. All entries above the diagonal are zero.
        /// </summary>
        /// <returns></returns>
        public bool IsLowerTriangular()
        {
            return Transpose().IsUpperTriangular();
        }

        /// <summary>
        /// Returns true if the matrix is upper triangular.
        /// Matrix must be square. All entries below the diagonal are zero.
        /// </summary>
        public bool IsUpperTriangular()
        {
            if (NRows != NColumns)
                throw new ArgumentException("Upper and lower triangular only defined for square matrices. Number of rows and columns must be equal.");
            int n = NRows;
            for (int r = 1; r < n; r++)
            {
                for (int c = 0; c < r; c++)
                    if (this[r, c] != 0)
                        return false;
            }
            return true;
        }

        /// <summary>
        /// Calculates the inverse of a lower triangular matrix. Used to calculate general inverses.
        /// </summary>
        /// <returns></returns>
        private Matrix LowerTriangularInverse()
        {
            if (!IsLowerTriangular())
                throw new Exception("Can only use with lower triangular matrix");
            int n = NRows;

            Matrix D = new Matrix(n, n);
            for (int i = 0; i < n; i++)
                D[i, i] = 1.0 / _elements[i][i];
            Matrix C = D.Multiply(this);

            //needed to prevent occasional rounding error
            for (int i = 0; i < n; i++)
                C[i, i] = 1.0; 
            
            Matrix Cinv = C.LowerTriangularInverseOneDiagonal();
            Matrix M = Cinv.Multiply(D);
            return M;
        }

        /// <summary>
        /// Used to calculate LowerTriangularInverse, which is used to calculate inverses.
        /// </summary>
        private Matrix LowerTriangularInverseOneDiagonal()
        {
            //The inverse of a lower triangular matrix with ones on the diagonal is also
            //lower triangular with ones on the diagonal.
            //This problem is relatively easy to solve and can be used to solve the more
            //general inversion problem.

            if (!IsLowerTriangular())
                throw new Exception("Can only use with lower triangular matrix");
            int n = NRows;
            for (int i = 0; i < n; i++)
                if (_elements[i][i] != 1)
                    throw new Exception("Can only use where all diagonal elements are 1");

            Matrix M = new Matrix(n, n);
            for (int i = 0; i < n; i++)
                M[i, i] = 1;

            for(int d = 1; d < n; d++)
                for (int r = d; r < n; r++)
                {
                    int c = r - d;
                    for (int i = 0; i < r; i++)
                        M._elements[r][c] -= _elements[r][i] * M._elements[i][c];
                }

            return M;
        }
        #endregion

        #region Eigenvales and Eigenvectors

        /// <summary>
        /// Returns the eigenvalues of a matrix in a vector.
        /// If the matrix is symmetric the eigenvectors are also retured as the columns of a matrix; otherwise this matrix is null.
        /// Based on Householder transformation and QR algorithm.
        /// </summary>
        /// <param name="eigenVectors"></param>
        public Matrix Eigenvalues(out Matrix eigenVectors)
        {
            if (!IsSymmetric())
                throw new ArgumentException("This algorithm only works for symmetric matrices.");

            Matrix Q;
            Matrix H = HouseholderTridiagonalization(out Q);

            Matrix P;
            Matrix B = H.QRAlgorithmWithEigenVectors(out P, this, Q.Transpose());

            eigenVectors = Q.Transpose().Multiply(P);

            Matrix E = new Matrix(NRows, 1);
            for (int i = 0; i < NRows; i++)
                E[i, 0] = B[i, i];

            OrderEigenVectorsAndValues(ref E, ref eigenVectors);

            return E;
        }

        /// <summary>
        /// Returns the eigenvalues of a matrix in a vector.
        /// If the matrix is symmetric the eigenvectors are also retured as the columns of a matrix; otherwise this matrix is null.
        /// Based only on QR algorithm.
        /// WARNING: This method is easy to follow, but may not be very accurate and may not converge.
        /// </summary>
        /// <param name="eigenVectors"></param>
        public Matrix EigenvaluesSimple(out Matrix eigenVectors)
        {
            if (NColumns != NRows)
                throw new ArgumentException("This algorithm currently only works for square matrices.");

            Matrix B;
            if (IsSymmetric())
            {
                B = QRAlgorithmWithEigenVectors(out eigenVectors);
            }
            else
            {
                eigenVectors = null;
                B = QRAlgorithm();
            }

            Matrix E = new Matrix(NRows, 1);
            for (int i = 0; i < NRows; i++)
                E[i, 0] = B[i, i];

            OrderEigenVectorsAndValues(ref E, ref eigenVectors);

            return E;
        }

        /// <summary>
        /// The QR algorithm for calculating eigenvalues. This version will work on any square matrix.
        /// </summary>
        public Matrix QRAlgorithm()
        {
            if (NColumns != NRows)
                throw new ArgumentException("This algorithm currently only works for square matrices.");

            //The algorithm will stop if:
            //1) It completes maxLoops number of loops.
            //2) The max absolute percentage change in any diagnonal element is less than maxChangeThreshold.
            const int maxLoops = 64;
            const double maxChangeThreshold = 0.001;

            Matrix D = Copy();
            double[] oldDiagonalElements = DiagonalElementsToArray(D);
            bool converged = false;
            for (int i = 0; i < maxLoops; i++)
            {
                Matrix Q, R;
                D.ModifiedGramSchmidtQrDecomposition(out Q, out R);
                D = R.Multiply(Q);

                //check for convergence
                double[] diagonalElements = DiagonalElementsToArray(D);
                double maxChange = MaxAbsolutePercentChangeInElements(oldDiagonalElements, diagonalElements);
                if(maxChange < maxChangeThreshold)
                {
                    converged = true;
                    break;
                }

                oldDiagonalElements = diagonalElements;
            }

            if (!converged)
                throw new Exception(string.Format("QR Algorithm failed to converge after {0} loops.", maxLoops));
            return D;
        }

        /// <summary>
        /// The QR algorithm for calculating eigenvalues and eigenvectors. This will only work for symmetric matrices.
        /// </summary>
        /// <param name="eigenVectors"></param>
        /// <param name="testTarget"></param>
        /// <param name="testEigenVectorPreMultiply"></param>
        public Matrix QRAlgorithmWithEigenVectors(out Matrix eigenVectors, Matrix testTarget = null, Matrix testEigenVectorPreMultiply = null)
        {
            if (!IsSymmetric())
                throw new ArgumentException("This algorithm only works for symmetric matrices.");

            //The algorithm will stop if:
            //1) It completes maxLoops number of loops.
            //2) The accuracy is less than accuracyThreshold.
            const int maxLoops = 128;
            const double accuracyThreshold = 0.01;

            Matrix E = Copy();
            eigenVectors = IdentityMatrix(E.NRows);
            bool converged = false;
            for (int i = 0; i < maxLoops; i++)
            {
                Matrix Q, R;
                E.ModifiedGramSchmidtQrDecomposition(out Q, out R);
                E = R.Multiply(Q);
                eigenVectors = eigenVectors.Multiply(Q);

                //check for convergence every 4th loop
                if(i % 4 != 0) continue;
                double accuracy;
                Matrix testEigenValues = E.DiagonalElements();
                if(testTarget != null && testEigenVectorPreMultiply != null)
                {
                    Matrix testEigenVectors = testEigenVectorPreMultiply.Multiply(eigenVectors);
                    accuracy = testTarget.EigenAccuracy(testEigenVectors, testEigenValues);
                }
                else
                {
                    accuracy = this.EigenAccuracy(eigenVectors, testEigenValues);
                }
                if(accuracy < accuracyThreshold)
                {
                    converged = true;
                    break;
                }
            }

            if (!converged)
                throw new Exception(string.Format("QR Algorithm failed to converge after {0} loops.", maxLoops));
            return E;
        }

        private Matrix DiagonalElements()
        {
            int n = Math.Min(NRows, NColumns);
            Matrix D = new Matrix(n, n);
            for(int r = 0; r < n; r++)
                D._elements[r][r] = _elements[r][r];
            return D;
        }

        private static double MaxAbsolutePercentChangeInElements(double[] oldArray, double[] newArray)
        {
            if (oldArray.Length != newArray.Length)
                throw new Exception("Arrays should be of equal length.");

            double maxChange = 0;
            for (int i = 0; i < oldArray.Length; i++)
            {
                double change = newArray[i] == oldArray[i] ? 0 : newArray[i] / oldArray[i] - 1;
                if (change > maxChange) maxChange = change;
            }
            return maxChange;
        }

        private static double[] DiagonalElementsToArray(Matrix M)
        {
            if (M.NColumns != M.NRows)
                throw new ArgumentException("This algorithm currently only works for square matrices.");

            double[] da = new double[M.NRows];
            for (int i = 0; i < M.NRows; i++)
                da[i] = M[i, i];
            
            return da;
        }

        private void OrderEigenVectorsAndValues(ref Matrix eigenValues, ref Matrix eigenVectors)
        {
            if (eigenVectors == null)
            {
                List<double> dList = new List<double>();
                dList.AddRange(eigenValues.ToVectorArray());
                dList.Sort();
                dList.Reverse();
                for (int r = 0; r < eigenValues.NRows; r++)
                    eigenValues[r, 0] = dList[r];
                return;
            }


            Matrix[] vectorColumns = eigenVectors.MatrixToColumMatrices();
            SortedDictionary<double, double[]> dict = new SortedDictionary<double, double[]>();
            for (int r = 0; r < eigenValues.NRows; r++)
            {
                dict.Add(eigenValues[r, 0], vectorColumns[r].ToVectorArray());
            }

            int c = eigenVectors.NColumns;
            foreach (KeyValuePair<double, double[]> pair in dict)
            {
                c--;
                for (int r = 0; r < eigenValues.NRows; r++)
                {
                    eigenVectors[r, c] = pair.Value[r];
                }
                eigenValues[c, 0] = pair.Key;
            }
        }

        /// <summary>
        /// For the current matrix, A, with an eigenvector matrix P and a diagonal matrix containing the eigenvectors, D, we should have: 
        /// A = PDP'. This method returns the average absolute percentage difference between elements of A and PDP'.
        /// </summary>
        /// <param name="P">Matrix containing eigenvectors in columns.</param>
        /// <param name="D">Matrix with corresponding eigenvalues on diagonal.</param>
        public double EigenAccuracy(Matrix P, Matrix D)
        {
            Matrix C = P.Multiply(D).Multiply(P.Transpose());
            return MeanAbsoluteDifference(this, C, true);
        }

        /// <summary>
        /// This method returns the average absolute difference between elements of T and E.
        /// If asPercent: |(E[i,j] - T[i,j])/T[i,j]|; otherwise: |E[i,j] - T[i,j]|.
        /// </summary>
        /// <param name="T">Target</param>
        /// <param name="E">Estimate</param>
        /// <param name="asPercent"></param>
        private static double MeanAbsoluteDifference(Matrix T, Matrix E, bool asPercent)
        {
            double meanAbsDiff = 0;
            if(asPercent)
            {
                for(int r = 0; r < T.NRows; r++)
                    for(int c = 0; c < T.NColumns; c++)
                    {
                        if (E._elements[r][c] != T._elements[r][c])
                            meanAbsDiff += Math.Abs(E._elements[r][c] / T._elements[r][c] - 1);
                    }          
            }
            else
            {
                for (int r = 0; r < T.NRows; r++)
                    for (int c = 0; c < T.NColumns; c++)
                        meanAbsDiff += Math.Abs(E._elements[r][c] - T._elements[r][c]);
            }
            
            meanAbsDiff /= (T.NRows * T.NColumns);
            return meanAbsDiff;
        }

        #endregion

        #region Gaussian Elimination

        /// <summary>
        /// Produces a matrix that is the result of Gaussian elimination on the current matrix.
        /// </summary>
        public Matrix GaussianElimination()
        {
            Matrix A = Copy();
            int nRows = A.NRows;

            for (int i = 0; i < nRows; i++)
            {
                //This swapping is not necessary in theory, but it should help reduce rounding error.
                int maxRow = i;
                double maxAbsPivot = Math.Abs(A[i, i]);
                for (int j = i + 1; j < nRows; j++)
                {
                    if (Math.Abs(A[j, i]) > maxAbsPivot)
                    {
                        maxAbsPivot = Math.Abs(A[j, i]);
                        maxRow = j;
                    }
                }
                A = A.SwapRows(i, maxRow);

                A = A.MultiplyRowByScalar(i, 1.0 / A[i, i]);
                for (int j = 0; j < nRows; j++)
                {
                    if (j == i) continue;
                    A = A.AddMultipleOfRowToAnotherRow(i, j, -A[j, i]);
                    A[j, i] = 0; //This is also not necessary in theory, but should help with rounding error.
                }
            }
            return A;
        }

        private Matrix SwapRows(int row1Index, int row2Index)
        {
            Matrix C = Copy();
            for (int i = 0; i < NColumns; i++)
            {
                C[row1Index, i] = _elements[row2Index][i];
                C[row2Index, i] = _elements[row1Index][i];
            }
            return C;
        }

        private Matrix MultiplyRowByScalar(int rowIndex, double scalar)
        {
            Matrix C = Copy();
            for (int i = 0; i < NColumns; i++)
                C[rowIndex, i] = scalar * _elements[rowIndex][i];
            return C;
        }

        private Matrix AddMultipleOfRowToAnotherRow(int rowIndexBeingAdded, int rowIndexToAddTo, double multiple)
        {
            Matrix C = Copy();
            for (int i = 0; i < NColumns; i++)
                C[rowIndexToAddTo, i] = _elements[rowIndexToAddTo][i] + multiple * _elements[rowIndexBeingAdded][i];
            return C;
        }
        #endregion
    }

}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.