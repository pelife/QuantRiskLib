using System;

namespace QuantRiskLib
{
    public class SerialCorrelationCorrection
    {
        /// <summary>
        /// Ratio of variance given serial correlation, to variance with no serial correlation.
        /// Multiply the uncorrected value by this factor.
        /// </summary>
        public static double SerialCorrelationVarianceCorrection(double serialCorrelation, int windowLength)
        {
            if (serialCorrelation < -1 || serialCorrelation > 1)
                throw new ArgumentException("Serial corrleation must be between -1 and +1");
            if (windowLength < 0)
                throw new ArgumentException("WindowLength must be positive");
            double d = 1.0 + 2 * serialCorrelation / (1 - serialCorrelation) - 2 * serialCorrelation * (1 - Math.Pow(serialCorrelation, windowLength)) / (windowLength * (1 - serialCorrelation) * (1 - serialCorrelation));
            return d;
        }

        /// <summary>
        /// Serial correlation correction for infinite window length.
        /// Ratio of variance given serial correlation, to variance with no serial correlation.
        /// </summary>
        public static double SerialCorrelationVarianceCorrection(double serialCorrelation)
        {
            if (serialCorrelation < -1 || serialCorrelation > 1)
                throw new ArgumentException("Serial corrleation must be between -1 and +1");
            double d = 1.0 + 2 * serialCorrelation / (1 - serialCorrelation);
            return d;
        }

        /// <summary>
        /// Premultiply the no serial correlation covariance matrix by this matrix to get the adjusted matrix.
        /// </summary>
        /// <param name="D">Matrix of decay factors from VAR(1) model.</param>
        /// <param name="n">Number of periods.</param>
        /// <param name="covarR">One-period covariance matrix.</param>
        public static Matrix SerialCorrelationVarianceCorrection(Matrix D, int n, Matrix covarR)
        {
            Matrix A = ABTransform(D, n);
            Matrix B = ABTransform(D.Transpose(), n);
            Matrix I2 = Matrix.IdentityMatrix(2);
            Matrix I4 = Matrix.IdentityMatrix(4);
            Matrix C = I4.Multiply(n).Add(Matrix.KroneckerProduct(I2, A)).Add(Matrix.KroneckerProduct(B.Transpose(), I2));

            Matrix vecCR = covarR.Vectorization();
            Matrix vecYn = C.Multiply(vecCR);
            Matrix Yn = vecYn.VectorizationInverse(2);
            Matrix M = Yn.Multiply(covarR.Inverse()).Multiply(1.0 / n);
            return M;
        }

        /// <summary>
        /// Premultiply the no serial correlation covariance matrix by this matrix to get the adjusted matrix.
        /// This version is the correction for a window of infinite length.
        /// </summary>
        /// <param name="D">Matrix of decay factors from VAR(1) model.</param>
        /// <param name="covarR">One-period covariance matrix.</param>
        public static Matrix SerialCorrelationVarianceCorrection(Matrix D, Matrix covarR)
        {
            Matrix A = ABTransform(D);
            Matrix B = ABTransform(D.Transpose());
            Matrix I2 = Matrix.IdentityMatrix(2);
            Matrix I4 = Matrix.IdentityMatrix(4);
            Matrix C = I4.Add(Matrix.KroneckerProduct(I2, A)).Add(Matrix.KroneckerProduct(B.Transpose(), I2));

            Matrix vecCR = covarR.Vectorization();
            Matrix vecYn = C.Multiply(vecCR);
            Matrix Yn = vecYn.VectorizationInverse(2);
            Matrix M = Yn.Multiply(covarR.Inverse());
            return M;
        }

        private static Matrix ABTransform(Matrix D, int n)
        {
            Matrix I = Matrix.IdentityMatrix(2);
            Matrix iMinusDInv = (I.Subtract(D)).Inverse();
            Matrix A = D.Multiply(iMinusDInv).Multiply(n);

            Matrix iMinusDn = I.Subtract(D.Power(n));
            Matrix B = D.Multiply(iMinusDn).Multiply(iMinusDInv).Multiply(iMinusDInv);
            return A.Subtract(B);
        }

        private static Matrix ABTransform(Matrix D)
        {
            Matrix I = Matrix.IdentityMatrix(2);
            Matrix iMinusDInv = (I.Subtract(D)).Inverse();
            Matrix A = D.Multiply(iMinusDInv);
            return A;
        }

    }

    public partial class Moments
    {
        /// <summary>
        /// Returns the sample standard deviation of an array, corrected for serial correlation.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the variance.  For time series, the last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in array is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        /// <param name="correctSerialCorrelation">If true, correction is for infinite return length.</param>
        public static double StandardDeviation(double[] array, double decayFactor, int length, bool correctSerialCorrelation)
        {
            double uncorrected = StandardDeviation(array, decayFactor, length);
            if (!correctSerialCorrelation)
                return uncorrected;

            double[] subArray = Tools.MostRecentValues(array, length);
            double serialCorrelation = SerialCorrelation(subArray, decayFactor);
            double varianceCorrectionFactor = SerialCorrelationCorrection.SerialCorrelationVarianceCorrection(serialCorrelation); //don't use version w/ windowLength. windowLenght != length. Length is how far back we go to get data for the calculation. windowLenght is the return lenght that we are correcting to (infinite, by default).
            return uncorrected * Math.Sqrt(varianceCorrectionFactor);
        }

        /// <summary>
        /// Returns the sample covariance between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// </summary>
        /// <param name="array1"></param>
        /// <param name="array2"></param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="correctSerialCorrelation">If true, correction is for infinite return length.</param>
        public static double Covariance(double[] array1, double[] array2, double decayFactor, bool correctSerialCorrelation)
        {
            if(!correctSerialCorrelation)
                return Covariance(array1, array2, decayFactor);

            if (Tools.ArrayAllEqual(array1) || Tools.ArrayAllEqual(array2))
                return 0.0;

            int len = Math.Min(array1.Length, array2.Length) - 1;
            Matrix C = CovarianceMatrixSerialCorrelationCorrected(array1, array2, decayFactor, len);
            return C[0, 1];
        }


        /// <summary>
        /// Returns the correlation between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// If the variance of either array is zero, returns zero.
        /// </summary>
        /// <param name="array1"></param>
        /// <param name="array2"></param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="correctSerialCorrelation">If true, correction is for infinite return length.</param>
        /// /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        public static double Correlation(double[] array1, double[] array2, double decayFactor, bool correctSerialCorrelation, int length)
        {
            double[] subArray1 = Tools.MostRecentValues(array1, length);
            double[] subArray2 = Tools.MostRecentValues(array2, length);
            return Correlation(subArray1, subArray2, decayFactor, correctSerialCorrelation);
        }
        

        /// <summary>
        /// Returns the correlation between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// If the variance of either array is zero, returns zero.
        /// </summary>
        /// <param name="array1"></param>
        /// <param name="array2"></param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="correctSerialCorrelation">If true, correction is for infinite return length.</param>
        public static double Correlation(double[] array1, double[] array2, double decayFactor, bool correctSerialCorrelation)
        {
            if (!correctSerialCorrelation)
                return Correlation(array1, array2, decayFactor);

            if (Tools.ArrayAllEqual(array1) || Tools.ArrayAllEqual(array2))
                return 0.0;
            if (Tools.ArraysEqual(array1, array2))
                return 1.0;

            int len = Math.Min(array1.Length, array2.Length) - 1;
            Matrix C = CovarianceMatrixSerialCorrelationCorrected(array1, array2, decayFactor, len);
            return C[0, 1] / Math.Sqrt(C[0, 0] * C[1, 1]);
        }

        internal static Matrix CovarianceMatrixSerialCorrelationCorrected(double[] array1, double[] array2, double decayFactor, int length)
        {
            double[] array1Len = Tools.MostRecentValues(array1, length);
            double[] array2Len = Tools.MostRecentValues(array2, length);
            double[] array1Lag = Tools.LagArray(array1, length, 1);
            double[] array2Lag = Tools.LagArray(array2, length, 1);

            double[,] lagArray = new double[length, 2];
            for (int i = 0; i < length; i++)
            {
                lagArray[i, 0] = array1Lag[i];
                lagArray[i, 1] = array2Lag[i];
            }
            WeightedLeastSquares wls1 = new WeightedLeastSquares(array1Len, lagArray, true, decayFactor);
            WeightedLeastSquares wls2 = new WeightedLeastSquares(array2Len, lagArray, true, decayFactor);
            Matrix D = new Matrix(2, 2);
            D[0, 0] = wls1.Betas[1].Value;
            D[0, 1] = wls1.Betas[2].Value;
            D[1, 0] = wls2.Betas[1].Value;
            D[1, 1] = wls2.Betas[2].Value;

            Matrix C = new Matrix(2, 2);
            C[0, 0] = Variance(array1Len, decayFactor);
            C[0, 1] = Covariance(array1Len, array2Len, decayFactor);
            C[1, 0] = C[0, 1];
            C[1, 1] = Variance(array2Len, decayFactor);

            Matrix M = SerialCorrelationCorrection.SerialCorrelationVarianceCorrection(D, C);
            Matrix CStar = M.Multiply(C);
            return CStar;
        }
    }

    public partial class ValueAtRisk
    {
        /// <summary>
        /// Hybrid VaR is based on historical data, but weights more recent data more heavily.
        /// </summary>
        /// <param name="returnArray">Historical returns from which VaR is to be calculated. The last value, with the highest index is assumed to be the most recent data point.</param>
        /// <param name="windowLength">Length of the VaR window. The number of historical returns that will be used to calculate the VaR.</param>
        /// <param name="confidenceLevel">VaR confidence level. 95% and 99% are typical values. If confidenceLevel is 95% then we expect 95% of days to be better (have a more positive return) than the VaR.</param>
        /// <param name="decayFactor">For a decay factor d, the most recent data point has weight 1, the second d, the third d^2, ...</param>
        /// <param name="correctSerialCorrelation">If true, correction is for infinite return length.</param>
        public static double HybridValueAtRisk(double[] returnArray, int windowLength, double confidenceLevel, double decayFactor, bool correctSerialCorrelation)
        {
            double uncorrected = HybridValueAtRisk(returnArray, windowLength, confidenceLevel, decayFactor);
            if (!correctSerialCorrelation)
                return uncorrected;

            double[] subArray = Tools.MostRecentValues(returnArray, windowLength);
            double serialCorrelation = Moments.SerialCorrelation(subArray, decayFactor);
            double varianceCorrectionFactor = SerialCorrelationCorrection.SerialCorrelationVarianceCorrection(serialCorrelation); //don't use version w/ windowLength. windowLenght != length. Length is how far back we go to get data for the calculation. windowLenght is the return lenght that we are correcting to (infinite, by default).
            return uncorrected * Math.Sqrt(varianceCorrectionFactor);
        }


        /// <summary>
        /// Returns the incremental incremental VaR for the subportfolio relative to the portfolio.
        /// </summary>
        /// <param name="portfolioArray">Return array of portfolio. The last element (index = n-1), is the most recent.</param>
        /// <param name="subPortfolioArray">Return array of the sub-portfolio for which incremental VaR is being measured. The last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1. Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="portfolioValutAtRisk">Value at risk of the portfolio.</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        /// <param name="correctSerialCorrelation">If true, correction is for infinite return length.</param>
        public static double IncrementalValueAtRisk(double[] portfolioArray, double[] subPortfolioArray, double decayFactor, double portfolioValutAtRisk, int length, bool correctSerialCorrelation)
        {
            double uncorrected = IncrementalValueAtRisk(portfolioArray, subPortfolioArray, decayFactor, portfolioValutAtRisk, length);
            if (!correctSerialCorrelation || uncorrected == 0.0)
                return uncorrected;

            int len = portfolioArray.Length;
            double[] otherArray = new double[len];
            for(int i = 0; i < len; i++)
                otherArray[i] = portfolioArray[i] - subPortfolioArray[i];

            Matrix CStar = Moments.CovarianceMatrixSerialCorrelationCorrected(portfolioArray, otherArray, decayFactor, length);
            double sdPort = Math.Sqrt(CStar[0, 0] + CStar[1, 1] + 2 * CStar[0, 1]);
            double iStandardDeviation = (CStar[0, 0] + CStar[0, 1]) / sdPort;

            return (portfolioValutAtRisk / sdPort) * iStandardDeviation;
        }
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.