using System;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///
    ///References:
    ///Miller, Michael B. 2012. Mathematics and Statistics for Financial Risk Management. New York: John Wiley & Sons.
    ///Chapters 2, 3 and 10

    public partial class Moments
    {
        #region Mean
        /// <summary>
        /// Returns the mean of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the mean.</param>
        public static double Mean(double[] array)
        {
            return Mean(array, 1.0);
        }

        /// <summary>
        /// Returns the mean of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the mean.  For time series, the last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in array is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        public static double Mean(double[] array, double decayFactor)
        {
            int n = array.Length;
            double m = 0.0;
            double d = 1.0;
            for (int i = 0; i < n; i++)
            {
                m += d * array[n - 1 - i];
                d *= decayFactor;
            }

            m *= GeometricSeries.InverseSumOfGeometricSeries(decayFactor, n);
            return m;
        }

        /// <summary>
        /// Returns the mean of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the mean.  For time series, the last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in array is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        public static double Mean(double[] array, double decayFactor, int length)
        {
            double[] subArray = Tools.MostRecentValues(array, length);
            return Mean(subArray, decayFactor);
        }

        /// <summary>
        /// Returns the weighted averages of the values in valueArray using the corresponding weights in weightArray.
        /// </summary>
        /// <param name="valueArray">array of values for which we are computing the weighted average</param>
        /// <param name="weightArray">array of weights used in computing the weighted average</param>
        public static double WeightedMean(double[] valueArray, double[] weightArray)
        {
            if (valueArray.Length != weightArray.Length)
                throw new ArgumentException("valueArray and weightArray must have the same number of elements");
            double m = 0, s = 0;
            for (int i = 0; i < valueArray.Length; i++)
            {
                m += valueArray[i] * weightArray[i];
                s += weightArray[i];
            }
            m /= s;

            return m;
        }
        #endregion

        #region StandardDeviation
        /// <summary>
        /// Returns the sample standard deviation of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the standard deviation.</param>
        public static double StandardDeviation(double[] array)
        {
            return StandardDeviation(array, 1.0);
        }

        /// <summary>
        /// Returns the sample standard deviation of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the standard deviation.  For time series, the last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in array is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        public static double StandardDeviation(double[] array, double decayFactor)
        {
            double v = Variance(array, decayFactor);
            return Math.Sqrt(v);
        }

        /// <summary>
        /// Returns the sample standard deviation of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the standard deviation.  For time series, the last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in array is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        public static double StandardDeviation(double[] array, double decayFactor, int length)
        {
            double[] subArray = Tools.MostRecentValues(array, length);
            return StandardDeviation(subArray, decayFactor);
        }

        /// <summary>
        /// Calculates the sample standard deviation of an array using a rolling window.
        /// The first element (index = 0) of the output array is the standard deviation of points [0] to [length - 1] in 'array', 
        /// the second element of the output array is the standard deviation of points [1] to [length] in 'array', and so on.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the standard deviation.  For time series, the last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in array is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="length">>Length of rolling window length. Must be less than or equal to the lenght of the array of data.</param>
        /// <returns></returns>
        public static double[] RollingStandardDeviation(double[] array, double decayFactor, int length)
        {
            double[] rollingSdArray = RollingVariance(array, decayFactor, length);
            for (int i = 0; i < rollingSdArray.Length; i++)
                rollingSdArray[i] = Math.Sqrt(rollingSdArray[i]);
            return rollingSdArray;
        }
        #endregion

        #region Incremental Standard Deviation
        /// <summary>
        /// Returns the incremental standard deviation of a variable. Incremental standard deviation is to standard deviation what incremental VaR is to VaR.
        /// </summary>
        /// <param name="portfolioArray">Return array of portfolio. The last element (index = n-1), is the most recent.</param>
        /// <param name="subPortfolioArray">Return array of the sub-portfolio for which incremental standard deviation is being measured. The last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1. Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        public static double IncrementalStandardDeviation(double[] portfolioArray, double[] subPortfolioArray, double decayFactor)
        {
            if (portfolioArray.Length != subPortfolioArray.Length)
                throw new ArgumentException("Arrays must be the same length");

            int len = portfolioArray.Length;
            double[] otherArray = new double[len];
            for (int i = 0; i < len; i++)
                otherArray[i] = portfolioArray[i] - subPortfolioArray[i];

            double varianceSubPort = Variance(subPortfolioArray, decayFactor);
            double covarSubOther = Covariance(subPortfolioArray, otherArray, decayFactor);
            //double sdPort = StandardDeviation(portfolioArray, decayFactor);
            double varianceOther = Variance(otherArray, decayFactor);
            double sdPort = Math.Sqrt(varianceSubPort + varianceOther + 2 * covarSubOther);
            double iStandardDeviation = (varianceSubPort + covarSubOther) / sdPort;
            return iStandardDeviation;
        }

        /// <summary>
        /// Returns the incremental standard deviation of a variable. Incremental standard deviation is to standard deviation what incremental VaR is to VaR.
        /// </summary>
        /// <param name="portfolioArray">Return array of portfolio. The last element (index = n-1), is the most recent.</param>
        /// <param name="subPortfolioArray">Return array of the sub-portfolio for which incremental standard deviation is being measured. The last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1. Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        public static double IncrementalStandardDeviation(double[] portfolioArray, double[] subPortfolioArray, double decayFactor, int length)
        {
            double[] portfolioArrayLen = Tools.MostRecentValues(portfolioArray, length);
            double[] subPortfolioArrayLen = Tools.MostRecentValues(subPortfolioArray, length);
            return IncrementalStandardDeviation(portfolioArrayLen, subPortfolioArrayLen, decayFactor);
        }
        #endregion

        #region Variance
        /// <summary>
        /// Returns the sample variance of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the variance.</param>
        public static double Variance(double[] array)
        {
            return Variance(array, 1.0);
        }

        /// <summary>
        /// Returns the sample variance of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the variance.  For time series, the last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in array is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        public static double Variance(double[] array, double decayFactor)
        {
            if(array.Length == 1)
                return 0.0;

            int n = array.Length;
            double v = 0.0;
            double d = 1.0;
            for (int i = 0; i < n; i++)
            {
                v += d * array[n - 1 - i] * array[n - 1 - i];
                d *= decayFactor;
            }

            double m = Mean(array, decayFactor);
            double s1 = GeometricSeries.SumOfGeometricSeries(decayFactor, n);
            double s2 = GeometricSeries.SumOfGeometricSeries(decayFactor * decayFactor, n);
            double a = s1 / (s1 * s1 - s2);
            double b = s1 * a;
            v = a * v - b * m * m;

            //Occasionally floating-point error will cause v to be less than zero when it should be exactly zero or very close to zero.
            if(v < 0)
            {
                bool allEqual = Tools.ArrayAllEqual(array);
                return allEqual ? 0 : double.Epsilon;
            }

            return v;
        }

        /// <summary>
        /// Returns the sample variance of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the variance.  For time series, the last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in array is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        /// <returns></returns>
        public static double Variance(double[] array, double decayFactor, int length)
        {
            double[] subArray = Tools.MostRecentValues(array, length);
            return Variance(subArray, decayFactor);
        }

        /// <summary>
        /// Calculates the sample variance of an array using a rolling window.
        /// The first element (index = 0) of the output array is the variance of points [0] to [length - 1] in 'array', 
        /// the second element of the output array is the variance of points [1] to [length] in 'array', and so on.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the variance.  For time series, the last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in array is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="length">>Length of rolling window length. Must be less than or equal to the lenght of the array of data.</param>
        /// <returns></returns>
        public static double[] RollingVariance(double[] array, double decayFactor, int length)
        {
            if (array.Length < length)
                throw new ArgumentException("Not enought points in array. Number of points in array must be greater than or equal to length.");

            double s1 = GeometricSeries.SumOfGeometricSeries(decayFactor, length);
            double s2 = GeometricSeries.SumOfGeometricSeries(decayFactor * decayFactor, length);
            double c = 1 / (s1 * s1 - s2);
            double a = s1 * c;
            
            double z1 = 0.0, z2 = 0.0;
            double d = 1.0;
            for (int i = 0; i < length; i++)
            {
                z1 += d * array[length - 1 - i] * array[length - 1 - i];
                z2 += d * array[length - 1 - i];
                d *= decayFactor;
            }

            double dToLength = Math.Pow(decayFactor, length);
            double[] vArray = new double[array.Length - length + 1];
            vArray[0] = a * z1 - c * z2 * z2;
            for(int i = 1; i < vArray.Length; i++)
            {
                z1 = decayFactor * z1 + array[length - 1 + i] * array[length - 1 + i] - dToLength * array[i - 1] * array[i - 1];
                z2 = decayFactor * z2 + array[length - 1 + i] - dToLength * array[i - 1];
                vArray[i] = a * z1 - c * z2 * z2;
            }
            return vArray;
        }
        #endregion

        #region Covariance and Correlation
        /// <summary>
        /// Returns the sample covariance between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// </summary>
        /// <param name="array1"></param>
        /// <param name="array2"></param>
        public static double Covariance(double[] array1, double[] array2)
        {
            if (array1.Length != array2.Length)
                throw new ArgumentException("Arrays must be the same length");

            int n = array1.Length;
            if (n == 1) return double.NaN;
            double m1 = Mean(array1);
            double m2 = Mean(array2);
            
            
            double covar = 0.0;
            for (int i = 0; i < n; i++)
                covar += (array1[i] - m1) * (array2[i] - m2);
            covar /= (n - 1);
            return covar;
        }

        /// <summary>
        /// Returns the sample covariance between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// </summary>
        /// <param name="array1"></param>
        /// <param name="array2"></param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        public static double Covariance(double[] array1, double[] array2, double decayFactor)
        {
            if (array1.Length != array2.Length)
                throw new ArgumentException("Arrays must be the same length");

            if (Tools.ArrayAllEqual(array1) || Tools.ArrayAllEqual(array2))
                return 0.0;

            int n = array1.Length;
            if (n == 1) return double.NaN;
            double covar = 0.0;
            double d = 1.0;
            for (int i = 0; i < n; i++)
            {
                covar += d * array1[n - 1 - i] * array2[n - 1 - i];
                d *= decayFactor;
            }

            double m1 = Mean(array1, decayFactor);
            double m2 = Mean(array2, decayFactor);
            double S1 = GeometricSeries.SumOfGeometricSeries(decayFactor, n);
            double S2 = GeometricSeries.SumOfGeometricSeries(decayFactor * decayFactor, n);
            double A = S1 / (S1 * S1 - S2);
            double B = S1 * A;
            covar = A * covar - B * m1 * m2;
            return covar;
        }

        /// <summary>
        /// Returns the correlation between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// </summary>
        /// <param name="array1"></param>
        /// <param name="array2"></param>
        public static double Correlation(double[] array1, double[] array2)
        {
            double covar = Covariance(array1, array2);
            double sd1 = StandardDeviation(array1);
            double sd2 = StandardDeviation(array2);
            double corr = covar / (sd1 * sd2);

            //In theory we shouldn't get values less than -1 or greater than +1, but with precision errors we can get values slightly outside these bounds.
            if (corr < -1) return -1;
            if (corr > 1) return 1;
            return corr;
        }

        /// <summary>
        /// Returns the correlation between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// </summary>
        /// <param name="array1"></param>
        /// <param name="array2"></param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        public static double Correlation(double[] array1, double[] array2, double decayFactor)
        {
            double covar = Covariance(array1, array2, decayFactor);
            double sd1 = StandardDeviation(array1, decayFactor);
            double sd2 = StandardDeviation(array2, decayFactor);
            double corr = covar / (sd1 * sd2);

            //In theory we shouldn't get values less than -1 or greater than +1, but with precision errors we can get values slightly outside these bounds.
            if (corr < -1) return -1;
            if (corr > 1) return 1;
            return corr;
        }

        /// <summary>
        /// Returns the serial covariance of an array of data.
        /// The serial covariance is the covariance of the data, with its own first lag.
        /// </summary>
        public static double SerialCovariance(double[] array)
        {
            int n = array.Length - 1;
            double[] sub1 = new double[n];
            double[] sub2 = new double[n];
            for (int i = 0; i < n; i++)
            {
                sub1[i] = array[i];
                sub2[i] = array[i + 1];
            }
            return Covariance(sub1, sub2);
        }

        /// <summary>
        /// Returns the serial covariance of an array of data.
        /// The serial covariance is the covariance of the data, with its own first lag.
        /// </summary>
        /// <param name="array"></param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        public static double SerialCovariance(double[] array, double decayFactor)
        {
            int n = array.Length - 1;
            double[] sub1 = new double[n];
            double[] sub2 = new double[n];
            for (int i = 0; i < n; i++)
            {
                sub1[i] = array[i];
                sub2[i] = array[i + 1];
            }
            return Covariance(sub1, sub2, decayFactor);
        }

        /// <summary>
        /// Returns the serial correlation of an array of data.
        /// The serial correlation is the correlation of the data, with its own first lag.
        /// </summary>
        public static double SerialCorrelation(double[] array)
        {
            int n = array.Length - 1;
            double[] sub1 = new double[n];
            double[] sub2 = new double[n];
            for (int i = 0; i < n; i++)
            {
                sub1[i] = array[i];
                sub2[i] = array[i + 1];
            }
            return Correlation(sub1, sub2);
        }

        /// <summary>
        /// Returns the serial correlation of an array of data.
        /// The serial correlation is the correlation of the data, with its own first lag.
        /// </summary>
        /// <param name="array"></param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1.  Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        public static double SerialCorrelation(double[] array, double decayFactor)
        {
            int n = array.Length - 1;
            double[] sub1 = new double[n];
            double[] sub2 = new double[n];
            for (int i = 0; i < n; i++)
            {
                sub1[i] = array[i];
                sub2[i] = array[i + 1];
            }
            return Correlation(sub1, sub2, decayFactor);
        }
        #endregion

        #region Skewness
        /// <summary>
        /// Returns the sample skewness of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the skewness.</param>
        public static double Skewness(double[] array)
        {
            int n = array.Length;
            if(n < 3) return double.NaN;
            double m = Mean(array);
            double sd = StandardDeviation(array);

            double skew = 0.0;
            for (int i = 0; i < n; i++)
                skew += Math.Pow(array[i] - m, 3.0);
            skew /= Math.Pow(sd, 3.0);
            skew = n * skew / ((n - 1.0) * (n - 2.0));
            return skew;
        }
        #endregion

        #region Kurtosis
        /// <summary>
        /// Returns the sample kurtosis of an array.
        /// Note: This is not excess kurtosis.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the kurtosis.</param>
        public static double Kurtosis(double[] array)
        {
            int n = array.Length;
            if (n < 4) return double.NaN;
            double m = Mean(array);
            double variance = Variance(array);

            double kurt = 0.0;
            for (int i = 0; i < n; i++)
                kurt += Math.Pow(array[i] - m, 4.0);
            kurt /= (variance * variance);
            kurt *= (n * (n + 1.0)) / ((n - 1.0) * (n - 2.0) * (n - 3.0));
            return kurt;
        }

        /// <summary>
        /// Returns the sample excess kurtosis of an array.
        /// </summary>
        /// <param name="array">Array of data for which we are calculating the excess kurtosis.</param>
        public static double ExcessKurtosis(double[] array)
        {
            double kurt = Kurtosis(array);
            int n = array.Length;
            kurt -= 3 * (n - 1.0) * (n - 1.0) / ((n - 2.0) * (n - 3.0));
            return kurt;
        }
        #endregion

        #region Coskewness and Cokurtosis
        public enum CoskewnessType
        {
            AAB, ABB
        }

        public enum CokurtosisType
        {
            AAAB, AABB, ABBB
        }

        /// <summary>
        /// Returns the sample coskewness between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// </summary>
        /// <param name="arrayA"></param>
        /// <param name="arrayB"></param>
        /// <param name="coskewnessType"></param>
        public static double Coskewness(double[] arrayA, double[] arrayB, CoskewnessType coskewnessType)
        {
            double m3 = ThirdCentralCrossMoment(arrayA, arrayB, coskewnessType);
            double sd1 = StandardDeviation(arrayA);
            double sd2 = StandardDeviation(arrayB);
            if (coskewnessType == CoskewnessType.AAB)
                return m3/(sd1*sd1*sd2);
            return m3 / (sd1 * sd2 * sd2);
        }

        /// <summary>
        /// Returns the sample cokurtosis between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// </summary>
        /// <param name="arrayA"></param>
        /// <param name="arrayB"></param>
        /// <param name="cokurtosisType"></param>
        public static double Coskurtosis(double[] arrayA, double[] arrayB, CokurtosisType cokurtosisType)
        {
            double m4 = FourthCentralCrossMoment(arrayA, arrayB, cokurtosisType);
            double sd1 = StandardDeviation(arrayA);
            double sd2 = StandardDeviation(arrayB);
            if (cokurtosisType == CokurtosisType.AAAB)
                return m4 / (sd1 * sd1 * sd1 * sd2);
            if(cokurtosisType == CokurtosisType.AABB)
               return m4 / (sd1 * sd1 * sd2 * sd2);
            return m4 / (sd1 * sd2 * sd2 * sd2);
        }

        /// <summary>
        /// Returns the sample third central cross moment between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// </summary>
        /// <param name="arrayA"></param>
        /// <param name="arrayB"></param>
        /// <param name="coskewnessType"></param>
        private static double ThirdCentralCrossMoment(double[] arrayA, double[] arrayB, CoskewnessType coskewnessType)
        {
            if(arrayA.Length != arrayB.Length)
                throw new ArgumentException("arrayA and arrayB should be the same length.");
            int n = arrayA.Length;

            double meanA = Mean(arrayA);
            double meanB = Mean(arrayB);
            double m3 = 0;
            if (coskewnessType == CoskewnessType.AAB)
            {
                for (int i = 0; i < n; i++)
                    m3 += (arrayA[i] - meanA) * (arrayA[i] - meanA) * (arrayB[i] - meanB);
            }
            else
            {
                for (int i = 0; i < n; i++)
                    m3 += (arrayA[i] - meanA) * (arrayB[i] - meanB) * (arrayB[i] - meanB);
            }
            m3 *= n / ((n - 1.0) * (n - 2.0)); //because this is the sample coskew multiply by n/(n-1)(n-2) rather than 1/n
            return m3;
        }

        /// <summary>
        /// Returns the sample fourth central cross moment between two arrays.
        /// Arrays should be of equal length, and contain more than one element.
        /// </summary>
        /// <param name="arrayA"></param>
        /// <param name="arrayB"></param>
        /// <param name="cokurtosisType"></param>
        private static double FourthCentralCrossMoment(double[] arrayA, double[] arrayB, CokurtosisType cokurtosisType)
        {
            if (arrayA.Length != arrayB.Length)
                throw new ArgumentException("arrayA and arrayB should be the same length.");
            int n = arrayA.Length;

            double meanA = Mean(arrayA);
            double meanB = Mean(arrayB);
            double m4 = 0;
            if (cokurtosisType == CokurtosisType.AAAB)
            {
                for (int i = 0; i < n; i++)
                    m4 += (arrayA[i] - meanA) * (arrayA[i] - meanA) * (arrayA[i] - meanA) * (arrayB[i] - meanB);
            }
            else if (cokurtosisType == CokurtosisType.AABB)
            {
                for (int i = 0; i < n; i++)
                    m4 += (arrayA[i] - meanA)*(arrayA[i] - meanA)*(arrayB[i] - meanB)*(arrayB[i] - meanB);
            }
            else
            {
                for (int i = 0; i < n; i++)
                    m4 += (arrayA[i] - meanA) * (arrayB[i] - meanB) * (arrayB[i] - meanB) * (arrayB[i] - meanB);
            }
            m4 *= n * (n - 1) / ((n - 1.0) * (n - 2.0) * (n - 3.0)); //because this is the sample coskew multiply by n(n+1)/(n-1)(n-2)(n-3) rather than 1/n
            return m4;
        }
        #endregion
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.