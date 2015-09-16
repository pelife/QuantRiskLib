using System;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///
    ///References:
    ///Miller, Michael B. 2012. Mathematics and Statistics for Financial Risk Management. New York: John Wiley & Sons.
    ///Chapters 8 and 10

    public class Regression
    {
        #region Public Assessors
        public Matrix Y { get; private set; }
        public Matrix X { get; private set; }

        /// <summary>
        /// Array of parameter estimates.
        /// </summary>
        public ParameterEstimate[] Betas { get; private set; }

        /// <summary>
        /// Total sum of squares.
        /// </summary>
        public double TSS { get { return ESS + RSS; } }

        /// <summary>
        /// Explained sum of squares.
        /// </summary>
        public double ESS { get; private set; }

        /// <summary>
        /// Residual sum of squares.
        /// </summary>
        public double RSS { get; private set; }

        /// <summary>
        /// R-squared is also known as the coefficient of determination.
        /// </summary>
        public double rSquared { get; private set; }

        /// <summary>
        /// Adjusted R-squared.
        /// </summary>
        public double AdjustedRSquared { get; private set; }

        /// <summary>
        /// F-statistic.
        /// </summary>
        public double fStatistic { get; private set; }

        /// <summary>
        /// Probability that the F-statistic is not significant. Values closer to 0% indicate that the regression parameters are more statistically significant.
        /// </summary>
        public double fValue { get; private set; }

        /// <summary>
        /// The standard deviation of the residuals.
        /// TrackingError = (RSS/(n-1))^0.5
        /// </summary>
        public double TrackingError { get; private set; }

        /// <summary>
        /// Residuals = Y - E[Y]
        /// </summary>
        public Matrix Residuals { get; private set; }
        #endregion

        #region Constuctors
        /// <summary>
        /// Populates the data and calculates statistics for an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand.</param>
        /// <param name="xData">Data for the independent variable, or regressor.</param>
        /// <param name="includeConstant">If true a column of ones is added to the start of the xData array.</param>
        /// <param name="calculateStatistics">If true, t-statistics, R-squared, F-statistics and related values are calculated.</param>
        public Regression(double[] yData, double[] xData, bool includeConstant, bool calculateStatistics = false)
        {
            double[,] xDataTwoD = new double[xData.Length, 1];
            for (int i = 0; i < xData.Length; i++)
                xDataTwoD[i, 0] = xData[i];

            PopulateMatrices(yData, xDataTwoD, includeConstant);
            CalculateParametersAndStatistic(calculateStatistics);
        }

        /// <summary>
        /// Populates the data and calculates statistics for an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand. Vector is Tx1, where T is the number of time periods.</param>
        /// <param name="xData">Data for the independent variables, or regressors. Matrix is TxN, where T is the number of time periods, and N is the number of independent variables.</param>
        /// <param name="includeConstant">If true a column of ones is added to the start of the xData array.</param>
        /// <param name="calculateStatistics">If true, t-statistics, R-squared, F-statistics and related values are calculated.</param>
        public Regression(double[] yData, double[,] xData, bool includeConstant, bool calculateStatistics = false)
        {
            PopulateMatrices(yData, xData, includeConstant);
            CalculateParametersAndStatistic(calculateStatistics);
        }

        /// <summary>
        /// Populates the data and calculates statistics for an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand.</param>
        /// <param name="xData">Data for the independent variable, or regressor.</param>
        /// <param name="includeConstant">If true a column of ones is added to the start of the xData array.</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        /// <param name="calculateStatistics">If true, t-statistics, R-squared, F-statistics and related values are calculated.</param>
        public Regression(double[] yData, double[] xData, bool includeConstant, int length, bool calculateStatistics = false)
        {
            double[] yDataSubSet = Tools.MostRecentValues(yData, length);
            double[] xDataSubSet = Tools.MostRecentValues(xData, length);

            double[,] xDataTwoDSubSet = new double[length, 1];
            for (int i = 0; i < length; i++)
                xDataTwoDSubSet[i, 0] = xDataSubSet[i];

            PopulateMatrices(yDataSubSet, xDataTwoDSubSet, includeConstant);
            CalculateParametersAndStatistic(calculateStatistics);
        }

        /// <summary>
        /// Populates the data and calculates statistics for an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand. Vector is Tx1, where T is the number of time periods.</param>
        /// <param name="xData">Data for the independent variables, or regressors. Matrix is TxN, where T is the number of time periods, and N is the number of independent variables.</param>
        /// <param name="includeConstant">If true a column of ones is added to the start of the xData array.</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        /// <param name="calculateStatistics">If true, t-statistics, R-squared, F-statistics and related values are calculated.</param>
        public Regression(double[] yData, double[,] xData, bool includeConstant, int length, bool calculateStatistics = false)
        {
            double[] yDataSubSet = Tools.MostRecentValues(yData, length);
            double[,] xDataSubSet = Tools.MostRecentValues(xData, length);

            PopulateMatrices(yDataSubSet, xDataSubSet, includeConstant);
            CalculateParametersAndStatistic(calculateStatistics);
        }
        #endregion

        private void PopulateMatrices(double[] yData, double[,] xData, bool includeConstant)
        {
            Y = new Matrix(yData);
            if (includeConstant)
            {
                double[,] xDataPlus = new double[xData.GetLength(0), xData.GetLength(1) + 1];
                for (int r = 0; r < xDataPlus.GetLength(0); r++)
                {
                    xDataPlus[r, 0] = 1.0;
                    for (int c = 1; c < xDataPlus.GetLength(1); c++)
                        xDataPlus[r, c] = xData[r, c - 1];
                }
                X = new Matrix(xDataPlus);
            }
            else
            {
                X = new Matrix(xData);
            }
            if (Y.NRows != X.NRows)
                throw new ArgumentException("Dependent and independent data do not have the same number of rows.");

        }

        private void CalculateParametersAndStatistic(bool calculateStatistics)
        {
            Matrix XprimeXinv = X.Transpose().Multiply(X).PseudoInverse();
            Matrix betaValues = XprimeXinv.Multiply(X.Transpose()).Multiply(Y);

            int t = Y.NRows;
            int n = betaValues.NRows;
            if (calculateStatistics == false)
            {
                Betas = new ParameterEstimate[n];
                for (int i = 0; i < n; i++)
                    Betas[i] = new ParameterEstimate(betaValues[i, 0]);
                return;
            }

            Matrix expectedValueY = X.Multiply(betaValues);
            Matrix Z = expectedValueY.DeviationsFromMean();
            ESS = Z.Transpose().Multiply(Z)[0, 0];
            Residuals = expectedValueY.Subtract(Y);
            RSS = Residuals.Transpose().Multiply(Residuals)[0, 0];
            TrackingError = Math.Sqrt(RSS / (t - 1));
            rSquared = 1 - RSS / TSS;
            AdjustedRSquared = 1 - (1 - rSquared) * (t - 1) / (t - n);

            
            double varianceOfErrorTerm = Residuals.Transpose().Multiply(Residuals)[0, 0] / (t - n);
            Betas = new ParameterEstimate[n];
            for (int i = 0; i < n; i++)
            {
                double sigma = Math.Sqrt(varianceOfErrorTerm * XprimeXinv[i, i]);
                double tStat = betaValues[i, 0] == 0 ? 0 : Math.Abs(betaValues[i, 0] - 0) / sigma;
                double pValue = 1.0 - Distributions.StudentsTCumulativeDensityFunction(tStat, t - n);
                Betas[i] = new ParameterEstimate(betaValues[i, 0], sigma, tStat, pValue);
            }

            fStatistic = (rSquared / (n - 1)) / ((1 - rSquared) / (t - n));

            if(n > 1)
                fValue = 1.0 - Distributions.FCumulativeDensityFunction(fStatistic, n - 1, t - n);
        }
    }

    public class WeightedLeastSquares
    {
        #region Private fields
        private readonly double _decayFactor;
        private Matrix _Xstar;
        private Matrix _Ystar;
        #endregion

        #region Public Assessors
        public Matrix Y { get; private set; }
        public Matrix X { get; private set; }

        /// <summary>
        /// Array of parameter estimates.
        /// </summary>
        public ParameterEstimate[] Betas { get; private set; }

        /// <summary>
        /// Total sum of squares.
        /// </summary>
        public double TSS { get { return ESS + RSS; } }

        /// <summary>
        /// Explained sum of squares.
        /// </summary>
        public double ESS { get; private set; }

        /// <summary>
        /// Residual sum of squares.
        /// </summary>
        public double RSS { get; private set; }

        /// <summary>
        /// R-squared is also known as the coefficient of determination.
        /// </summary>
        public double rSquared { get; private set; }

        /// <summary>
        /// Adjusted R-squared.
        /// </summary>
        public double AdjustedRSquared { get; private set; }

        /// <summary>
        /// F-statistic.
        /// </summary>
        public double fStatistic { get; private set; }

        /// <summary>
        /// Probability that the F-statistic is not significant. Values closer to 0% indicate that the regression parameters are more statistically significant.
        /// </summary>
        public double fValue { get; private set; }

        /// <summary>
        /// The standard deviation of the residuals.
        /// TrackingError = (RSS/(n-1))^0.5
        /// </summary>
        public double TrackingError { get; private set; }

        /// <summary>
        /// Residuals = Y - E[Y]
        /// </summary>
        public Matrix Residuals { get; private set; }
        #endregion

        #region Constuctors
        /// <summary>
        /// Populates the data and calculates statistics for an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand.</param>
        /// <param name="xData">Data for the independent variable, or regressor.</param>
        /// <param name="includeConstant">If true a column of ones is added to the start of the xData array.</param>
        /// <param name="decayFactor">More recent data (higher index value) is weighted more heavily. Most recent data point has weight d, then d^2, then d^3, ... . If decay factor is 1.0 then weighted least squares is equivalent to ordinary least squares.</param>
        /// <param name="calculateStatistics">If true, t-statistics, R-squared, F-statistics and related values are calculated.</param>
        public WeightedLeastSquares(double[] yData, double[] xData, bool includeConstant, double decayFactor, bool calculateStatistics = false)
        {
            _decayFactor = decayFactor;

            double[,] xDataTwoD = new double[xData.Length, 1];
            for (int i = 0; i < xData.Length; i++)
                xDataTwoD[i, 0] = xData[i];

            PopulateMatrices(yData, xDataTwoD, includeConstant);
            CalculateParametersAndStatistics(calculateStatistics);
        }

        /// <summary>
        /// Populates the data and calculates statistics for an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand. Vector is Tx1, where T is the number of time periods.</param>
        /// <param name="xData">Data for the independent variables, or regressors. Matrix is TxN, where T is the number of time periods, and N is the number of independent variables.</param>
        /// <param name="includeConstant">If true a column of ones is added to the start of the xData array.</param>
        /// <param name="decayFactor">More recent data (higher index value) is weighted more heavily. Most recent data point has weight d, then d^2, then d^3, ... . If decay factor is 1.0 then weighted least squares is equivalent to ordinary least squares.</param>
        /// <param name="calculateStatistics">If true, t-statistics, R-squared, F-statistics and related values are calculated.</param>
        public WeightedLeastSquares(double[] yData, double[,] xData, bool includeConstant, double decayFactor, bool calculateStatistics = false)
        {
            _decayFactor = decayFactor;

            PopulateMatrices(yData, xData, includeConstant);
            CalculateParametersAndStatistics(calculateStatistics);
        }

        /// <summary>
        /// Populates the data and calculates statistics for an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand.</param>
        /// <param name="xData">Data for the independent variable, or regressor.</param>
        /// <param name="includeConstant">If true a column of ones is added to the start of the xData array.</param>
        /// <param name="decayFactor">More recent data (higher index value) is weighted more heavily. Most recent data point has weight d, then d^2, then d^3, ... . If decay factor is 1.0 then weighted least squares is equivalent to ordinary least squares.</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        /// <param name="calculateStatistics">If true, t-statistics, R-squared, F-statistics and related values are calculated.</param>
        public WeightedLeastSquares(double[] yData, double[] xData, bool includeConstant, double decayFactor, int length, bool calculateStatistics = false)
        {
            _decayFactor = decayFactor;

            double[] yDataSubSet = Tools.MostRecentValues(yData, length);
            double[] xDataSubSet = Tools.MostRecentValues(xData, length);

            double[,] xDataTwoDSubSet = new double[length, 1];
            for (int i = 0; i < length; i++)
                xDataTwoDSubSet[i, 0] = xDataSubSet[i];

            PopulateMatrices(yDataSubSet, xDataTwoDSubSet, includeConstant);
            CalculateParametersAndStatistics(calculateStatistics);
        }

        /// <summary>
        /// Populates the data and calculates statistics for an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand. Vector is Tx1, where T is the number of time periods.</param>
        /// <param name="xData">Data for the independent variables, or regressors. Matrix is TxN, where T is the number of time periods, and N is the number of independent variables.</param>
        /// <param name="includeConstant">If true a column of ones is added to the start of the xData array.</param>
        /// <param name="decayFactor">More recent data (higher index value) is weighted more heavily. Most recent data point has weight d, then d^2, then d^3, ... . If decay factor is 1.0 then weighted least squares is equivalent to ordinary least squares.</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        /// <param name="calculateStatistics">If true, t-statistics, R-squared, F-statistics and related values are calculated.</param>
        public WeightedLeastSquares(double[] yData, double[,] xData, bool includeConstant, double decayFactor, int length, bool calculateStatistics = false)
        {
            _decayFactor = decayFactor;

            double[] yDataSubSet = Tools.MostRecentValues(yData, length);
            double[,] xDataSubSet = Tools.MostRecentValues(xData, length);

            PopulateMatrices(yDataSubSet, xDataSubSet, includeConstant);
            CalculateParametersAndStatistics(calculateStatistics);

        }
        #endregion

        private void PopulateMatrices(double[] yData, double[,] xData, bool includeConstant)
        {
            Y = new Matrix(yData);
            if (includeConstant)
            {
                double[,] xDataPlus = new double[xData.GetLength(0), xData.GetLength(1) + 1];
                for (int r = 0; r < xDataPlus.GetLength(0); r++)
                {
                    xDataPlus[r, 0] = 1.0;
                    for (int c = 1; c < xDataPlus.GetLength(1); c++)
                        xDataPlus[r, c] = xData[r, c - 1];
                }
                X = new Matrix(xDataPlus);
            }
            else
            {
                X = new Matrix(xData);
            }
            if (Y.NRows != X.NRows)
                throw new ArgumentException("Dependent and independent data do not have the same number of rows.");

            //The following is equivalent to the matrix operation _Xstar = W.Multiply(X) and _Ystar = W.Multiply(Y);
            //where W is a diagonal matrix with exponentially declining weights along the diagonal.
            _Xstar = new Matrix(X.NRows, X.NColumns);
            _Ystar = new Matrix(Y.NRows, Y.NColumns);
            double d = 1.0;
            double sqrtDecayFactor = Math.Sqrt(_decayFactor);   //We use the square root of the decay factor to be consistent with our definition of variance.
            for (int r = Y.NRows - 1; r >= 0; r--)
            {
                _Ystar[r, 0] = d * Y[r, 0];
                for (int c = 0; c < X.NColumns; c++)
                {
                    _Xstar[r, c] = d * X[r, c];
                }
                d *= sqrtDecayFactor;
            }

        }

        private void CalculateParametersAndStatistics(bool calculateStatistics)
        {
            Matrix XprimeXinv = _Xstar.Transpose().Multiply(_Xstar).PseudoInverse();
            Matrix betaValues = XprimeXinv.Multiply(_Xstar.Transpose()).Multiply(_Ystar);

            int n = betaValues.NRows;
            if (calculateStatistics == false)
            {
                Betas = new ParameterEstimate[n];
                for (int i = 0; i < n; i++)
                    Betas[i] = new ParameterEstimate(betaValues[i, 0]);
                return;
            }

            int t = Y.NRows;
            Matrix expectedValueY = _Xstar.Multiply(betaValues);
            Matrix Z = expectedValueY.DeviationsFromMean();
            ESS = Z.Transpose().Multiply(Z)[0, 0];
            Residuals = expectedValueY.Subtract(_Ystar);
            RSS = Residuals.Transpose().Multiply(Residuals)[0, 0];
            TrackingError = Math.Sqrt(RSS / (t - 1));
            rSquared = ESS / TSS;
            AdjustedRSquared = 1 - (1 - rSquared) * (t - 1) / (t - n);

            double varianceOfErrorTerm = Residuals.Transpose().Multiply(Residuals)[0, 0] / (t - n);
            Betas = new ParameterEstimate[n];
            for (int i = 0; i < n; i++)
            {
                double sigma = Math.Sqrt(varianceOfErrorTerm * XprimeXinv[i, i]);
                double tStat = betaValues[i, 0] == 0 ? 0 : Math.Abs(betaValues[i, 0] - 0) / sigma;
                double pValue = 1.0 - Distributions.StudentsTCumulativeDensityFunction(tStat, t - n);
                Betas[i] = new ParameterEstimate(betaValues[i, 0], sigma, tStat, pValue);
            }

            fStatistic = (rSquared / (n - 1)) / ((1 - rSquared) / (t - n));
            fValue = 1.0 - Distributions.FCumulativeDensityFunction(fStatistic, n - 1, t - n);
        }
    }

    public class ParameterEstimate
    {
        public double Value { get; private set; }
        public double StandardDeviation { get; private set; }

        /// <summary>
        /// This is an indication of how statistically significant the parameter value is. It is tested against the null hypothesis that the parameter is zero. The value is specified to always be positive. |[beta]|/[standard deviation of beta].
        /// </summary>
        public double tStatistic { get; private set; }
        
        /// <summary>
        /// The probability of observing a parameter value this different than zero by chance. Always less than 50%. Closer to 0% is more statistically significant.
        /// </summary>
        public double pValue { get; private set; }

        internal ParameterEstimate(double value, double standardDeviation, double tStatistic, double pValue)
        {
            Value = value;
            StandardDeviation = standardDeviation;
            this.tStatistic = tStatistic;
            this.pValue = pValue;
        }

        internal ParameterEstimate(double value)
        {
            Value = value;
            StandardDeviation = double.NaN;
            tStatistic = double.NaN;
            pValue = double.NaN;
        }
    }

    /// <summary>
    /// This class will calculate the slope, and only the slope for a univariate regression.
    /// While limited, this calculation is much quicker than a full regression analysis.
    /// </summary>
    public class RegressionParametersOnly
    {
        #region beta only
        /// <summary>
        /// Returns the slope of an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand.</param>
        /// <param name="xData">Data for the independent variable, or regressor.</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        public static double Beta(double[] yData, double[] xData, int length)
        {
            double[] xDataSubSet = Tools.MostRecentValues(xData, length);
            double[] yDataSubSet = Tools.MostRecentValues(yData, length);

            double mX = Moments.Mean(xData);
            double mY = Moments.Mean(yData);
            double sumXY = 0, sumXX = 0;
            for (int i = 0; i < length; i++)
            {
                sumXY += xDataSubSet[i] * yDataSubSet[i];
                sumXX += xDataSubSet[i] * xDataSubSet[i];
            }

            return (sumXY - length * mY * mX) / (sumXX - length * mX * mX);
        }

        /// <summary>
        /// Returns the slope of an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand.</param>
        /// <param name="xData">Data for the independent variable, or regressor.</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        /// <param name="decayFactor">More recent data (higher index value) is weighted more heavily. Most recent data point has weight d, then d^2, then d^3, ... . If decay factor is 1.0 then weighted least squares is equivalent to ordinary least squares.</param>
        public static double Beta(double[] yData, double[] xData, double decayFactor, int length)
        {
            double[] xDataSubSet = Tools.MostRecentValues(xData, length);
            double[] yDataSubSet = Tools.MostRecentValues(yData, length);

            return Beta(yDataSubSet, xDataSubSet, decayFactor);
        }

        /// <summary>
        /// Returns the slope of an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand.</param>
        /// <param name="xData">Data for the independent variable, or regressor.</param>
        /// <param name="decayFactor">More recent data (higher index value) is weighted more heavily. Most recent data point has weight d, then d^2, then d^3, ... . If decay factor is 1.0 then weighted least squares is equivalent to ordinary least squares.</param>
        public static double Beta(double[] yData, double[] xData, double decayFactor)
        {
            if (yData.Length != xData.Length)
                throw new ArgumentException("y array and x array must be the same length");

            int n = yData.Length;

            double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumD = 0;
            double d = 1.0;
            for (int i = 0; i < n; i++)
            {
                sumX += d * d * xData[n - 1 - i];
                sumY += d * d * xData[n - 1 - i];
                sumXX += (d * xData[n - 1 - i]) * (d * xData[n - 1 - i]);
                sumXY += (d * xData[n - 1 - i]) * (d * yData[n - 1 - i]);
                sumD += d * d;
                d *= decayFactor;
            }

            return (sumD * sumXY - sumY * sumX) / (sumD * sumXX - sumX * sumX);
        }
        #endregion

        #region alpha and beta

        /// <summary>
        /// Returns the intercept and slope of an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand.</param>
        /// <param name="xData">Data for the independent variable, or regressor.</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        /// <param name="alpha"></param>
        public static double Beta(double[] yData, double[] xData, int length, out double alpha)
        {
            double[] xDataSubSet = Tools.MostRecentValues(xData, length);
            double[] yDataSubSet = Tools.MostRecentValues(yData, length);

            double mX = Moments.Mean(xData);
            double mY = Moments.Mean(yData);
            double sumXY = 0, sumXX = 0;
            for (int i = 0; i < length; i++)
            {
                sumXY += xDataSubSet[i] * yDataSubSet[i];
                sumXX += xDataSubSet[i] * xDataSubSet[i];
            }

            double beta = (sumXY - length * mY * mX) / (sumXX - length * mX * mX);
            alpha = mY - beta * mX;
            return beta;
        }

        /// <summary>
        /// Returns the intercept and slope of an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand.</param>
        /// <param name="xData">Data for the independent variable, or regressor.</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        /// <param name="decayFactor">More recent data (higher index value) is weighted more heavily. Most recent data point has weight d, then d^2, then d^3, ... . If decay factor is 1.0 then weighted least squares is equivalent to ordinary least squares.</param>
        /// <param name="alpha"></param>
        public static double Beta(double[] yData, double[] xData, double decayFactor, int length, out double alpha)
        {
            double[] xDataSubSet = Tools.MostRecentValues(xData, length);
            double[] yDataSubSet = Tools.MostRecentValues(yData, length);

            return Beta(yDataSubSet, xDataSubSet, decayFactor, out alpha);
        }

        /// <summary>
        /// Returns the intercept and slope of an OLS regression analysis.
        /// </summary>
        /// <param name="yData">Data for the dependent variable, or regressand.</param>
        /// <param name="xData">Data for the independent variable, or regressor.</param>
        /// <param name="decayFactor">More recent data (higher index value) is weighted more heavily. Most recent data point has weight d, then d^2, then d^3, ... . If decay factor is 1.0 then weighted least squares is equivalent to ordinary least squares.</param>
        /// <param name="alpha"></param>
        public static double Beta(double[] yData, double[] xData, double decayFactor, out double alpha)
        {
            if (yData.Length != xData.Length)
                throw new ArgumentException("y array and x array must be the same length");

            int n = yData.Length;

            double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0, sumD = 0;
            double d = 1.0;
            for (int i = 0; i < n; i++)
            {
                sumX += d * d * xData[n - 1 - i];
                sumY += d * d * xData[n - 1 - i];
                sumXX += (d * xData[n - 1 - i]) * (d * xData[n - 1 - i]);
                sumXY += (d * xData[n - 1 - i]) * (d * yData[n - 1 - i]);
                sumD += d * d;
                d *= decayFactor;
            }

            alpha = (sumXX * sumY - sumX * sumXY) / (sumD * sumXX - sumX * sumX);
            return (sumD * sumXY - sumY * sumX) / (sumD * sumXX - sumX * sumX);
        }
        #endregion
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.