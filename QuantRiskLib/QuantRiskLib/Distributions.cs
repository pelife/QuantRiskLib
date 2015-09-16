using System;
using System.Collections.Generic;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///
    ///References:
    ///Miller, Michael B. 2012. Mathematics and Statistics for Financial Risk Management. New York: John Wiley & Sons.
    ///Chapter 4

    public partial class Distributions
    {
        #region Beta Distribution

        /// <summary>
        /// Returns the CDF of a beta distribution with shape parameters a and b.
        /// PDF is: f(x) = x^(a-1) * (1-x)^(b-1) / B(a,b)
        /// </summary>
        /// <param name="x">Value of the random variable for which the CDF is beign evaluated. x is between 0 and 1.</param>
        /// <param name="a">First shape parameter.</param>
        /// <param name="b">Second shape parameter.</param>
        /// <returns></returns>
        public static double BetaCumulativeDistributionFunction(double x, double a, double b)
        {
            if(x < 0 || x > 1)
                throw new ArgumentException("x must be between 0 and 1");

            return MMath.RegularizedIncompleteBetaFunction(x, a, b);
        }

        /// <summary>
        /// Returns the PDF of a beta distribution with shape parameters a and b.
        /// f(x) = x^(a-1) * (1-x)^(b-1) / B(a,b)
        /// </summary>
        /// <param name="x">Value of the random variable for which the PDF is beign evaluated. x is between 0 and 1.</param>
        /// <param name="a">Number of trials.</param>
        /// <param name="b">Number of times the event occurs in n trials.</param>
        /// <returns></returns>
        public static double BetaProbabilityDensityFunction(double x, double a, double b)
        {
            if (x < 0 || x > 1)
                throw new ArgumentException("x must be between 0 and 1");

            double B = MMath.BetaFunction(a, b);
            double d = Math.Pow(x, a - 1) * Math.Pow(1-x, b-1) / B;
            return d;
        }
        #endregion

        #region Binomial Distribution
        /// <summary>
        /// Returns the probability of an event occuring in k or less times out of n trials, given 
        /// that the probability of the event occuring in any one trial is p.
        /// </summary>
        /// <param name="p">Probability that the event occurs in any one trial.</param>
        /// <param name="n">Number of trials.</param>
        /// <param name="k">Number of times the event occurs in n trials.</param>
        /// <returns></returns>
        public static double BinomialCumulativeDistributionFunction(double p, int n, int k)
        {
            if(k >= n)
                return 1.0;

            double sum = 0;
            for (int i = 0; i <= k; i++)
                sum += BinomialProbabilityDensityFunction(p, n, i);
            return sum;
        }

        /// <summary>
        /// Returns the probability of an event occuring in k out of n trials, given 
        /// that the probability of the event occuring in any one trial is p.
        /// </summary>
        /// <param name="p">Probability that the event occurs in any one trial.</param>
        /// <param name="n">Number of trials.</param>
        /// <param name="k">Number of times the event occurs in n trials.</param>
        /// <returns></returns>
        public static double BinomialProbabilityDensityFunction(double p, int n, int k)
        {
            //# of combinations can exceed max double for n > 1020
            if (n > 1020)
            {
                double logB = MMath.LogCombin(n, k) + k * Math.Log(p) + (n - k) * Math.Log(1 - p);
                return Math.Exp(logB);
            }

            return MMath.Combinations(n, k) * Math.Pow(p, k) * Math.Pow(1 - p, n - k);
        }
        #endregion

        #region Chi-Squared Distribution
        /// <summary>
        /// Chi-squared probability density function.
        /// </summary>
        /// <param name="x">The value at which the PDF is evaluated.</param>
        /// <param name="k">Degress of freedom, or number independent standard normal distributions.</param>
        /// <returns></returns>
        public static double ChiSquaredProbabilityDensityFunction(double x, int k)
        {
            double g = MMath.GammaFunction(0.5 * k);
            double a = 1.0 / (Math.Pow(2, 0.5 * k) * g);
            return a * Math.Pow(x, 0.5 * k - 1.0) * Math.Exp(-0.5 * x);
        }
        #endregion

        #region F Distribution
        /// <summary>
        /// Returns the PDF of the F distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="k1">Degrees of freedom for numerator chi-sqared distribution. k1 > 0.</param>
        /// <param name="k2">Degrees of freedom for denominator chi-sqared distribution. k2 > 0.</param>
        public static double FProbabilityDensityFunction(double x, int k1, int k2)
        {
            if (k1 <= 0 || k2 <= 0) throw new ArgumentException("k1 and k2 must be greater than 0.");
            if (x == 0) return 0.0;

            double a1 = Math.Pow(k1 * x, 0.5 * k1);
            double a2 = Math.Pow(k2, 0.5 * k2);
            double a3 = Math.Pow(k1 * x + k2, 0.5 * (k1 + k2));
            double a = a1 * a2 / a3;
            double b = MMath.BetaFunction(0.5 * k1, 0.5 * k2);
            double c = x * b;

            return a / c;
        }

        /// <summary>
        /// Returns the CDF of the F distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="k1">Degrees of freedom for numerator chi-sqared distribution.</param>
        /// <param name="k2">Degrees of freedom for denominator chi-sqared distribution.</param>
        public static double FCumulativeDensityFunction(double x, int k1, int k2)
        {
            //Tested against Excel for 
            //k1 & k2 = 1, 2, 3, ... , 9, 10, 20, ... , 100
            //x = 0.0, 0.1, 0.2, ... , 6.0
            //maximum discrepancy was less than 0.16% for (x, k1, k2) = (6, 1, 9)

            if (k1 <= 0 || k2 <= 0) throw new ArgumentException("k1 and k2 must be greater than 0.");
            if (x < 0) throw new ArgumentException("x must be greater than 0.");

            if (x == 0)
                return 0.0;
            double x2 = (k1 * x) / (k1 * x + k2);
            double p = MMath.RegularizedIncompleteBetaFunction(x2, 0.5 * k1, 0.5 * k2);
            return Math.Min(1, p);
        }

        /// <summary>
        /// Returns the inverse of the CDF of the F distribution.
        /// For k = 3, and 5+ the solution is an approximation. 
        /// </summary>
        /// <param name="p">Cumulative probability of the distribution. p is between 0 and 1.</param>
        /// <param name="k1">Degrees of freedom for numerator chi-sqared distribution.</param>
        /// <param name="k2">Degrees of freedom for denominator chi-sqared distribution.</param>
        public static double FCumulativeDistributionFunctionInverse(double p, int k1, int k2)
        {
            if (k1 <= 4 || k2 <= 4) throw new ArgumentException("k1 and k2 must be greater than 4.");

            //Tested against Excel for k1 & k2 = 5, 6, ... 19, 20, 30, 40, ... , 100 and p = 0.01, 0.02, ..., 0.99.
            //Max absolute difference was 0.0313 for (p, k1, k2) = (0.99, 7, 5).

            const int maxIterations = 16;
            const double minAccuracy = 0.0001; //min absolute difference between p and p implied by return value.

            if (p < 0 || p > 1.0) throw new ArgumentException("Probability must be between 0 and 1");
            if (p == 0) return 0.0;
            if (p == 1) return double.PositiveInfinity;

            double mu = k2 / (k2 - 2.0);
            double sigma = 2 * k2 * k2 * (k1 + k2 - 2.0) / (k1 * (k2 - 2.0) * (k2 - 2.0) * (k2 - 4.0));
            double normAprox = Math.Max(0, NormalCumulativeDistributionFunctionInverse(p, mu, sigma));

            double x;
            if (k1 < 10 && k2 < 10)
                if (p < 0.03)
                    x = 0.1592355;
                else
                    x = 0.1819022 * normAprox + 6.0340494 * Math.Log(p) / Math.Sqrt(k1 + k2) - 10.4102348 / Math.Sqrt(k1 + k2) - 3.8062818 * Math.Log(p) / k2 - 2.1284338 * Math.Log(p) / k1 - 0.0446516 * k2 - 0.1014539 * k1 - 7.980049829 * p + 7.988783845 * p * p + 6.864968931;
            else if (k1 < 10)
                x = 0.8100214 * normAprox + 3.1485113 * Math.Log(p) / Math.Sqrt(k1 + k2) - 5.63992350 / Math.Sqrt(k1 + k2) - 6.8494933 * Math.Log(p) / k2 - 0.5725057 * Math.Log(p) / k1 - 0.0102242 * k2 - 0.0099565 * k1 - 2.820859600 * p + 3.036058600 * p * p + 2.103689000;
            else if (k2 < 10)
                x = 0.1944974 * normAprox + 5.5846823 * Math.Log(p) / Math.Sqrt(k1 + k2) - 2.6254457 / Math.Sqrt(k1 + k2) - 0.7057866 * Math.Log(p) / k2 - 6.49649350 * Math.Log(p) / k1 + 0.0048548 * k2 - 0.0088253 * k1 - 6.523519381 * p + 6.710305803 * p * p + 3.305158009;
            else
                x = 0.8717865 * normAprox + 0.9139534 * Math.Log(p) / Math.Sqrt(k1 + k2) - 0.0864391 / Math.Sqrt(k1 + k2) - 1.5438103 * Math.Log(p) / k2 - 0.14952740 * Math.Log(p) / k1 + 0.0004922 * k2 - 0.0003965 * k1 - 0.636234600 * p + 1.177650500 * p * p + 0.103934500;

            x = Math.Max(0.00001, x);
            double oldAbsAccuracy = double.MaxValue;
            double step = 0;
            for (int i = 0; i < maxIterations; i++)
            {
                double p0 = FCumulativeDensityFunction(x, k1, k2);
                double absAccuracy = Math.Abs(p - p0);
                if (absAccuracy < minAccuracy) return x;
                if (absAccuracy > oldAbsAccuracy)
                {
                    x = x - 0.5 * step;
                    p0 = FCumulativeDensityFunction(x, k1, k2);
                    absAccuracy = Math.Abs(p - p0);
                }
                oldAbsAccuracy = absAccuracy;
                double oldX = x;

                double slope = FProbabilityDensityFunction(x, k1, k2);
                step = (p - p0) / slope;
                if (step < -x) step = -0.5 * x;
                if (step > 1.0) step = 1.0;
                x = x + step;

                if (x < 0) x = 0.5 * oldX;
            }
            throw new Exception("Solution did not converge");
        }

        //currently not used ...might be useful for future version of FCumulativeDistributionFunctionInverse()
        private static double DerivativeOfFProbabilityDensityFunction(double x, double k1, double k2, double pdf)
        {
            //double pdf = FProbabilityDensityFunction(x, k1, k2);
            double a = (k1 * k2 * (1.0 - x)) / (k1 * x + k2);
            return pdf * (1.0 / x) * (0.5 * a - 1.0);
        }
        #endregion

        #region Non-Central Chi-Squared Distribution
        /// <summary>
        /// Non-central chi-squared cumulative distribution function.
        /// Approximation based on Sankaran (1963). 
        /// </summary>
        /// <param name="x">The value at which the CDF is evaluated. Must be non-negative.</param>
        /// <param name="k">degrees of freedom; cannot be negative</param>
        /// <param name="lambda">non-centrality parameter; cannot be negative</param>
        /// <returns></returns>
        public static double NonCentralChiSquaredCumulativeDistributionFunction(double x, double k, double lambda)
        {
            if (k < 0) throw new ArgumentException("degrees of freedom cannot be negative");
            if (lambda < 0) throw new ArgumentException("non-centrality parameter cannot be negative");
            if (x < 0) throw new ArgumentException("x cannot be negative");

            if (x == 0) return 0.0;
            if (x == double.PositiveInfinity) return 1.0;

            double h = 1.0 - (2.0/3.0)*(k + lambda)*(k + 3.0*lambda)/Math.Pow(k + 2.0*lambda, 2);
            double p = (k + 2.0 * lambda) / Math.Pow(k + lambda, 2);
            double m = (h - 1.0)*(1.0 - 3.0*h);
            double f = Math.Pow(x/(k + lambda), h);
            double s = 1.0 + h*p*(h - 1.0 - 0.5*(2.0 - h)*m*p);
            double d = h*(1.0 + 0.5*m*p)*Math.Sqrt(2.0*p);
            double z = (f - s)/d;
            return StandardNormalCumulativeDistributionFunction(z);
        }
        #endregion

        #region Normal Distribution
        /// <summary>
        /// Returns the PDF of the normal distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="mean">Mean of the distribution.</param>
        /// <param name="standardDeviation">Standard deviation of the distribution.</param>
        public static double NormalProbabilityDensityFunction(double x, double mean, double standardDeviation)
        {
            const double sqrtTwoPiInv = 0.398942280401433;
            double z = (x - mean) / standardDeviation;
            return sqrtTwoPiInv * Math.Exp(-0.5 * z * z) / standardDeviation;
        }

        /// <summary>
        /// Returns the CDF of the normal distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="mean">Mean of the distribution.</param>
        /// <param name="standardDeviation">Standard deviation of the distribution.</param>
        public static double NormalCumulativeDistributionFunction(double x, double mean, double standardDeviation)
        {
            double z = (x - mean) / standardDeviation;
            return StandardNormalCumulativeDistributionFunction(z);
        }

        /// <summary>
        /// Returns the inverse of the CDF of the normal distribution.
        /// </summary>
        /// <param name="p">Cumulative probability of the distribution. 0 &lt;= p &gt;= 1.</param>
        ///<param name="mean">Mean of the distribution.</param>
        /// <param name="standardDeviation">Standard deviation of the distribution.</param>
        public static double NormalCumulativeDistributionFunctionInverse(double p, double mean, double standardDeviation)
        {
            double x = StandardNormalCumulativeDistributionFunctionInverse(p);
            return standardDeviation * x + mean;
        }

        /// <summary>
        /// Returns the PDF of the standard normal distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        public static double StandardNormalProbabilityDensityFunction(double x)
        {
            const double SqrtTwoPiInv = 0.398942280401433;
            return SqrtTwoPiInv * Math.Exp(-0.5 * x * x);
        }

        /// <summary>
        /// Returns the CDF of the standard normal distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        public static double StandardNormalCumulativeDistributionFunction(double x)
        {
            //Approimation based on Abramowitz & Stegun (1964)

            if (x < 0)
                return 1.0 - StandardNormalCumulativeDistributionFunction(-x);
            const double b0 = 0.2316419; 
            const double b1 = 0.319381530;
            const double b2 = -0.356563782; 
            const double b3 = 1.781477937;
            const double b4 = -1.821255978;
            const double b5 = 1.330274429;
            double pdf = StandardNormalProbabilityDensityFunction(x);
            double a = 1.0 / (1.0 + b0 * x);
            return 1.0 - pdf * (b1 * a + b2 * Math.Pow(a, 2) + b3 * Math.Pow(a, 3) + b4 * Math.Pow(a, 4) + b5 * Math.Pow(a, 5));
        }

        /// <summary>
        /// Returns the inverse of the CDF of the standard normal distribution.
        /// </summary>
        /// <param name="p">Cumulative probability of the distribution. p is between 0 and 1.</param>
        public static double StandardNormalCumulativeDistributionFunctionInverse(double p)
        {
            //based on pseudo-code from http://home.online.no/~pjacklam/notes/invnorm/

            if (p < 0 || p > 1.0) throw new ArgumentException("Probability must be between 0 and 1");
            if (p == 0) return double.NegativeInfinity;
            if (p == 1) return double.PositiveInfinity;
            if (p == 0.5) return 0.0;

            const double a1 = -3.969683028665376e+01;
            const double a2 =  2.209460984245205e+02;
            const double a3 = -2.759285104469687e+02;
            const double a4 =  1.383577518672690e+02;
            const double a5 = -3.066479806614716e+01;
            const double a6 =  2.506628277459239e+00;

            const double b1 = -5.447609879822406e+01;
            const double b2 =  1.615858368580409e+02;
            const double b3 = -1.556989798598866e+02;
            const double b4 =  6.680131188771972e+01;
            const double b5 = -1.328068155288572e+01;

            const double c1 = -7.784894002430293e-03;
            const double c2 = -3.223964580411365e-01;
            const double c3 = -2.400758277161838e+00;
            const double c4 = -2.549732539343734e+00;
            const double c5 =  4.374664141464968e+00;
            const double c6 =  2.938163982698783e+00;

            const double d1 =  7.784695709041462e-03;
            const double d2 =  3.224671290700398e-01;
            const double d3 =  2.445134137142996e+00;
            const double d4 =  3.754408661907416e+00;

            const double pLow  = 0.02425;
            const double pHigh = 1 - pLow;

            double q, x;
            if (0 < p && p < pLow)
            {
                q = Math.Sqrt(-2*Math.Log(p));
                x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
                return x;
            }
      
            if (pLow <= p && p <= pHigh)
            {
                q = p - 0.5;
                double r = q*q;
                x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
                return x;
            }

            //if(p_high < p && p < 1)
            q = Math.Sqrt(-2*Math.Log(1-p));
            x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
            return x;
        }
        #endregion

        #region Student's t Distribution
        /// <summary>
        /// Returns the PDF of Student's t distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="k">Degrees of freedom.</param>
        public static double StudentsTProbabilityDensityFunction(double x, double k)
        {
            double a = MMath.GammaFunction(0.5 * k + 0.5);
            double b = Math.Pow(1.0 + (x * x) / k, -0.5 * k - 0.5);
            double c = Math.Sqrt(k * Math.PI) * MMath.GammaFunction(0.5 * k);
            
            return a * b / c;
        }

        /// <summary>
        /// Returns the CDF of Student's t distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="k">Degrees of freedom.</param>
        public static double StudentsTCumulativeDensityFunction(double x, double k)
        {
            //Tested against Excel for 
            //df = 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200, 30, 400, 500 and
            //x = -7.0, -6.9, -6.8, ..., 7.0
            //maximum discrepancy was less than 0.04%
            if(double.IsNaN(x))
                return double.NaN;
            if (x == 0)
                return 0.5;
            if (x > 0)
            {
                double x2 = k / (x * x + k);
                return 1.0 - 0.5 * MMath.RegularizedIncompleteBetaFunction(x2, 0.5 * k, 0.5);
            }
            return 1.0 - StudentsTCumulativeDensityFunction(-x, k);
        }

        /// <summary>
        /// Returns the inverse of the CDF of the Student's t distribution.
        /// For k = 3, and 5+ the solution is an approximation. 
        /// </summary>
        /// <param name="p">Cumulative probability of the distribution. p is between 0 and 1.</param>
        /// <param name="k">Degrees of freedom.</param>
        public static double StudentsTCumulativeDistributionFunctionInverse(double p, double k)
        {
            //Tested against Excel for k = 1, 2, ... 20, 30, ..., 100 and p = 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.02, ..., 0.99.
            //Max absolute difference was 0.0088 for k = 3 and p = 0.000001 (-103.2906 vs -103.2995).

            const int maxIterations = 32;
            double maxProbError = Math.Min(0.0001, 0.001 * Math.Min(p, 1 - p)); //max absolute difference between p and p implied by return value.
            const double maxXDiff = 0.00000001; //If bisection method is called, it will stop when the difference between the two values it is testing is less than maxXDiff.

            if (p < 0 || p > 1.0) throw new ArgumentException("Probability must be between 0 and 1");
            if (k <= 0) throw new ArgumentException("Degrees of freedom must be greater than 0.");
            if (p == 0) return double.NegativeInfinity;
            if (p == 1) return double.PositiveInfinity;
            if (p == 0.5) return 0.0;

            if (k == 1)
            {
                return Math.Tan(Math.PI*(p - 0.5));
            }
            else if (k == 2)
            {
                double alpha = 4*p*(1 - p);
                return 2*(p - 0.5)*Math.Sqrt(2.0/alpha);
            }
            else if (k == 4)
            {
                double alpha = 4*p*(1 - p);
                double sqrtAlpha = Math.Sqrt(alpha);
                double q = Math.Cos((1.0/3.0)*Math.Acos(sqrtAlpha))/sqrtAlpha;
                return Math.Sign(p - 0.5)*2*Math.Sqrt(q - 1);
            }
            else
            {
                double x;

                if (k > 5.5 && (p < 0.017 || p > 0.983))
                    x = StudentsTCumulativeDistributionFunctionInverseApproximation(p, k); //this method is more accurate for extreme values pf p.
                else if (k > 7.875)
                    x = StandardNormalCumulativeDistributionFunctionInverse(p);
                else if (k < 1.375)
                    x = StudentsTCumulativeDistributionFunctionInverse(p, 1);
                else if (k < 2.625)
                    x = StudentsTCumulativeDistributionFunctionInverse(p, 2);
                else
                    x = StudentsTCumulativeDistributionFunctionInverse(p, 4);

                double oldAbsProbError = double.MaxValue;
                double oldX = double.NaN;
                for (int i = 0; i < maxIterations; i++)
                {
                    double p0 = StudentsTCumulativeDensityFunction(x, k);
                    double absProbError = Math.Abs(p - p0);
                    if (absProbError < maxProbError) return x;

                    if (absProbError > oldAbsProbError)
                    {
                        bool reachedDiffLimit = false;
                        //x = BisectionSearch(x, oldX, k, p, i, maxIterations, maxXDiff, ref reachedDiffLimit, null, null);
                        x = SecantSearch(x, oldX, k, p, i, maxIterations, maxXDiff, ref reachedDiffLimit, null, null);
                        if (reachedDiffLimit) return x;
                        p0 = StudentsTCumulativeDensityFunction(x, k);
                        absProbError = Math.Abs(p - p0);
                        if (absProbError < maxProbError) return x;
                        throw new Exception("Solution did not converge");
                    }
                    oldAbsProbError = absProbError;
                    oldX = x;

                    double slope = StudentsTProbabilityDensityFunction(x, k);
                    x = x + (p - p0)/slope;
                    if ((p < 0.50 && x > 0) || (p > 0.50 && x < 0))
                        x = 0.5*oldX; //if you overshoot so much that you end up with the wrong sign for x, then take the average of the last x and 0.
                }
                throw new Exception("Solution did not converge");
            }
        }

        private static double BisectionSearch(double x1, double x2, double k, double p, int splitCount, int maxSplits, double maxXDiff, ref bool reachedDiffLimit, double? e1, double? e2)
        {
            if (Math.Abs(x1 - x2) < maxXDiff)
            {
                reachedDiffLimit = true;
                return 0.5 * (x1 + x2);
            }

            double xm = 0.5*(x1 + x2);
            e1 = e1 ?? StudentsTCumulativeDensityFunction(x1, k) - p;
            double em = StudentsTCumulativeDensityFunction(xm, k) - p;
            e2 = e2 ?? StudentsTCumulativeDensityFunction(x2, k) - p;

            splitCount++;
            if (splitCount > maxSplits)
                return xm;

            if (Math.Sign(e1.Value) != Math.Sign((em)))
            {
                xm = BisectionSearch(x1, xm, k, p, splitCount, maxSplits, maxXDiff, ref reachedDiffLimit, e1, em);
            }
            else if (Math.Sign(e2.Value) != Math.Sign(em))
            {
                xm = BisectionSearch(x2, xm, k, p, splitCount, maxSplits, maxXDiff, ref reachedDiffLimit, e2, em);
            }
            return xm;
        }

        private static double SecantSearch(double x1, double x2, double k, double p, int splitCount, int maxSplits, double maxXDiff, ref bool reachedDiffLimit, double? e1, double? e2)
        {
            if (Math.Abs(x1 - x2) < maxXDiff)
            {
                reachedDiffLimit = true;
                return 0.5 * (x1 + x2);
            }

            double xm = e1 == null ? 0.5 * (x1 + x2) : (x1 * e2.Value - x2 * e1.Value) / (e2.Value - e1.Value);
            e1 = e1 ?? StudentsTCumulativeDensityFunction(x1, k) - p;
            double em = StudentsTCumulativeDensityFunction(xm, k) - p;
            e2 = e2 ?? StudentsTCumulativeDensityFunction(x2, k) - p;

            splitCount++;
            if (splitCount > maxSplits)
                return xm;

            if (Math.Sign(e1.Value) != Math.Sign((em)))
            {
                xm = BisectionSearch(x1, xm, k, p, splitCount, maxSplits, maxXDiff, ref reachedDiffLimit, e1, em);
            }
            else if (Math.Sign(e2.Value) != Math.Sign(em))
            {
                xm = BisectionSearch(x2, xm, k, p, splitCount, maxSplits, maxXDiff, ref reachedDiffLimit, e2, em);
            }
            return xm;
        }

        /// <summary>
        /// This approximation is more accurate for very high and very low values of p.
        /// </summary>
        /// <param name="p"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        private static double StudentsTCumulativeDistributionFunctionInverseApproximation(double p, double k)
        {
            if(p < 0.5)
                return -StudentsTCumulativeDistributionFunctionInverseApproximation(1 - p, k);

            double[] paramArray = StcdfiaParams(p);
            return Math.Exp(paramArray[0] + paramArray[1] / k + paramArray[2] / (k * k));
        }

        /// <summary>
        /// These constants, {A, B, C}, are based on regressions of the form x = A + B/k + C/k^2 for various value of p.
        /// </summary>
        private static double[] StcdfiaParams(double p)
        {
            if(p < 0.994500000)
                return new[] {0.844289344098942, 1.60247403675763, 1.21204210950434};       //p = 99%
            if(p < 0.999450000)
                return new[] {1.12763432956672, 2.6724693160194, 2.75261355460215};         //p = 99.9%
            if(p < 0.999945000)
                return new[] {1.3109905400331, 3.85826102579765, 4.54414522746746};         //p = 99.99%
            if(p < 0.999994500)
                return new[] {1.44451562705915, 5.1673200281153, 6.35867233089288};         //p = 99.999%
            if(p < 0.999999450)
                return new[] {1.54796181658583, 6.59696293772306, 8.07465724847696};        //p = 99.9999%
            if(p < 0.999999945)
                return new[] {1.63123569060444, 8.13645869984937, 9.64357980688031};        //p = 99.99999%
            return new[] { 1.70006767767707, 9.77199585223768, 11.0556716544646 };          //p = 99.999999%
        }
        #endregion

        #region Uniform Distribution
        /// <summary>
        /// Returns the PDF of the uniform distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="min">Minimum value of the distribution.</param>
        /// <param name="max">Maximum value of the distribution.</param>
        public static double UniformProbabilityDensityFunction(double x, double min, double max)
        {
            if (max <= min) throw new ArgumentException("max must be greater than min");
            if (x < min || x > max) return 0.0;
            return 1.0 / (max - min);
        }

        /// <summary>
        /// Returns the CDF of the uniform distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="min">Minimum value of the distribution.</param>
        /// <param name="max">Maximum value of the distribution.</param>
        public static double UniformCumulativeDensityFunction(double x, double min, double max)
        {
            if (max <= min) throw new ArgumentException("max must be greater than min");
            return (x - min) / (max - min);
        }
        #endregion

        /// <summary>
        /// The Irwin-Hall distribution results from the sum on n independent standard uniform variables
        /// </summary>
        /// <param name="x">The value at which to evaluate the distribution.</param>
        /// <param name="n">The number of standard uniform variables.</param>
        public static double IrwinHallProbabilityDensityFunction(double x, int n)
        {
            if (x < 0 || x > n) return 0.0;
            double d = 0;
            for (int i = 0; i <= n; i++)
                d += Math.Pow(-1, i) * MMath.Combinations(n, i) * Math.Pow(x - i, n - 1) * Math.Sign(x - i);
            d *= 0.5;
            d /= MMath.Factorial(n - 1);
            return d;
        }

        /// <summary>
        /// Returns the PDF of the lognormal distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="mu">This is *not* the mean. See Miller (2012).</param>
        /// <param name="sigma">This is *not* the standard deviation. See Miller (2012).</param>
        public static double LogNormalProbabilityDensityFunction(double x, double mu, double sigma)
        {
            const double sqrtTwoPiInv = 0.398942280401433;
            double z = (Math.Log(x) - mu) / sigma;
            return sqrtTwoPiInv * Math.Exp(-0.5 * z * z) / (x * sigma);
        }

        /// <summary>
        /// The PDF of the Poisson distribution.
        /// The probability of n iid events occuring in an interval, if the expected number of event is lambda.
        /// </summary>
        /// <param name="lambda">Equals both the mean and variance of the distribution.</param>
        /// <param name="n">The number of events.</param>
        public static double PoissonProbabilityDensityFunction(double lambda, int n)
        {
            return Math.Pow(lambda, n) * Math.Exp(-lambda) / MMath.Factorial(n);
        }

        /// <summary>
        /// Maps the input array to a standard normal distribution. 
        /// The points in the output array have the same rank order as the input array, but are normally distributed.
        /// Note, if two points in the input array have the same value, the values in the output array will be different.
        /// </summary>
        public static double[] NormalizeArray(double[] array)
        {
            List<NormalizePoint> points = new List<NormalizePoint>();
            int n = array.Length;
            for(int i = 0; i < n; i++)
                points.Add(new NormalizePoint(i, array[i]));
            points.Sort((a,b) => a.OriginalValue.CompareTo(b.OriginalValue));

            double interval = 1.0 / n;
            for(int i = 0; i < n; i++)
                points[i].NormalizedValue = StandardNormalCumulativeDistributionFunctionInverse(interval * (i + 0.5));

            points.Sort((a, b) => a.OriginalIndex.CompareTo(b.OriginalIndex));
            double[] outArray = new double[n];
            for(int i = 0; i < n; i++)
                outArray[i] = points[i].NormalizedValue;

            return outArray;
        }

        private class NormalizePoint
        {
            public readonly int OriginalIndex;
            public readonly double OriginalValue;
            public double NormalizedValue;

            public NormalizePoint(int originalIndex, double value)
            {
                OriginalIndex = originalIndex;
                OriginalValue = value;
            }
        }
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.
