using System;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    
    public partial class Distributions
    {
        #region Pareto Distribution
        /// <summary>
        /// Returns the PDF of the Pareto distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="min">Minimum value of the distribution. Also known as the scale parameter.</param>
        /// <param name="alpha">Shape parameter.</param>
        public static double ParetoProbabilityDensityFunction(double x, double min, double alpha)
        {
            if (alpha <= 0) throw new ArgumentException("alpha must be greater than zero.");
            if (x < min) return 0.0;
            return alpha * Math.Pow(min / x, alpha) / x;
        }

        /// <summary>
        /// Returns the CDF of the Pareto distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="min">Minimum value of the distribution. Also known as the scale parameter.</param>
        /// <param name="alpha">Shape parameter. Must be greater than 0.</param>
        public static double ParetoCumulativeDensityFunction(double x, double min, double alpha)
        {
            if (alpha <= 0) throw new ArgumentException("alpha must be greater than zero.");
            if(x < min) return 0.0;
            return 1 - Math.Pow(min / x, alpha);
        }

        /// <summary>
        /// Returns the inverse of the CDF of the Pareto distribution.
        /// </summary>
        /// <param name="p">Cumulative probability of the distribution. 0 &lt;= p &gt;= 1.</param>
        /// <param name="min">Minimum value of the distribution. Also known as the scale parameter.</param>
        /// <param name="alpha">Shape parameter. Must be greater than 0.</param>
        public static double ParetoCumulativeDensityFunctionInverse(double p, double min, double alpha)
        {
            if (alpha <= 0) throw new ArgumentException("alpha must be greater than zero.");
            if (p < 0 || p > 1) throw new ArgumentException("p is a probability and must be between 0 and 1, inclusive.");
            return min * Math.Pow(1 - p, -1.0 / alpha);
        }
        #endregion

        #region Gumbel Distribution
        /// <summary>
        /// Returns the PDF of the Gumbel distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="mu">Location parameter.</param>
        /// <param name="sigma">Scale parameter. Must be greater than 0.</param>
        public static double GumbelProbabilityDensityFunction(double x, double mu, double sigma)
        {
            if (sigma <= 0) throw new ArgumentException("sigma must be greater than zero.");
            double z = (x - mu) / sigma;
            return Math.Exp(-(z + Math.Exp(-z))) / sigma;
        }

        /// <summary>
        /// Returns the CDF of the Gumbel distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="mu">Location parameter.</param>
        /// <param name="sigma">Scale parameter. Must be greater than 0.</param>
        public static double GumbelCumulativeDensityFunction(double x, double mu, double sigma)
        {
            if (sigma <= 0) throw new ArgumentException("sigma must be greater than zero.");
            double z = (x - mu) / sigma;
            return Math.Exp(-Math.Exp(-z));
        }

        /// <summary>
        /// Returns the inverse of the CDF of the Gumbel distribution.
        /// </summary>
        /// <param name="p">Cumulative probability of the distribution. 0 &lt;= p &gt;= 1.</param>
        /// <param name="mu">Location parameter.</param>
        /// <param name="sigma">Scale parameter. Must be greater than 0.</param>
        public static double GumbelCumulativeDensityFunctionInverse(double p, double mu, double sigma)
        {
            if (sigma <= 0) throw new ArgumentException("sigma must be greater than zero.");
            if (p < 0 || p > 1) throw new ArgumentException("p is a probability and must be between 0 and 1, inclusive.");
            double z = -Math.Log(-Math.Log(p));
            return sigma * z + mu;
        }
        #endregion

        #region Fréchet Distribution
        /// <summary>
        /// Returns the PDF of the Fréchet distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="mu">Location parameter.</param>
        /// <param name="sigma">Scale parameter. Must be greater than 0.</param>
        /// <param name="alpha">Shape parameter.</param>
        public static double FrechetProbabilityDensityFunction(double x, double mu, double sigma, double alpha)
        {
            if (sigma <= 0) throw new ArgumentException("sigma must be greater than zero.");
            if(x <= mu) return 0.0;

            double z = (x - mu) / sigma;
            return (alpha / sigma) * Math.Pow(z, -1 - alpha) * Math.Exp(-Math.Pow(z, -alpha));
        }

        /// <summary>
        /// Returns the CDF of the Fréchet distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="mu">Location parameter.</param>
        /// <param name="sigma">Scale parameter. Must be greater than 0.</param>
        /// <param name="alpha">Shape parameter.</param>
        public static double FrechetCumulativeDensityFunction(double x, double mu, double sigma, double alpha)
        {
            if (sigma <= 0) throw new ArgumentException("sigma must be greater than zero.");
            if (x <= mu) return 0.0;

            double z = (x - mu) / sigma;
            return Math.Exp(-Math.Pow(z, -alpha));
        }

        /// <summary>
        /// Returns the inverse of the CDF of the Fréchet distribution.
        /// </summary>
        /// <param name="p">Cumulative probability of the distribution. 0 &lt;= p &gt;= 1.</param>
        /// <param name="mu">Location parameter.</param>
        /// <param name="sigma">Scale parameter. Must be greater than 0.</param>
        /// <param name="alpha">Shape parameter.</param>
        public static double FrechetCumulativeDensityFunctionInverse(double p, double mu, double sigma, double alpha)
        {
            if (sigma <= 0) throw new ArgumentException("sigma must be greater than zero.");
            if (p < 0 || p > 1) throw new ArgumentException("p is a probability and must be between 0 and 1, inclusive.");
            double z = Math.Pow(-Math.Log(p), -1.0 / alpha);
            return sigma * z + mu;
        }
        #endregion

        #region Reversed Weibull Distribution
        /// <summary>
        /// Returns the PDF of the Reversed Weibull distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="mu">Location parameter.</param>
        /// <param name="sigma">Scale parameter. Must be greater than 0.</param>
        /// <param name="alpha">Shape parameter.</param>
        public static double ReversedWeibullProbabilityDensityFunction(double x, double mu, double sigma, double alpha)
        {
            if (sigma <= 0) throw new ArgumentException("sigma must be greater than zero.");
            if (x >= mu) return 0.0;

            double z = (x - mu) / sigma;
            return (alpha / sigma) * Math.Exp(-Math.Pow(-z, alpha)) * Math.Pow(-z, alpha - 1);
        }

        /// <summary>
        /// Returns the CDF of the Reversed Weibull distribution.
        /// </summary>
        /// <param name="x">Value at which the distribution is evaluated.</param>
        /// <param name="mu">Location parameter.</param>
        /// <param name="sigma">Scale parameter. Must be greater than 0.</param>
        /// <param name="alpha">Shape parameter.</param>
        public static double ReversedWeibullCumulativeDensityFunction(double x, double mu, double sigma, double alpha)
        {
            if (sigma <= 0) throw new ArgumentException("sigma must be greater than zero.");
            if (x >= mu) return 1.0;

            double z = (x - mu) / sigma;
            return Math.Exp(-Math.Pow(-z, alpha));
        }

        /// <summary>
        /// Returns the inverse of the CDF of the Reversed Weibull distribution.
        /// </summary>
        /// <param name="p">Cumulative probability of the distribution. 0 &lt;= p &gt;= 1.</param>
        /// <param name="mu">Location parameter.</param>
        /// <param name="sigma">Scale parameter. Must be greater than 0.</param>
        /// <param name="alpha">Shape parameter.</param>
        public static double ReversedWeibullCumulativeDensityFunctionInverse(double p, double mu, double sigma, double alpha)
        {
            if (sigma <= 0) throw new ArgumentException("sigma must be greater than zero.");
            if (p < 0 || p > 1) throw new ArgumentException("p is a probability and must be between 0 and 1, inclusive.");
            double z = -Math.Pow(-Math.Log(p), 1.0 / alpha);
            return sigma * z + mu;
        }
        #endregion
    }
}


//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.
