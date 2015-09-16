using System;
using System.Collections.Generic;
using System.Linq;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///
    ///References:
    ///Multivariate Distributions and Copulas
    ///available at www.risk256.com/writing

    public class Copulas
    {
        #region Rank based measures of dependence
        /// <summary>
        /// Returns Kendall's tau.
        /// P[concordance]-P[discordance]
        /// Ties are ignored.
        /// </summary>
        public static double KendallsTau(double[] array1, double[] array2) 
        {
            if(array1.Length != array2.Length)
                throw new ArgumentException("Arrays must be the same length");

            int n = array1.Length;
            double tau = 0;
            for(int i = 1; i < n; i++)
                for (int j = 0; j < i; j++)
                {
                    tau += Math.Sign(array1[i] - array1[j]) * Math.Sign(array2[i] - array2[j]);
                }
            int combin = n * (n-1) / 2;
            return tau / combin;
        }

        /// <summary>
        /// Calculates Spearman's Rho.
        /// Correlation of ranks of array elements.
        /// </summary>
        public static double SpearmansRho(double[] array1, double[] array2)
        {
            if (array1.Length != array2.Length)
                throw new ArgumentException("Arrays must be the same length");

            double[] rank1 = RankArray(array1);
            double[] rank2 = RankArray(array2);
            return Moments.Correlation(rank1, rank2);
        }

        /// <summary>
        /// Returns an array of the rank/order of the numbers in an array.
        /// The rank of the greastest number equals 1.
        /// Ties are given a rank equal to the average of their position (1 based index) in an ordered array.
        /// e.g. [30, 10, 10, 20] => [1, 3.5, 3.5, 2]
        /// </summary>
        /// <param name="array"></param>
        /// <returns></returns>
        public static double[] RankArray(double[] array)
        {
            List<double> list = array.ToList();
            List<double> sortList = array.ToList();
            sortList.Sort();
            sortList.Reverse();

            int[] simpleRankArray = new int[array.Length];
            for (int i = 0; i < simpleRankArray.Length; i++)
                simpleRankArray[i] = sortList.IndexOf(list[i]) + 1;

            double[] rankArray = new double[array.Length];
            for (int i = 0; i < simpleRankArray.Length; i++)
                rankArray[i] = FractionalRank(simpleRankArray[i], simpleRankArray.Count(lmb => lmb == simpleRankArray[i]));
            return rankArray;
        }

        private static double FractionalRank(int simpleRank, int n)
        {
            return simpleRank + 0.5 * (n - 1);
        }
        #endregion

        public enum CopulaType
        {
            Clayton, FGM, Frank, Gumbel, Independent, Joe
        }

        /// <summary>
        /// Returns the value of a copula, given two input cumulative distribution fucntions, u and v.
        /// </summary>
        public static double CopulaCumulativeDistributionFunction(CopulaType type, double u, double v, double alpha)
        {
            CheckCopulaAlpha(type, alpha);

            switch (type)
            {
                case CopulaType.Clayton:
                    return Math.Pow(Math.Pow(u, -alpha) + Math.Pow(v, -alpha) - 1.0, -1.0 / alpha);
                case CopulaType.FGM:
                    return u * v * (1.0 + alpha * (1.0 - u) * (1.0 - v));
                case CopulaType.Frank:
                    double d1 = 1 + (FrankPart(u, alpha) * FrankPart(v, alpha)) / FrankPart(1.0, alpha);
                    return -(1.0/alpha) * Math.Log(d1);
                case CopulaType.Gumbel:
                    double d2 = -Math.Pow(Math.Pow(-Math.Log(u), alpha) + Math.Pow(-Math.Log(v), alpha), 1.0 / alpha);
                    return Math.Exp(d2);
                case CopulaType.Independent:
                    return u * v;
                case CopulaType.Joe:
                    double ju = Math.Pow(1.0 - u, alpha);
                    double jv = Math.Pow(1.0 - v, alpha);
                    double d3 = ju + jv - ju * jv;
                    return 1 - Math.Pow(d3, 1.0 / alpha);
                default:
                    throw new ArgumentException("Copula type not expected.");
            }
        }

        /// <summary>
        /// Returns the value of the density function of a copula, given two input cumulative distribution fucntions, u and v.
        /// </summary>
        public static double CopulaDensityFunction(CopulaType type, double u, double v, double alpha)
        {
            CheckCopulaAlpha(type, alpha);

            switch (type)
            {
                case CopulaType.Clayton:
                    return (1.0 + alpha) * Math.Pow(u * v, -alpha - 1.0) * Math.Pow(Math.Pow(u, -alpha) + Math.Pow(v, -alpha) - 1.0, -1.0 / alpha - 2.0);
                case CopulaType.FGM:
                    return 1.0 + alpha * (1.0 - 2.0 * u) * (1.0 - 2.0 * v);
                case CopulaType.Frank:
                    double d1 = FrankPart(u, alpha) * FrankPart(v, alpha) + FrankPart(1, alpha);
                    return -alpha * FrankPart(1, alpha) * Math.Exp(-alpha * (u + v)) / (d1 * d1);
                case CopulaType.Gumbel:
                    double lnu = Math.Log(u);
                    double lnv = Math.Log(v);
                    double g1 = Math.Pow(lnu * lnv, alpha - 1.0) / (u * v);
                    double g2 = CopulaCumulativeDistributionFunction(type, u, v, alpha);
                    double g3 = Math.Pow(-lnu + -lnv, 1.0 / alpha - 2.0);
                    double g4 = alpha - 1.0 + Math.Pow(-lnu + -lnv, 1.0 / alpha);
                    return g1 * g2 * g3 * g4;
                case CopulaType.Independent:
                    return 1.0;
                case CopulaType.Joe:
                    double ju = Math.Pow(1.0 - u, alpha);
                    double jv = Math.Pow(1.0 - v, alpha);
                    double d3 = ju + jv - ju * jv;
                    return Math.Pow(1.0 - u, alpha - 1.0) * Math.Pow(1.0 - v, alpha) * Math.Pow(d3, (1.0 / alpha) - 2.0) * (1.0 - alpha - d3);
                default:
                    throw new ArgumentException("Copula type not expected.");
            }
        }

        /// <summary>
        /// Returns the inverse of the first marginal CDF of a copula, for two cumulative distribution fucntions, u and v.
        /// Used in Monte Carlo simulation.
        /// </summary>
        /// <param name="type"></param>
        /// <param name="u">random number [0-1]</param>
        /// <param name="C1">random number [0-1], independent of u</param>
        /// <param name="alpha"></param>
        /// <returns>v, random number [0-1]</returns>
        public static double CopulaFirstMarginalInverse(CopulaType type, double u, double C1, double alpha)
        {
            CheckCopulaAlpha(type, alpha);

            switch (type)
            {
                case CopulaType.Clayton:
                    double d = Math.Pow(C1, -alpha / (1.0 + alpha)) + Math.Pow(u, alpha) - 1.0;
                    return u * Math.Pow(d, -1.0 / alpha);
                case CopulaType.Frank:
                    double f1 = C1 * (Math.Exp(-alpha) - 1.0);
                    double f2 = 1.0 + (Math.Exp(-alpha * u) - 1.0) * (1.0 - C1);
                    return -(1.0 / alpha) * Math.Log(1.0 + f1 / f2);
                case CopulaType.Independent:
                    return C1;
                default:
                    throw new ArgumentException("Copula type not expected.");
            }
        }

        private static void CheckCopulaAlpha(CopulaType type, double alpha)
        {
            switch (type)
            {
                case CopulaType.Clayton:
                    if (alpha <= 0.0) throw new ArgumentException("Invalid alpha. Alpha should be greater than zero.");
                    return;
                case CopulaType.FGM:
                    if (alpha < -1.0 || alpha > 1.0) throw new ArgumentException("Invalid alpha. Should be: -1 <= alpha <= +1.");
                    return;
                case CopulaType.Frank:
                    if (alpha == 0.0) throw new ArgumentException("Invalid alpha. Alpha cannot equal zero.");
                    return;
                case CopulaType.Gumbel:
                    if (alpha < 1.0) throw new ArgumentException("Invalid alpha. Should be: alpha >= 1.");
                    return;
                case CopulaType.Independent:
                    return; //no alpha for independent
                case CopulaType.Joe:
                    if (alpha < 0.0) throw new ArgumentException("Invalid alpha. Should be: alpha >= 0.");
                    return;
                default:
                    throw new ArgumentException("Copula type not expected.");
            }
        }

        private static double FrankPart(double d, double alpha)
        {
            return Math.Exp(-alpha * d) - 1;
        }
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.
