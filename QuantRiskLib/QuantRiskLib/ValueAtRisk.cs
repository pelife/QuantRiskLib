using System;
using System.Collections.Generic;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///
    ///References:
    ///Miller, Michael B. 2012. Mathematics and Statistics for Financial Risk Management. New York: John Wiley & Sons.
    ///Chapters 5, and 10
    
    public partial class ValueAtRisk
    {
        /// <summary>
        /// Assumes returns are normally distributed.
        /// WARNING: Not very accurate. Intended for educational use and as a baseline for comparison to other models.
        /// </summary>
        /// <param name="returnArray">Historical returns from which VaR is to be calculated. The last value, with the highest index is assumed to be the most recent data point.</param>
        /// <param name="windowLength">Length of the VaR window. The number of historical returns that will be used to calculate the VaR.</param>
        /// <param name="confidenceLevel">VaR confidence level. 95% and 99% are typical values. If confidenceLevel is 95% then we expect 95% of days to be better (have a more positive return) than the VaR.</param>
        public static double NormalValueAtRisk(double[] returnArray, int windowLength, double confidenceLevel)
        {
            double[] da = new double[windowLength];
            for(int i=0; i<windowLength; i++)
                da[i] = returnArray[returnArray.Length - windowLength + i];
            double mean = Moments.Mean(da);
            double standardDeviation = Moments.StandardDeviation(da);
            double var = Distributions.NormalCumulativeDistributionFunctionInverse(1 - confidenceLevel, mean, standardDeviation);
            return -var; //By convention VaR is quoted as a positive value, even though it corresponds to a loss.
        }


        /// <summary>
        /// Nonparametric method. Distribution of expected returns is assumed to correspond to historical returns.
        /// WARNING: Not very accurate. Intended for educational use and as a baseline for comparisons.
        /// </summary>
        /// <param name="returnArray">Historical returns from which VaR is to be calculated. The last value, with the highest index is assumed to be the most recent data point.</param>
        /// <param name="windowLength">Length of the VaR window. The number of historical returns that will be used to calculate the VaR.</param>
        /// <param name="confidenceLevel">VaR confidence level. 95% and 99% are typical values. If confidenceLevel is 95% then we expect 95% of days to be better (have a more positive return) than the VaR.</param>
        public static double HistoricalValueAtRisk(double[] returnArray, int windowLength, double confidenceLevel)
        {
            return HybridValueAtRisk(returnArray, windowLength, confidenceLevel, 1.00);
        }

        /// <summary>
        /// Hybrid VaR is based on historical data, but weights more recent data more heavily.
        /// </summary>
        /// <param name="returnArray">Historical returns from which VaR is to be calculated. The last value, with the highest index is assumed to be the most recent data point.</param>
        /// <param name="windowLength">Length of the VaR window. The number of historical returns that will be used to calculate the VaR.</param>
        /// <param name="confidenceLevel">VaR confidence level. 95% and 99% are typical values. If confidenceLevel is 95% then we expect 95% of days to be better (have a more positive return) than the VaR.</param>
        /// <param name="decayFactor">For a decay factor d, the most recent data point has weight 1, the second d, the third d^2, ...</param>
        public static double HybridValueAtRisk(double[] returnArray, int windowLength, double confidenceLevel, double decayFactor)
        {
            if (returnArray.Length < windowLength)
                throw new ArgumentException(string.Format("returnArray does not have enough data points to calcualte VaR for nDay = {0}.", windowLength));

            if (confidenceLevel < 0.50)
            {
                double[] negReturnArray = new double[returnArray.Length];
                for (int j = 0; j < returnArray.Length; j++)
                    negReturnArray[j] = -returnArray[j];
                return -HybridValueAtRisk(negReturnArray, windowLength, 1.0 - confidenceLevel, decayFactor);
            }

            //initialize VaR to the minimum return
            double var = double.MaxValue;
            for (int j = returnArray.Length - windowLength; j < returnArray.Length; j++)
                if (returnArray[j] < var)
                    var = returnArray[j];

            //Loop through all data to find the minimum return, then the second worst, ...
            //This is generally faster than sorting when confidenceLevel is high, becaus we only need to identify a few minimums.
            double totalWt = GeometricSeries.SumOfGeometricSeries(decayFactor, windowLength);
            double varLimitWt = (1 - confidenceLevel) * totalWt;
            double cummWt = 0;
            List<int> alreadyHaveMinIndexList = new List<int>();
            for (int i = 0; i < windowLength; i++)
            {
                double min = double.MaxValue;
                int minIndex = -1;
                for (int j = returnArray.Length - windowLength; j < returnArray.Length; j++)
                {
                    if (returnArray[j] >= var && returnArray[j] < min && alreadyHaveMinIndexList.Contains(j) == false)
                    {
                        min = returnArray[j];
                        minIndex = j;
                    }
                }
                alreadyHaveMinIndexList.Add(minIndex);

                double wt = Math.Pow(decayFactor, returnArray.Length - minIndex - 1);
                cummWt += wt;
                if (cummWt >= varLimitWt) break;
                var = min;
            }
            return -var;    //By convention VaR is quoted as a positive value, even though it corresponds to a loss.
        }

        #region Incremental VaR
        /// <summary>
        /// Returns the incremental incremental VaR for the subportfolio relative to the portfolio.
        /// </summary>
        /// <param name="portfolioArray">Return array of portfolio. The last element (index = n-1), is the most recent.</param>
        /// <param name="subPortfolioArray">Return array of the sub-portfolio for which incremental VaR is being measured. The last element (index = n-1), is the most recent.</param>
        /// <param name="decayFactor">In most applications, the decay factor is between 0 and 1. Weigth on the last element in arrays is 1.0, the 2nd to last element d, 3rd to last d^2, ...</param>
        /// <param name="portfolioValutAtRisk">Value at risk of the portfolio.</param>
        /// <param name="length">Window length. Method uses the most recent n points, n = length.</param>
        public static double IncrementalValueAtRisk(double[] portfolioArray, double[] subPortfolioArray, double decayFactor, double portfolioValutAtRisk, int length)
        {
            double[] portfolioArrayLen = Tools.MostRecentValues(portfolioArray, length);
            double[] subPortfolioArrayLen = Tools.MostRecentValues(subPortfolioArray, length);
            double portfolioStandardDeviation = Moments.StandardDeviation(portfolioArray, decayFactor, length);

            if (Tools.ArrayAllEqual(subPortfolioArray)) return 0.0;

            double iStandardDeviation = Moments.IncrementalStandardDeviation(portfolioArrayLen, subPortfolioArrayLen, decayFactor);
            return (portfolioValutAtRisk / portfolioStandardDeviation) * iStandardDeviation;
        }
        #endregion

    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.