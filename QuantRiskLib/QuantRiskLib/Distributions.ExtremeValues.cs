using System;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    
    public partial class Distributions
    {
        /// <summary>
        /// For n independent draws from a standard uniform distribution (range = 0-1), the probability that the minimum is less than or equal to x.
        /// </summary>
        public static double IndependentStandardUniformMinimumCumulativeDistributionFunction(double x, int nDistributions)
        {
            double p = 0;
            int sign = 1;
            for (int i = 1; i < nDistributions; i++)
            {
                p += MMath.Combinations(nDistributions, i) * Math.Pow(x, i) * sign;
                sign *= -1;
            }
            return p;
        }

        /// <summary>
        /// For n independent draws from a standard uniform distribution (range = 0-1), the probability that the maximum is less than or equal to x.
        /// </summary>
        public static double IndependentStandardUniformMaximumCumulativeDistributionFunction(double x, int nDistributions)
        {
            return Math.Pow(x, nDistributions);
        }
    }
}


//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.