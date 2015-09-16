using System;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///
    ///References:
    ///Miller, Michael B. 2012. Mathematics and Statistics for Financial Risk Management. New York: John Wiley & Sons.
    ///Chapters 4, 6, 9
    
    public class MonteCarlo
    {
        public class RandomPlus : Random
        {
            public RandomPlus()
            {
                
            }

            public RandomPlus(int seed) : base(seed)
            {
                
            }

            /// <summary>
            /// Returns a random number from a standard normal distribution.
            /// </summary>
            public double NextStandardNormal()
            {
                return NextNormal(0.0, 1.0);
            }

            /// <summary>
            /// Returns a random number from a normal distribution with the given mean and standard deviation.
            /// </summary>
            public double NextNormal(double mean, double standardDeviation)
            {
                double U = NextDouble();
                double N = Distributions.NormalCumulativeDistributionFunctionInverse(U, mean, standardDeviation);
                return N;
            }

            /// <summary>
            /// Returns a random number from a Student's t distribution.
            /// </summary>
            public double NextStudentsT(double degreesOfFreedom)
            {
                double U = NextDouble();
                double S = Distributions.StudentsTCumulativeDistributionFunctionInverse(U, degreesOfFreedom);
                return S;
            }

            /// <summary>
            /// Returns a random number from a Poisson distribution with the given mean.
            /// </summary>
            public int NextPoisson(double mean)
            {
                double L = Math.Exp(-mean);
                int k = 0;
                double p = 1.0;
                do
                {
                    k++;
                    p *= NextDouble();
                } while (p > L);
                return k - 1;
            }
        }

        /// <summary>
        /// Correlation is not preserved when transforming two normal variables to two uniform variables.
        /// Returns the required correlation between the two normal variables needed to produce the desired correlation 
        /// between the two uniform variables.
        /// </summary>
        /// <param name="uniformCorr">The desired correlation between the two uniform variables</param>
        /// <returns></returns>
        public static double UniformToNormalCorrelation(double uniformCorr)
        {
            //Without correction max error would be ~1.81% and mean absolute error would be ~1.17%.
            //This approximation will reduce the max error to within ~0.06% and the mean absolute error to withing ~0.02%.
            double uCorrPi = uniformCorr * Math.PI;
            double normalCorr = uniformCorr + 0.01818 * Math.Sin(uCorrPi) - 0.00221 * Math.Sin(2 * uCorrPi) + 0.000647 * Math.Sin(3 * uCorrPi);
            return normalCorr;
        }
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.
