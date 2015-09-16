using System;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///
    ///References:
    ///Miller, Michael B. 2012. Mathematics and Statistics for Financial Risk Management. New York: John Wiley & Sons.
    ///Chapters 1 and 10
    
    public class GeometricSeries
    {
        /// <summary>
        /// Returns the sum of a geometric series of length n, who's first element is 1.
        /// For decay factor d, S = 1 + d + d^2 + ... + d^(n-1)
        /// </summary>
        /// <param name="decayFactor">Decay factor Typically between -1 adn +1.</param>
        /// <param name="length">Number of elements in the geometric series, must be positive.</param>
        /// <returns></returns>
        public static double SumOfGeometricSeries(double decayFactor, int length)
        {
            if (length < 1) return double.NaN;
            if (decayFactor == 1.0) return length;
            return (1.0 - Math.Pow(decayFactor, length)) / (1.0 - decayFactor);
        }

        /// <summary>
        /// Returns the inverse of the sum of a geometric series of length n, who's first element is 1.
        /// For decay factor d, S = 1 + d + d^2 + ... + d^(n-1). Return 1/S.
        /// </summary>
        /// <param name="decayFactor">Decay factor Typically between -1 adn +1.</param>
        /// <param name="length">Number of elements in the geometric series, must be positive.</param>
        /// <returns></returns>
        public static double InverseSumOfGeometricSeries(double decayFactor, int length)
        {
            if (length < 1) return double.NaN;
            if (decayFactor == 1.0) return 1.0 / length;
            return (1.0 - decayFactor) / (1.0 - Math.Pow(decayFactor, length));
        }

        /// <summary>
        /// Returns the sum of an infinite geometric series who's first element is 1.
        /// For decay factor d, S = 1 + d + d^2 + ...
        /// </summary>
        /// <param name="decayFactor">Decay factor. Typically between -1 adn +1.</param>
        /// <returns></returns>
        public static double SumOfInfiniteGeometricSeries(double decayFactor)
        {
            if (decayFactor >= 1) return double.PositiveInfinity;
            if (decayFactor <= 1) return double.NaN;
            return 1.0  / (1.0 - decayFactor);
        }

        /// <summary>
        /// Returns the half-life of a geometric series of length n, who's first element is 1.
        /// For decay factor d,  1 + d + d^2 + ... + d^(h-1) = 0.5 * [1 + d + d^2 + ... + d^(n-1)]
        /// </summary>
        /// <param name="decayFactor">Decay factor Typically between -1 adn +1.</param>
        /// <param name="length">Number of elements in the geometric series, must be positive.</param>
        /// <returns></returns>
        public static double HalfLifeOfGeometricSeries(double decayFactor, int length)
        {
            if (length < 1) return double.NaN;
            return Math.Log(0.5 + 0.5 * Math.Pow(decayFactor, length)) / Math.Log(decayFactor);
        }
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.