using System;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///
    ///A great reference for most of these formulas:
    ///Haug, Espen G. 1998. The Complete Guide to Option Pricing Formulas. New York: McGraw Hill.
    
    public class Options
    {
        public enum PutCallFlag
        {
            Put, Call
        }

        public enum OptionGreeks
        {
            Delta, Gamma, Vega, Theta, Rho
        }

        /// <summary>
        /// Gereralized Black, Scholes, Merton formula (1973) for European options with dividend
        /// </summary>
        /// <returns>price of the option</returns>
        public static double BlackScholesPrice(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag)
        {
            //setting dividendYield = 0 gives the classic Black Scholes model
            //setting dividendYield = foreign risk-free rate gives a model for European currency options, see Garman and Kohlhagen (1983)

            double sqrtT = Math.Sqrt(yearsToExpiry);
            double d1 = (Math.Log(underlyingPrice / strike) + (riskFreeRate - dividendYield + 0.5 * vol * vol) * yearsToExpiry) / (vol * sqrtT);
            double d2 = d1 - vol * sqrtT;
            if(putCallFlag == PutCallFlag.Call)
            {
                double N1 = Distributions.StandardNormalCumulativeDistributionFunction(d1);
                double N2 = Distributions.StandardNormalCumulativeDistributionFunction(d2);
                return N1 * underlyingPrice * Math.Exp(-dividendYield * yearsToExpiry) - N2 * strike * Math.Exp(-riskFreeRate * yearsToExpiry);
            }
            double Nn1 = Distributions.StandardNormalCumulativeDistributionFunction(-d1);
            double Nn2 = Distributions.StandardNormalCumulativeDistributionFunction(-d2);
            return Nn2 * strike * Math.Exp(-riskFreeRate * yearsToExpiry) - Nn1 * underlyingPrice * Math.Exp(-dividendYield * yearsToExpiry);
        }

        #region Black Scholes Greeks

        /// <summary>
        /// Black, Scholes, Merton formula (1973) for European option.
        /// Returns the price, and outputs the Greeks.
        /// Not as elegant looking, but faster than calculating each separately.
        /// </summary>
        public static double BlackScholesPriceAndGreeks(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag, out double delta, out double gamma, out double vega, out double theta, out double rho, out double convexity)
        {
            //setting dividendYield = 0 gives the classic Black Scholes model
            //setting dividendYield = foreign risk-free rate gives a model for European currency options, see Garman and Kohlhagen (1983)

            double sqrtT = Math.Sqrt(yearsToExpiry);
            double d1 = (Math.Log(underlyingPrice / strike) + (riskFreeRate - dividendYield + 0.5 * vol * vol) * yearsToExpiry) / (vol * sqrtT);
            double d2 = d1 - vol * sqrtT;
            double N1 = Distributions.StandardNormalCumulativeDistributionFunction(d1);
            double N2 = Distributions.StandardNormalCumulativeDistributionFunction(d2);
            double n1 = Distributions.StandardNormalProbabilityDensityFunction(d1);
            double n2 = Distributions.StandardNormalProbabilityDensityFunction(d2);

            double eNegRiskFreeRateTimesYearsToExpiry = Math.Exp(-riskFreeRate * yearsToExpiry);
            double eNegDivYieldYearsToExpiry = Math.Exp(-dividendYield * yearsToExpiry);

            double price;
            if(putCallFlag == PutCallFlag.Call)
                price = N1 * underlyingPrice * eNegDivYieldYearsToExpiry - N2 * strike * eNegRiskFreeRateTimesYearsToExpiry;
            else
                price = (1 - N2) * strike * eNegRiskFreeRateTimesYearsToExpiry - (1 - N1) * underlyingPrice * eNegDivYieldYearsToExpiry;

            if (putCallFlag == PutCallFlag.Call)
                delta = eNegDivYieldYearsToExpiry * N1;
            else
                delta = eNegDivYieldYearsToExpiry * (N1 - 1.0);

            gamma = eNegDivYieldYearsToExpiry * n1 / (underlyingPrice * vol * sqrtT);
            vega = underlyingPrice * Math.Exp(-dividendYield * yearsToExpiry) * n1 * sqrtT;

            if(putCallFlag == PutCallFlag.Call)
                rho = yearsToExpiry * strike * eNegRiskFreeRateTimesYearsToExpiry * N2;
            else
                rho = -yearsToExpiry * strike * eNegRiskFreeRateTimesYearsToExpiry * (1 - N2);

            double a = -underlyingPrice * eNegDivYieldYearsToExpiry * n1 * vol / (2.0 * sqrtT);
            double b = dividendYield * underlyingPrice * eNegDivYieldYearsToExpiry;
            if(putCallFlag == PutCallFlag.Call)
                theta = a + b * N1 - riskFreeRate * strike * eNegRiskFreeRateTimesYearsToExpiry * N2;
            else
                theta = a - b * (1 - N1) + riskFreeRate * strike * eNegRiskFreeRateTimesYearsToExpiry * (1 - N2);

            if (putCallFlag == PutCallFlag.Call)
                convexity = rho * ((n2 * sqrtT) / (N2 * vol) - yearsToExpiry);
            else
                convexity = -rho * ((n2 * sqrtT) / ((1 - N2) * vol) + yearsToExpiry);

            return price;
        }

        /// <summary>
        /// Black, Scholes, Merton formula (1973) for European option Greeks.
        /// </summary>
        public static double BlackScholesGreek(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag, OptionGreeks optionGreek)
        {
            switch (optionGreek)
            {
                case OptionGreeks.Delta:
                    return BlackScholesDelta(strike, underlyingPrice, yearsToExpiry, vol, riskFreeRate, dividendYield, putCallFlag);
                case OptionGreeks.Gamma:
                    return BlackScholesGamma(strike, underlyingPrice, yearsToExpiry, vol, riskFreeRate, dividendYield, putCallFlag);
                case OptionGreeks.Rho:
                    return BlackScholesRho(strike, underlyingPrice, yearsToExpiry, vol, riskFreeRate, dividendYield, putCallFlag);
                case OptionGreeks.Theta:
                    return BlackScholesTheta(strike, underlyingPrice, yearsToExpiry, vol, riskFreeRate, dividendYield, putCallFlag);
                case OptionGreeks.Vega:
                    return BlackScholesVega(strike, underlyingPrice, yearsToExpiry, vol, riskFreeRate, dividendYield, putCallFlag);
                default:
                    return double.NaN;
            }
        }

        /// <summary>
        /// Black, Scholes, Merton formula (1973) for European options with dividend. Delta = d[price]/d[underlyingPrice].
        /// </summary>
        /// <returns>option delta</returns>
        public static double BlackScholesDelta(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag)
        {
            double sqrtT = Math.Sqrt(yearsToExpiry);
            double d1 = (Math.Log(underlyingPrice / strike) + (riskFreeRate - dividendYield + 0.5 * vol * vol) * yearsToExpiry) / (vol * sqrtT);
            double N1 = Distributions.StandardNormalCumulativeDistributionFunction(d1);
            if (putCallFlag == PutCallFlag.Call)
                return Math.Exp(-dividendYield * yearsToExpiry) * N1;
            return Math.Exp(-dividendYield * yearsToExpiry) * (N1 - 1.0);
        }

        /// <summary>
        /// Black, Scholes, Merton formula (1973) for European options with dividend. Gamma = second derivative of price with respect to underlying price.
        /// </summary>
        /// <returns>option gamma</returns>
        public static double BlackScholesGamma(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag)
        {
            //Put and call values are the same. I've left putCallFlag in parameters so that all the Black Scholes formulas have the same argument list.
            double sqrtT = Math.Sqrt(yearsToExpiry);
            double d1 = (Math.Log(underlyingPrice / strike) + (riskFreeRate - dividendYield + 0.5 * vol * vol) * yearsToExpiry) / (vol * sqrtT);
            double n1 = Distributions.StandardNormalProbabilityDensityFunction(d1);
            return Math.Exp(-dividendYield * yearsToExpiry) * n1 / (underlyingPrice * vol * sqrtT);
        }

        /// <summary>
        /// Black, Scholes, Merton formula (1973) for European options with dividend. Vega = d[price]/d[vol].
        /// </summary>
        /// <returns>option vega</returns>
        public static double BlackScholesVega(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag)
        {
            //Put and call values are the same. I've left putCallFlag in parameters so that all the Black Scholes formulas have the same argument list.
            double sqrtT = Math.Sqrt(yearsToExpiry);
            double d1 = (Math.Log(underlyingPrice / strike) + (riskFreeRate - dividendYield + 0.5 * vol * vol) * yearsToExpiry) / (vol * sqrtT);
            double n1 = Distributions.StandardNormalProbabilityDensityFunction(d1);
            return underlyingPrice * Math.Exp(-dividendYield * yearsToExpiry) * n1 * sqrtT;
        }

        /// <summary>
        /// Black, Scholes, Merton formula (1973) for European options with dividend. Rho = d[price]/d[riskFreeRate].
        /// </summary>
        /// <returns>option rho</returns>
        public static double BlackScholesRho(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag)
        {
            double sqrtT = Math.Sqrt(yearsToExpiry);
            double d2 = (Math.Log(underlyingPrice / strike) + (riskFreeRate - dividendYield - 0.5 * vol * vol) * yearsToExpiry) / (vol * sqrtT);
            if (putCallFlag == PutCallFlag.Call)
            {
                double N2 = Distributions.StandardNormalCumulativeDistributionFunction(d2);
                return yearsToExpiry * strike * Math.Exp(-riskFreeRate * yearsToExpiry) * N2;
            }
            double Nn2 = Distributions.StandardNormalCumulativeDistributionFunction(-d2);
            return -yearsToExpiry * strike * Math.Exp(-riskFreeRate * yearsToExpiry) * Nn2;
        }

        /// <summary>
        /// Black, Scholes, Merton formula (1973) for European options with dividend. Theta = -d[price]/d[yearsToExpiry].
        /// </summary>
        /// <returns>option theta</returns>
        public static double BlackScholesTheta(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag)
        {
            double sqrtT = Math.Sqrt(yearsToExpiry);
            double d1 = (Math.Log(underlyingPrice / strike) + (riskFreeRate - dividendYield + 0.5 * vol * vol) * yearsToExpiry) / (vol * sqrtT);
            double d2 = d1 - vol * sqrtT;
            double a = -underlyingPrice * Math.Exp(-dividendYield * yearsToExpiry) * Distributions.StandardNormalProbabilityDensityFunction(d1) * vol / (2.0 * sqrtT);
            double b = dividendYield * underlyingPrice * Math.Exp(-dividendYield * yearsToExpiry);
            if (putCallFlag == PutCallFlag.Call)
                return a + b * Distributions.StandardNormalCumulativeDistributionFunction(d1) - riskFreeRate * strike * Math.Exp(-riskFreeRate * yearsToExpiry) * Distributions.StandardNormalCumulativeDistributionFunction(d2);
            return a - b * Distributions.StandardNormalCumulativeDistributionFunction(-d1) + riskFreeRate * strike * Math.Exp(-riskFreeRate * yearsToExpiry) * Distributions.StandardNormalCumulativeDistributionFunction(-d2);
        }

        /// <summary>
        /// Black, Scholes, Merton formula (1973) for European options with dividend. Convexity = second derivative of price with respect to interest rate.
        /// Not considered a standard Greek, but can be useful for portfolios with fixed income instruments and equity options.
        /// </summary>
        /// <returns>option interest rate convexity</returns>
        public static double BlackScholesRateConvexity(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag)
        {
            double sqrtT = Math.Sqrt(yearsToExpiry);
            double d2 = (Math.Log(underlyingPrice / strike) + (riskFreeRate - dividendYield - 0.5 * vol * vol) * yearsToExpiry) / (vol * sqrtT);
            double N2 = Distributions.StandardNormalCumulativeDistributionFunction(d2); 
            double n2 = Distributions.StandardNormalProbabilityDensityFunction(d2);

            double rho, convexity;
            if(putCallFlag == PutCallFlag.Call)
            {
                rho = yearsToExpiry * strike * Math.Exp(-riskFreeRate * yearsToExpiry) * N2;
                convexity = rho * ((n2 * sqrtT) / (N2 * vol) - yearsToExpiry);
            }
            else
            {
                rho = -yearsToExpiry * strike * Math.Exp(-riskFreeRate * yearsToExpiry) * (1 - N2);
                convexity = -rho * ((n2 * sqrtT) / ((1 - N2) * vol) + yearsToExpiry);
            }
            return convexity;
        }

        /// <summary>
        /// Calculates the implied vol for a European option using  Newton-Raphson approach.
        /// </summary>
        public static double BlackScholesImpliedVol(double price, double strike, double underlyingPrice, double yearsToExpiry, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag)
        {
            const double tolerance = 0.001;
            const int maxLoops = 16;

            double vol = Math.Sqrt(2 * Math.Abs(Math.Log(underlyingPrice / strike) / yearsToExpiry + riskFreeRate));    //Manaster and Koehler intial vol value
            vol = Math.Max(0.01, vol);
            double vega;
            double impliedPrice = BlackScholesPriceAndVega(strike, underlyingPrice, yearsToExpiry, vol, riskFreeRate, dividendYield, putCallFlag, out vega);

            int nLoops = 0;
            while (Math.Abs(impliedPrice - price) > tolerance)
            {
                nLoops++;
                if(nLoops > maxLoops)
                    throw new Exception("BlackScholesImpliedVol did not converge.");
                
                vol = vol - (impliedPrice - price) / vega;
                if(vol <= 0)
                    vol = 0.5 * (vol + (impliedPrice - price) / vega); //half way btwn previous estimate and zero
                impliedPrice = BlackScholesPriceAndVega(strike, underlyingPrice, yearsToExpiry, vol, riskFreeRate, dividendYield, putCallFlag, out vega);
            }
            return vol;
        }
        
        /// <summary>
        /// Used for BlackScholesImpliedVol
        /// </summary>
        /// <returns>price of option</returns>
        private static double BlackScholesPriceAndVega(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag, out double vega)
        {
            double sqrtT = Math.Sqrt(yearsToExpiry);
            double d1 = (Math.Log(underlyingPrice / strike) + (riskFreeRate - dividendYield + 0.5 * vol * vol) * yearsToExpiry) / (vol * sqrtT);
            double d2 = d1 - vol * sqrtT;
            if (putCallFlag == PutCallFlag.Call)
            {
                double N1 = Distributions.StandardNormalCumulativeDistributionFunction(d1);
                double N2 = Distributions.StandardNormalCumulativeDistributionFunction(d2);
                double nn1 = Distributions.StandardNormalProbabilityDensityFunction(d1);

                vega = underlyingPrice * Math.Exp(-dividendYield * yearsToExpiry) * nn1 * sqrtT;
                return N1 * underlyingPrice * Math.Exp(-dividendYield * yearsToExpiry) - N2 * strike * Math.Exp(-riskFreeRate * yearsToExpiry);
            }
            double Nn1 = Distributions.StandardNormalCumulativeDistributionFunction(-d1);
            double Nn2 = Distributions.StandardNormalCumulativeDistributionFunction(-d2);
            double n1 = Distributions.StandardNormalProbabilityDensityFunction(d1);

            vega = underlyingPrice * Math.Exp(-dividendYield * yearsToExpiry) * n1 * sqrtT;
            return Nn2 * strike * Math.Exp(-riskFreeRate * yearsToExpiry) - Nn1 * underlyingPrice * Math.Exp(-dividendYield * yearsToExpiry);
        }
        #endregion

        /// <summary>
        /// Barone Adesi and Whaley price approximation for an American option
        /// </summary>
        /// <returns>price of the option</returns>
        public static double BaroneAdesiWhaley(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag)
        {
            if (dividendYield <= 0 && putCallFlag == PutCallFlag.Call)
                return BlackScholesPrice(strike, underlyingPrice, yearsToExpiry, vol, riskFreeRate, dividendYield, putCallFlag);

            const double accuracy = 0.00001;
            const double maxIterations = 500;

            double volSqrd = vol * vol;
            double sqrtT = Math.Sqrt(yearsToExpiry);
            double m = 2.0 * riskFreeRate / volSqrd; 
            double n = 2.0 * (riskFreeRate - dividendYield) / volSqrd;
            double k = 1.0 - Math.Exp(-riskFreeRate * yearsToExpiry);
            double a = Math.Sqrt(Math.Pow((n - 1), 2.0) + (4.0 * m / k));
            double q;
            if(putCallFlag == PutCallFlag.Call)
                q = 0.5 * (1 - n + a);
            else
                q = 0.5 * (1 - n - a);

            //calculate seed value
            double undPriceSeed;
            if (putCallFlag == PutCallFlag.Call)
            {
                double b = -(n - 1.0) + Math.Sqrt(Math.Pow((n - 1), 2.0) + 4.0 * m);
                double undPriceStarInf = strike / (1.0 - 2.0 / b);
                double h = -((riskFreeRate - dividendYield) * yearsToExpiry + 2.0 * vol * sqrtT) * (strike / (undPriceStarInf - strike));
                undPriceSeed = strike + (undPriceStarInf - strike) * (1.0 - Math.Exp(h));
            }
            else
            {
                double b = -(n - 1.0) - Math.Sqrt(Math.Pow((n - 1), 2.0) + 4.0 * m);
                double undPriceStarInf = strike / (1.0 - 2.0 / b);
                double h = ((riskFreeRate - dividendYield) * yearsToExpiry - 2.0 * vol * sqrtT) * (strike / (strike - undPriceStarInf));
                undPriceSeed = undPriceStarInf + (strike - undPriceStarInf) * Math.Exp(h);
            }

            //Newton-Raphson
            int nIterations = 0; 
            double undPriceI = undPriceSeed;
            double g = 1.0;
            double gPrime = 1.0;
            
            while ((Math.Abs(g) > accuracy) && (Math.Abs(gPrime) > accuracy)  && (nIterations++ < maxIterations) && (undPriceI > 0.0))
            {
                double bsPrice = BlackScholesPrice(strike, undPriceI, yearsToExpiry, vol, riskFreeRate, dividendYield, putCallFlag);
                double d1 = (Math.Log(undPriceI / strike) + (riskFreeRate - dividendYield + 0.5 * volSqrd) * yearsToExpiry) / (vol * sqrtT);
                if (putCallFlag == PutCallFlag.Call)
                {
                    g = (1.0 - 1.0 / q) * undPriceI - strike - bsPrice + (1.0 / q) * undPriceI * Math.Exp(-dividendYield * yearsToExpiry) * Distributions.StandardNormalCumulativeDistributionFunction(d1);
                    gPrime = (1.0 - 1.0 / q) * (1.0 - Math.Exp(-dividendYield * yearsToExpiry) * Distributions.StandardNormalCumulativeDistributionFunction(d1))
                        + (1.0 / q) * Math.Exp(-dividendYield * yearsToExpiry) * Distributions.StandardNormalProbabilityDensityFunction(d1) * (1.0 / (vol * sqrtT));                    
                }
                else
                {
                    g = strike - undPriceI - bsPrice + (undPriceI / q) * (1.0 - Math.Exp(-dividendYield * yearsToExpiry));
                    gPrime = (1.0 / q - 1.0) * (1.0 - Math.Exp(-dividendYield * yearsToExpiry) * Distributions.StandardNormalCumulativeDistributionFunction(-d1))
                        + (1.0 / q) * Math.Exp(-dividendYield * yearsToExpiry) * Distributions.StandardNormalProbabilityDensityFunction(-d1) * (1.0 / (vol * sqrtT));                    
                }
                undPriceI = undPriceI - (g / gPrime);
            }

            double undPriceStar;
            if (Math.Abs(g) > accuracy)
            {
                //undPriceStar = undPriceSeed; //could continue with this value, currently throws exception
                throw new Exception("Newton-Raphson did not converge");
            }
            undPriceStar = undPriceI;

            double price;
            double bsPrice2 = BlackScholesPrice(strike, underlyingPrice, yearsToExpiry, vol, riskFreeRate, dividendYield, putCallFlag);
            if (putCallFlag == PutCallFlag.Call && underlyingPrice >= undPriceStar)
            {
                price = underlyingPrice - strike;
            }
            else if (putCallFlag == PutCallFlag.Put && underlyingPrice <= undPriceStar)
            {
                price = strike - underlyingPrice;
            }
            else
            {
                double d1 = (Math.Log(undPriceStar / strike) + (riskFreeRate - dividendYield + 0.5 * volSqrd) * yearsToExpiry) / (vol * sqrtT);
                if (putCallFlag == PutCallFlag.Put) d1 *= -1.0;
                double A = (1.0 - Math.Exp(-dividendYield * yearsToExpiry) * Distributions.StandardNormalCumulativeDistributionFunction(d1)) * (undPriceStar / q);
                if (putCallFlag == PutCallFlag.Put) A *= -1.0;
                price = bsPrice2 + A * Math.Pow((underlyingPrice / undPriceStar), q);
            }
            return Math.Max(price, bsPrice2); // know value will never be less than BS value
        }

        /// <summary>
        /// Bjerksund and Stensland price approximation for an American option
        /// </summary>
        /// <returns>price of the option</returns>
        public static double BjerksundStensland(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag)
        {            
            double b = riskFreeRate - dividendYield;
            if (putCallFlag == PutCallFlag.Call)
                return BjerksundStenslandCall(strike, underlyingPrice, yearsToExpiry, vol, riskFreeRate, b);
            return BjerksundStenslandCall(underlyingPrice, strike, yearsToExpiry, vol, riskFreeRate - b, -b);
        }


        /// <summary>
        /// Bjerksund and Stensland price approximation for an American call option
        /// </summary>
        /// <param name="x">strike</param>
        /// <param name="s">underlying price</param>
        /// <param name="t">time to expiry in years</param>
        /// <param name="v">vol</param>
        /// <param name="r">risk free rate</param>
        /// <param name="b">cost of carry = r - dividend yield</param>
        /// <returns>price of the option</returns>
        private static double BjerksundStenslandCall(double x, double s, double t, double v, double r, double b)
        {
            if (b >= r)
                return BlackScholesPrice(x, s, t, v, r, r - b, PutCallFlag.Call);

            double beta = 0.5 - b / (v * v) + Math.Sqrt(Math.Pow(b / (v * v) - 0.5, 2) + 2 * r / (v * v));
            double bInf = beta / (beta - 1.0) * x;
            double b0 = Math.Max(x, r / (r - b) * x);
            double ht = -(b * t + 2 * v * Math.Sqrt(t)) * b0 / (bInf - b0);
            double i = b0 + (bInf - b0) * (1 - Math.Exp(ht));
            double alpha = (i - x) * Math.Pow(i, -beta);
            
            if (s >= i)
                return s - x;

            double phi1 = BjerksundStenslandPhi(s, t, beta, i, i, r, b, v);
            double phi2 = BjerksundStenslandPhi(s, t, 1, i, i, r, b, v);
            double phi3 = BjerksundStenslandPhi(s, t, 1, x, i, r, b, v);
            double phi4 = BjerksundStenslandPhi(s, t, 0, i, i, r, b, v);
            double phi5 = BjerksundStenslandPhi(s, t, 0, x, i, r, b, v);
            if(alpha == 0)
                return phi2 - phi3 - x * (phi4 - phi5);
            return alpha * Math.Pow(s, beta) - alpha * phi1 + phi2 - phi3 - x * (phi4 - phi5);
        }

        /// <summary>
        /// Used in Bjerksund and Stensland formula
        /// </summary>
        private static double BjerksundStenslandPhi(double s, double t, double gamma, double h, double i, double r, double b, double v)
        {
            double lambda = (-r + gamma * b + 0.5 * gamma * (gamma - 1) * v * v) * t;
            double d1a = Math.Log(s / h);
            double d1b = (b + (gamma - 0.5) * v * v) * t;
            double d1 = -(d1a + d1b); //-(Math.Log(s / h) + (b + (gamma - 0.5) * v * v) * t) ;
            double d2 = (v * Math.Sqrt(t));
            double d = d1 / d2; //-(Math.Log(s / h) + (b + (gamma - 0.5) * v * v) * t) / (v * Math.Sqrt(t));
            double kappa = 2 * b / (v * v) + 2 * gamma - 1;
            double n1 = Distributions.StandardNormalCumulativeDistributionFunction(d);
            double n2 = Distributions.StandardNormalCumulativeDistributionFunction(d - 2.0 * Math.Log(i/s)/(v*Math.Sqrt(t)));
            double p1 = Math.Exp(lambda);
            double p2 = Math.Pow(s, gamma);
            double p3 = n2 == 0 ? n1 : n1 - Math.Pow(i / s, kappa) * n2;
            double phi = p1 * p2 * p3;
            return phi;
        }

        /// <summary>
        /// Cox-Ross-Rubinstein Binomial Tree
        /// </summary>
        /// <returns>price of the option</returns>
        public static double BinomialTree(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, double dividendYield, PutCallFlag putCallFlag, bool isEuropean, int nSteps)
        {
            double b = riskFreeRate - dividendYield;
            double[] optionValueArray = new double[nSteps + 1];
            
            int z;
            if (putCallFlag == PutCallFlag.Call)
                z = 1;
            else z = -1;

            double dt = yearsToExpiry / nSteps;
            double u = Math.Exp(vol * Math.Sqrt(dt));
            double d = 1.0 / u;
            double a = Math.Exp(b * dt);
            double p = (a - d) / (u - d);
            double df = Math.Exp(-riskFreeRate * dt);

            for(int i=0; i<=nSteps; i++)
                optionValueArray[i] = Math.Max(0, z * (underlyingPrice*Math.Pow(u, i)*Math.Pow(d, nSteps - i)-strike));

            for(int j=nSteps-1; j>=0; j--)
            {
                for(int i=0; i<=j; i++)
                {
                    double e = df *(p*optionValueArray[i+1] + (1-p)*optionValueArray[i]);
                    if(isEuropean)
                    {
                        optionValueArray[i] = e;
                    }
                    else
                    {
                        double intrinsic = z * (underlyingPrice * Math.Pow(u, i) * Math.Pow(d, j-i) - strike);
                        optionValueArray[i] = Math.Max(intrinsic, e);
                    }
                }
            }
            return optionValueArray[0];
        }

        
        /// <summary>
        /// Finite difference method for European or American option without dividend.
        /// </summary>
        public static double FiniteDifference(double strike, double underlyingPrice, double yearsToExpiry, double vol, double riskFreeRate, PutCallFlag putCallFlag, bool isEuropean, int nPriceSteps, out double currentDelta, out double currentGamma, out double currentTheta)
        {
            //Adapted from Wilmott, Paul. Paul Wilmott Introduces Quantitative Finance. Wiley, 2001.

            double[] undPriceArray = new double[nPriceSteps + 1]; //underlying asset prices
            double[] payoffArray = new double[nPriceSteps + 1];  //payoff array

            double dS = 2 * underlyingPrice / nPriceSteps; //max price = twice the strike or underlying
            double dt = 0.9 / (vol * vol) / (nPriceSteps * nPriceSteps); 
            int nTimeSteps = (int)(yearsToExpiry / dt) + 1; 
            dt = yearsToExpiry / nTimeSteps;

            double[,] optionValueArray = new double[nPriceSteps + 1, nTimeSteps + 1]; //option value array

            int q = 1;
            if(putCallFlag == PutCallFlag.Put) q = -1;

            for(int i=0; i<=nPriceSteps; i++)
            {
                undPriceArray[i] = i * dS; 
                optionValueArray[i, 0] = Math.Max(q * (undPriceArray[i] - strike), 0); 
                payoffArray[i] = optionValueArray[i, 0]; 
            }
            
            for(int k = 1; k<= nTimeSteps; k++)
            {
                for(int i = 1; i< nPriceSteps; i++)  
                {
                    double delta = (optionValueArray[i + 1, k - 1] - optionValueArray[i - 1, k - 1]) / 2.0 / dS; 
                    double gamma = (optionValueArray[i + 1, k - 1] - 2 * optionValueArray[i, k - 1] + optionValueArray[i - 1, k - 1]) / dS / dS; 
                    double theta = -0.5 * vol * vol * undPriceArray[i] * undPriceArray[i] * gamma - riskFreeRate * undPriceArray[i] * delta + riskFreeRate * optionValueArray[i, k - 1]; 
                    optionValueArray[i, k] = optionValueArray[i, k - 1] - dt * theta;
                }
                 
                optionValueArray[0, k] = optionValueArray[0, k - 1] * (1 - riskFreeRate * dt); //Boundary condition at S=0
                optionValueArray[nPriceSteps, k] = 2 * optionValueArray[nPriceSteps - 1, k] - optionValueArray[nPriceSteps - 2, k]; //Boundary condition at S=infinity

                //Check for early exercise
                if(!isEuropean)
                { 
                    for(int i=0; i<=nPriceSteps; i++)
                        optionValueArray[i, k] = Math.Max(optionValueArray[i, k], payoffArray[i]);
                }

            }

            int underlyingIndex = (int)(underlyingPrice / dS);
            double currentValue = optionValueArray[underlyingIndex, nTimeSteps];
            currentDelta = (optionValueArray[underlyingIndex + 1, nTimeSteps] - optionValueArray[underlyingIndex - 1, nTimeSteps]) / 2.0 / dS; //Central difference
            currentGamma = (optionValueArray[underlyingIndex + 1, nTimeSteps] - 2 * optionValueArray[underlyingIndex, nTimeSteps] + optionValueArray[underlyingIndex - 1, nTimeSteps]) / dS / dS; //Central difference
            currentTheta = -0.5 * vol * vol * undPriceArray[underlyingIndex] * undPriceArray[underlyingIndex] * currentGamma - riskFreeRate * undPriceArray[underlyingIndex] * currentDelta + riskFreeRate * optionValueArray[underlyingIndex, nTimeSteps]; //Black-Scholes

            return currentValue;
        }

        
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.