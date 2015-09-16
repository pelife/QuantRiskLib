using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///
    ///References:
    ///Miller, Michael B. 2012. Mathematics and Statistics for Financial Risk Management. New York: John Wiley & Sons.
    ///Chapter 1
    
    public class MMath
    {
        #region Permutations
        /// <summary>
        /// Returns the number of possible permutations of k objects from a set of n object. The order of the chosen objects does matter.
        /// </summary>
        /// <param name="n">Number of objects</param>
        /// <param name="k">Number of objects chosen</param>
        /// <returns></returns>
        public static double Permutations(int n, int k)
        {
            if (n < 0 || k < 0) throw new ArgumentException("Combin not defined for negative n or k");
            //long permut = 1;
            //for (int i = n; i > (n - k); i--)
            //    permut = checked(permut * i);\
            double logPermut = LogPermutations(n, k);
            return Math.Exp(logPermut);
        }

        internal static double LogPermutations(int n, int k)
        {
            if (n < 0 || k < 0) throw new ArgumentException("Combin not defined for negative n or k");
            double permut = 0;
            for (int i = n; i > (n - k); i--)
                permut += Math.Log(i);
            return permut;
        }

        /// <summary>
        /// Use this class to generate all possible permutations of elements from an array of doubles.
        /// </summary>
        public class PermutationGenerator
        {
            private readonly double[] _data;

            /// <summary>
            /// Initialize the class with the array of doubles for which you wish to make permutations.
            /// </summary>
            /// <param name="data">The array of doubles for which you wish to make permutations.</param>
            public PermutationGenerator(double[] data)
            {
                _data = data;
            }

            /// <summary>
            /// Generate all possible permutations of length nChosen.
            /// </summary>
            /// <param name="nChosen">nChosen must be ≥ 0 and ≤ data.Length</param>
            public List<double[]> GetPermutations(int nChosen)
            {
                if (nChosen < 0)
                    throw new ArgumentException("nChosen cannot be less than zero. nChosen must be ≥ 0 and ≤ data.Length");
                if (nChosen > _data.Length)
                    throw new ArgumentException("nChosen cannot be less than zero. nChosen must be ≥ 0 and ≤ data.Length");
                if (nChosen == 0)
                    return new List<double[]> { new double[0] };

                if (nChosen == _data.Length)
                    return GetPermut(_data);

                CombinationGenerator cg = new CombinationGenerator(_data);
                List<double[]> combins = cg.Combinations(nChosen);
                List<double[]> permuts = new List<double[]>();
                foreach (double[] da in combins)
                {
                    permuts.AddRange(GetPermut(da));
                }
                return permuts;
            }

            /// <summary>
            /// This iterative algorithm generates all of the permutations of the array da[] (i.e. all of the permutation of da[].Length).
            /// </summary>
            private static List<double[]> GetPermut(double[] da)
            {
                if (da.Length == 1)
                    return new List<double[]> { da };

                List<double[]> permuts = new List<double[]>();
                for (int i = 0; i < da.Length; i++)
                {
                    List<double[]> subPermuts = GetPermut(ArrayMinusElementAtX(da, i));
                    foreach (double[] subPermut in subPermuts)
                    {
                        double[] permut = AddValueToStartOfArray(subPermut, da[i]);
                        permuts.Add(permut);
                    }
                }

                return permuts;
            }

            /// <summary>
            /// Returns a copy of the array da[], minus the element at index x. 
            /// </summary>
            private static double[] ArrayMinusElementAtX(double[] da, int x)
            {
                double[] oa = new double[da.Length - 1];
                int j = 0;
                for (int i = 0; i < da.Length; i++)
                {
                    if (i == x) continue;
                    oa[j] = da[i];
                    j++;
                }
                return oa;
            }

            /// <summary>
            /// Returns a copy of da[] with d added to the start.
            /// </summary>
            private static double[] AddValueToStartOfArray(double[] da, double d)
            {
                double[] oa = new double[da.Length + 1];
                oa[0] = d;
                for (int i = 1; i < oa.Length; i++)
                    oa[i] = da[i - 1];
                return oa;
            }
        }
        #endregion

        #region Combinations
        /// <summary>
        /// Returns the number of possible combinations of k objects from a set of n object. The order of the chosen objects does not matter.
        /// </summary>
        /// <param name="n">Number of objects</param>
        /// <param name="k">Number of objects chosen</param>
        public static double Combinations(int n, int k)
        {
            double logCombin = LogCombin(n, k);
            double combin = Math.Exp(logCombin);
            return Math.Round(combin, 0);
        }

        internal static double LogCombin(int n, int k)
        {
            return LogFactorial(n) - LogFactorial(k) - LogFactorial(n - k);
        }

        private static double LogFactorial(int x)
        {
            double lf = 0;
            int i = 2;

            if(x >= 100000)
            {
                lf = 1051299.2218991187;
                i = 100001;
            }
            else if(x >= 50000)
            {
                lf = 490995.24304985348;
                i = 50001;
            }
            else if(x >= 10000)
            {
                lf = 82108.927836814153;
                i = 10001;
            }
            else if(x >= 5000)
            {
                lf = 37591.143508876841;
                i = 5001;
            }
            else if(x >= 4000)
            {
                lf = 29181.264544594731;
                i = 4001;
            }
            else if(x >= 3000)
            {
                lf = 21024.024853045572;
                i = 3001;
            }
            else if(x >= 2000)
            {
                lf = 13206.52435051381;
                i = 2001;
            }
            else if(x >= 1000)
            {
                lf = 5912.1281784881712;
                i = 1001;
            }
            else if (x >= 500)
            {
                lf = 2611.3304584601597;
                i = 501;
            }
            else if(x >= 100)
            {
                lf = 363.73937555556358;
                i = 101;
            }

            for (; i <= x; i++)
                lf += Math.Log(i);
            return lf;
        }

        /// <summary>
        /// Use this class to generate all possible combinations of elements of an array of doubles.
        /// </summary>
        public class CombinationGenerator
        {
            private readonly double[] _data;

            /// <summary>
            /// Initialize the class with the array of doubles for which you wish to make combinations.
            /// </summary>
            /// <param name="data">The array of doubles for which you wish to make combinations.</param>
            public CombinationGenerator(double[] data)
            {
                if (data.Length > 30)
                    throw new ArgumentException("data.Length must be ≤ 30"); //prevents overflow in Combinations(...)
                _data = data;
            }

            /// <summary>
            /// Generate all possible combinations of length nChosen.
            /// </summary>
            /// <param name="nChosen">nChosen must be ≥ 0 and ≤ data.Length</param>
            public List<double[]> Combinations(int nChosen)
            {
                if (nChosen < 0)
                    throw new ArgumentException("nChosen cannot be less than zero. nChosen must be ≥ 0 and ≤ data.Length");
                if (nChosen > _data.Length)
                    throw new ArgumentException("nChosen cannot be less than zero. nChosen must be ≥ 0 and ≤ data.Length");
                if (nChosen == 0)
                    return new List<double[]> { new double[0] };
                if (nChosen == _data.Length)
                    return new List<double[]> { _data };

                //This algorithm may not be super efficient, but it is very easy to follow.
                //It runs through all of the binary numbers from 0 to (2^n - 1) then uses the bit arrays where the number of 1's is equal
                //to nChosen. These bit arrays correspond to all of the possible combinations of the data array.
                List<double[]> combinations = new List<double[]>();
                int twoToTheN = (int)Math.Pow(2, _data.Length);
                for (int i = 0; i < twoToTheN; i++)
                {
                    BitArray ba = new BitArray(new[] { i });
                    if (CountOnes(ba) != nChosen) continue;

                    double[] combin = new double[nChosen];
                    int combinIndex = 0;
                    for (int j = 0; j < _data.Length; j++)
                    {
                        if (ba[j])
                        {
                            combin[combinIndex] = _data[j];
                            combinIndex++;
                        }
                    }
                    combinations.Add(combin);
                }
                return combinations;
            }

            /// <summary>
            /// Counts the number of bytes equal to 1 in a bit array.
            /// </summary>
            private static int CountOnes(BitArray ba)
            {
                int c = 0;
                for (int i = 0; i < ba.Length; i++)
                {
                    if (ba[i])
                        c++;
                }
                return c;
            }

        }
        #endregion


        /// <summary>
        /// Returns n! 
        /// 0! = 1; otherwise, n! = n * (n-1) * (n-2) * ... * 2 * 1,
        /// </summary>
        public static long Factorial(int n)
        {
            if (n < 0) throw new ArgumentException("Factorial not defined for negative n");
            if (n > 20) throw new ArgumentException("Answer will exceed max long");
            long fact = 1;
            for (int i = n; i > 0; i--)
                fact *= i;
            return fact;
        }

        #region Gamma and Beta functions
        /// <summary>
        /// Integral from 0 to infinity of e^(-t) * t^(z-1) dt
        /// </summary>
        /// <param name="z">If z > 143 the return value will exceed the double.MaxValue. The function will throw an exception.</param>
        public static double GammaFunction(double z)
        {
            //approximation based on Stirling's formula
            //based on approximation from NIST Handbook of Mathematical Function. 2010. Cambridge University Press.

            if (z > 143)
                throw new ArgumentException("Cannot currently compute gamma function for z > 143");

            //need z<21 or z<11 becaus Factorial(x) will not work for x>20
            if (z > 0 && z < 21 && z % 1.0 == 0)    
                return Factorial((int)Math.Round(z - 1, 0));
            if (z > 0 && z < 11 && z % 0.5 == 0)
            {
                int n = (int)z;
                const double sqrtPi = 1.77245385090552;
                return sqrtPi * Factorial(2 * n) / (Math.Pow(4, n) * Factorial(n));
            }

            const double sqrtTwoPi = 2.5066282746310002;
            double d = 1.0 + 1.0 / (12 * z) + 1.0 / (288 * z * z) - 139.0 / (51840 * Math.Pow(z, 3)) - 571.0 / (2488320 * Math.Pow(z, 4))
                + 163879.0 / (209018880 * Math.Pow(z, 5)) + 5246819.0 / (75246796800 * Math.Pow(z, 6));
            double g = Math.Pow(z, z - 0.5) * Math.Exp(-z) * sqrtTwoPi * d;
            return g;
        }

        /// <summary>
        /// Integral from 0 to infinity of e^(-t) * t^(z-1) dt
        /// </summary>
        /// <param name="z">If Real[z] > 143 the return value will exceed the double.MaxValue. The function will throw an exception.</param>
        public static Complex GammaFunction(Complex z)
        {
            if (z.Imaginary == 0)
                return GammaFunction(z.Real);
            if (z.Real > 143)
                throw new ArgumentException("Cannot currently compute gamma function for Real[z] > 143");

            Complex logGamma = LogGammaFunction(z);
            return Complex.Exp(logGamma);
        }

        /// <summary>
        /// Returns the log of the gamma function.
        /// </summary>
        public static double LogGammaFunction(double z)
        {
            //approximation based on Stirling's formula
            //based on approximation from NIST Handbook of Mathematical Function. 2010. Cambridge University Press.

            if(z > 0 && z % 1.0 == 0)
                return LogFactorial((int) Math.Round(z - 1, 0));
            if (z > 0 && z % 0.5 == 0)
            {
                int n = (int)z;
                const double halfLogPi = 0.5723649429247;
                const double ln4 = 1.38629436111989;
                return halfLogPi + LogFactorial(2 * n) - n * ln4 - LogFactorial(n);
            }

            const double halfLogTwoPi = 0.918938533204673; //0.5 * ln(2*pi)
            double lg = (z - 0.5) * Math.Log(z) - z + halfLogTwoPi + 1.0 / (12 * z) - 1.0 / (360 * Math.Pow(z, 3)) + 1.0 / (1260 * Math.Pow(z, 5));
            return lg;
        }

        /// <summary>
        /// Returns the log of the gamma function.
        /// </summary>
        public static Complex LogGammaFunction(Complex z)
        {
            if (z.Imaginary == 0)
                return LogGammaFunction(z.Real);

            const double halfLogTwoPi = 0.918938533204673; //0.5 * ln(2*pi)
            Complex lg = (z - 0.5) * Complex.Log(z) - z + halfLogTwoPi + 1.0 / (12 * z) - 1.0 / (360 * Complex.Pow(z, 3)) + 1.0 / (1260 * Complex.Pow(z, 5));
            return lg;
        }

        /// <summary>
        /// B(a, b) = Integral(0,1) of (t^(a-1))*((1-t)^(b-1)) dt. Also know as Euler integral. 
        /// </summary>
        public static double BetaFunction(double a, double b)
        {
            if (a + b > 143)
            {
                if (b > 20)
                    return 2.5066282746310002 * Math.Pow(a, a - 0.5) * Math.Pow(b, b - 0.5) / Math.Pow(a + b, a + b - 0.5);
                return GammaFunction(b) * Math.Pow(a, -b);
            }

            return GammaFunction(a) * GammaFunction(b) / GammaFunction(a + b);
        }

        /// <summary>
        /// B(x, a, b) = Integral(0,x) of (t^(a-1))*((1-t)^(b-1)) dt. Also know as Euler integral. 
        /// </summary>
        public static double IncompleteBetaFunction(double x, double a, double b)
        {
            if (x == 0) return 0;
            if (x == 1) return BetaFunction(a, b);

            //The numerical integration is less accurate for high values of x, so variable step sizes are used with 
            //progressively more steps for higher values of x.
            if (x <= 0.9)
            {
                return PartialIncompleteBetaFunction(0, x, a, b);
            }
            if (x <= 0.99)
            {
                double s1 = PartialIncompleteBetaFunction(0, 0.9, a, b);
                double s2 = PartialIncompleteBetaFunction(0.90, x, a, b);
                return s1 + s2;
            }
            else
            {
                double s1 = PartialIncompleteBetaFunction(0, 0.9, a, b);
                double s2 = PartialIncompleteBetaFunction(0.9, 0.99, a, b);
                double s3 = PartialIncompleteBetaFunction(0.99, x, a, b);
                return s1 + s2 + s3;
            }
        }

        private static double PartialIncompleteBetaFunction(double xL, double xU, double a, double b)
        {
            if (xU == xL) return 0.0;
            if (xU < xL) return Double.NaN;

            int nSteps = 1000;
            if (a < 1) nSteps = 80000;

            double dt = (xU - xL) / nSteps;
            double t = xL + 0.5 * dt;
            double sum = 0;
            for (int i = 0; i < nSteps; i++)
            {
                sum += Math.Pow(t, a - 1) * Math.Pow(1.0 - t, b - 1) * dt;
                t += dt;
            }
            return sum;
        }

        /// <summary>
        /// IncompleteBetaFunction(x, a, b) / BetaFunction(a, b).
        /// Also known as the CDF of the beta distribution.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        public static double RegularizedIncompleteBetaFunction(double x, double a, double b)
        {
            //parts of this function are based on:
            //A Method for Computing the Incomplete Beta Function by A.R. DiDonato and M.P. Jarnagin, Jr., 1996.

            if (x <= 0.0) return 0.0; 
            if (x >= 1.0) return 1.0;

            if(a % 1.0 == 0 && b % 1.0 == 0 && a + b > 0)
                return RegularizedIncompleteBetaFunction(x, (int) a, (int) b);
            if (b % 1.0 == 0 && a + b < 172)
                return RegularizedIncompleteBetaFunction(x, a, (int)b);
            if (a % 1.0 == 0 && a + b < 172)
                return 1.0 - RegularizedIncompleteBetaFunction(1.0 - x, b, (int)a);
            if (a == 0.5 && b == 0.5)
                return (2.0 / Math.PI) * Math.Atan(Math.Sqrt(x / (1.0 - x)));
            if (a == 0.5 && b % 0.5 == 0)
            {
                return 1.0 - RegularizedIncompleteBetaFunction(1.0 - x, b, 0.5);
            }
            if (a % 0.5 == 0 && b == 0.5)
            {
                if (a < 45)
                {
                    //CASE B (partial)
                    double s = 0;
                    int n = (int) Math.Round(a - 0.5, 0);
                    double logGammaHalf = LogGammaFunction(0.5);
                    for (int i = 0; i < n; i++)
                    {
                        s += Math.Exp(LogGammaFunction(i + 1) - LogGammaFunction(i + 1.5) - logGammaHalf)*Math.Pow(x, i);
                    }
                    return RegularizedIncompleteBetaFunction(x, 0.5, 0.5) - Math.Sqrt(x*(1 - x))*s;
                }
                else
                {
                    //CASE C (partial)
                    double M = BetaFunction(a, 0.5);
                    double lambda = Math.Sqrt(1.0 - Math.Pow(M*Math.Sqrt((a - 1.0)/Math.PI)*Double.Epsilon, 1.0/(a - 1.0)));

                    //Gaussian weights
                    double[] weights = new[]
                        {
                            0.066671344308688, 0.14945134915058, 0.21908636251598, 0.26926671931000, 0.29552422471475,
                            0.29552422471475, 0.26926671931000, 0.21908636251598, 0.14945134915058, 0.066671344308688
                        }; 
                    
                    //mgas = (1+ y)/2 where y is the Gaussian abscissa
                    double[] mgas = new[]
                        {
                            0.013046735791414, 0.067468316655507, 0.16029521585049, 0.28330230293538, 0.42556283050918,
                            0.57443716949081, 0.71669769706462, 0.83970478414951, 0.93253168334449, 0.98695326420859
                        };    

                    double sqrtOneMinusX = Math.Sqrt(1 - x);
                    double s = 0;
                    for (int i = 0; i < 10; i++)
                    {
                        s += weights[i] * Math.Pow(1.0 - Math.Pow((lambda - sqrtOneMinusX)*mgas[i] + sqrtOneMinusX, 2), a - 1.0);
                    }
                    double l = Math.Exp(LogGammaFunction(a + 0.5) - LogGammaFunction(a) - LogGammaFunction(0.5));
                    return (lambda - sqrtOneMinusX)*l*s;
                }
            }
            if (a % 0.5 == 0 && b % 0.5 == 0)
            {
                //CASE B and C
                double s = 0;
                int n = (int)Math.Round(b - 0.5, 0);
                double logGammaA = LogGammaFunction(a);
                double xToTheA = Math.Pow(x, a);
                for (int i = 0; i < n; i++)
                {
                    s += Math.Exp(LogGammaFunction(a + i + 0.5) - logGammaA - LogGammaFunction(i + 1.5)) * Math.Pow(1 - x, i) * xToTheA;
                }
                double d = RegularizedIncompleteBetaFunction(x, a, 0.5) + Math.Sqrt(1 - x) * s;
                return Math.Max(0.0, Math.Min(1.0, d));
            }
            if (x > 0.5)
                return 1.0 - RegularizedIncompleteBetaFunction(1.0 - x, b, a); 
   
            double ribf = IncompleteBetaFunction(x, a, b) / BetaFunction(a, b);
            return Math.Max(0.0, Math.Min(1.0, ribf));
        }

        private static double RegularizedIncompleteBetaFunction(double x, int a, int b)
        {
            double s = 0;
            int c = a + b - 1;
            if(c < 21)
            {
                for(int i = a; i < a + b; i++)
                {
                    s += Math.Pow(x, i) * Math.Pow(1.0 - x, c - i) / (Factorial(i) * Factorial(c - i));
                }
                s *= Factorial(c);
            }
            else
            {
                for (int i = a; i < a + b; i++)
                {
                    s += Math.Pow(x, i) * Math.Pow(1.0 - x, c - i)  * Combinations(c, i);
                }
            }

            return Math.Max(0.0, Math.Min(1.0, s));
        }

        private static double RegularizedIncompleteBetaFunction(double x, double a, int b)
        {
            if(a + b > 172)
                throw new ArgumentException("Cannot currently compute RegularizedIncompleteBetaFunction for a + b > 172");

            double s = 0;
            for (int i = 1; i < b + 1; i++)
            {
                s += Math.Pow(1 - x, i - 1) * Math.Exp(LogGammaFunction(a + i - 1) - LogGammaFunction(i)); //Math.Pow(1 - x, i - 1) * GammaFunction(a + i - 1) / GammaFunction(i)
            }
            s *= (Math.Pow(x, a) * Math.Exp(-LogGammaFunction(a)));   //s *= (Math.Pow(x, a) / GammaFunction(a))

            return s;
        }

        /// <summary>
        /// Returns x(x+1)(x+2)...(x+n-1) for n>0; for n=0 returns 1.
        /// </summary>
        public static double PochhammerFunction(double x, double n)
        {
            if (n < 0) return Double.NaN;
            if (n == 0) return 1;
            double p = 1;
            for (int i = 0; i < n; i++)
                p *= (x + i);
            return p;
        }
        #endregion

        /// <summary>
        /// Returns the real roots of a quadratic equation.
        /// The real roots are the values x such that ax^2+ bx + c = 0.
        /// </summary>
        public static double[] QuadraticEquationRealRoots(double a, double b, double c)
        {
            double d = b * b - 4.0 * a * c;
            if (d < 0)
                return new[] { Double.NaN, Double.NaN };

            double sqrtD = Math.Sqrt(d);
            double r1 = (-b + sqrtD) / (2 * a);
            double r2 = (-b - sqrtD) / (2 * a);
            return new[] { r1, r2 };
        }

        /// <summary>
        /// The Shannon entropy of a probability array.
        /// H = -Σ[p(i) x ln(p(i))].
        /// </summary>
        /// <param name="probabilityArray">The sum of this array should equal 1.</param>
        public static double Entropy(double[] probabilityArray)
        {
            double h = 0;
            for (int i = 0; i < probabilityArray.Length; i++)
            {
                if (probabilityArray[i] == 0 || probabilityArray[i] == 1) continue; // p(i) x ln(p(i)) goes to 0 in the limit, but will cause an NaN is we try to calculate. (MM, 2/24/2014)
                h -= probabilityArray[i] * Math.Log(probabilityArray[i]);
            }
            return h;
        }
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.