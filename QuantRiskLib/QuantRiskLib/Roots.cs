using System;

namespace QuantRiskLib
{
    ///Source: www.risk256.com
    ///

    public class Roots
    {
        /// <summary>
        /// Uses the secant method to find a root of a function. For a function y(x), the root, r, is the place where y(r) = 0. 
        /// Note if the function has more than one root, this function will only find one (and may not find any).
        /// </summary>
        /// <param name="function">Function for which we want to find the roots.</param>
        /// <param name="initialEstimate">The initial estimate of the root. This is the starting point for the search.</param>
        /// <param name="convergenceCriteria">The algorith will stop when the |y(r)| is less than convergenceCriteria or it has run through maxLoops.</param>
        /// <param name="maxLoops">The algorith will stop when the |y(r)| is less than convergenceCriteria or it has run through maxLoops.</param>
        /// <param name="delta">Used to calculate the slope of the function. Slope(x) ≈ [f(x+delta) - f(x)]/delta. Making too small can cause precision issues.</param>
        public static double SecantMethod(Func<double, double> function, double initialEstimate = 0, double convergenceCriteria = 0.000001, int maxLoops = 128, double delta = 0.0001)
        {
            double root = initialEstimate;
            for (int i = 0; i < maxLoops; i++)
            {
                double y = function(root);
                //Console.WriteLine("{0, 5} => {1}", i, y);
                if (Math.Abs(y) < convergenceCriteria) break;
                double slope = (function(root + delta) - y) / delta; //(function(root + delta) - function(root - delta)) / (2 * delta) could be more accurate, but is definitley slower.
                root -= y / slope;
            }
            return root;
        }

        /// <summary>
        /// Used to find a root of a function. For a function y(x), the root, r, is the place where y(r) = 0. 
        /// Not sure what the official name of this method is. It can be looked at as a finite difference version of Halley's method or an extended secant method.
        /// Note if the function has more than one root, this function will only find one (and may not find any).
        /// </summary>
        /// <param name="function">Function for which we want to find the roots.</param>
        /// <param name="initialEstimate">The initial estimate of the root. This is the starting point for the search.</param>
        /// <param name="convergenceCriteria">The algorith will stop when the |y(r)| is less than convergenceCriteria or it has run through maxLoops.</param>
        /// <param name="maxLoops">The algorith will stop when the |y(r)| is less than convergenceCriteria or it has run through maxLoops.</param>
        /// <param name="delta">Used to calculate the slope and curvature of the function. Slope(x) ≈ [f(x+delta) - f(x)]/delta. Making too small can cause precision issues.</param>
        public static double SecantMethodExtended(Func<double, double> function, double initialEstimate = 0, double convergenceCriteria = 0.000001, int maxLoops = 128, double delta = 0.0001)
        {
            double root = initialEstimate;
            for (int i = 0; i < maxLoops; i++)
            {
                double y = function(root);
                //Console.WriteLine("{0, 5} => {1}", i, y);
                if (Math.Abs(y) < convergenceCriteria) break;
                double yMnus = function(root - delta);
                double yPlus = function(root + delta);
                double slope = (yPlus - yMnus) / (2 * delta);
                double curve = ((yPlus - y) - (y - yMnus)) / (delta * delta);
                root -= 2 * y * slope / (2 * slope * slope - y * curve);
            }
            return root;
        }

        /// <summary>
        /// Uses Newton–Raphson method to find a root of a function. For a function y(x), the root, r, is the place where y(r) = 0. 
        /// Note if the function has more than one root, this function will only find one (and may not find any).
        /// </summary>
        /// <param name="function">Function for which we want to find the roots.</param>
        /// <param name="firstDerivative">First derivative of function.</param>
        /// <param name="initialEstimate">The initial estimate of the root. This is the starting point for the search.</param>
        /// <param name="convergenceCriteria">The algorith will stop when the |y(r)| is less than convergenceCriteria or it has run through maxLoops.</param>
        /// <param name="maxLoops">The algorith will stop when the |y(r)| is less than convergenceCriteria or it has run through maxLoops.</param>
        public static double NewtonRaphson(Func<double, double> function, Func<double, double> firstDerivative, double initialEstimate = 0, double convergenceCriteria = 0.000001, int maxLoops = 128)
        {
            double root = initialEstimate;
            for (int i = 0; i < maxLoops; i++)
            {
                double y = function(root);
                //Console.WriteLine("{0, 5} => {1}", i, y);
                if (Math.Abs(y) < convergenceCriteria) break;
                double slope = firstDerivative(root);
                root -= y / slope;
            }
            return root;
        }

        /// <summary>
        /// Uses Halley's method to find a root of a function. For a function y(x), the root, r, is the place where y(r) = 0. 
        /// Note if the function has more than one root, this function will only find one (and may not find any).
        /// </summary>
        /// <param name="function">Function for which we want to find the roots.</param>
        /// <param name="firstDerivative">First derivative of function.</param>
        /// <param name="secondDerivative">Second derivative of function.</param>
        /// <param name="initialEstimate">The initial estimate of the root. This is the starting point for the search.</param>
        /// <param name="convergenceCriteria">The algorith will stop when the |y(r)| is less than convergenceCriteria or it has run through maxLoops.</param>
        /// <param name="maxLoops">The algorith will stop when the |y(r)| is less than convergenceCriteria or it has run through maxLoops.</param>
        public static double HalleysMethod(Func<double, double> function, Func<double, double> firstDerivative, Func<double, double> secondDerivative, double initialEstimate = 0, double convergenceCriteria = 0.000001, int maxLoops = 128)
        {
            double root = initialEstimate;
            for (int i = 0; i < maxLoops; i++)
            {
                double y = function(root);
                //Console.WriteLine("{0, 5} => {1}", i, y);
                if (Math.Abs(y) < convergenceCriteria) break;
                double slope = firstDerivative(root);
                double curve = secondDerivative(root);
                root -= 2 * y * slope / (2 * slope * slope - y * curve);
            }
            return root;
        }
    
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.
