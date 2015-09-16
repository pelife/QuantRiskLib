using System;
using System.ComponentModel;
using System.Linq;

namespace QuantRiskLib
{

    public enum PaymentFrequency
    {
        Zero = 0, Annual = 1, SemiAnnual = 2, Quarterly = 4, Monthly = 12
    }

    ///Source: www.risk256.com
    ///
    public class BondMath
    {
        public static double CalendarDaysPerYear = 365.25;

        public enum DayCountConvention
        {
            [Description("Actual/Actual U.S.")] ActActUS, //U.S. Treasury Bonds and notes
        }

        /// <summary>
        /// Interest = Principal x Coupon x DayCountFactor
        /// </summary>
        public static double DayCountFactor(DayCountConvention dayCountConvention, DateTime date1, DateTime date2, DateTime date3, PaymentFrequency paymentFrequency)
        {
            switch(dayCountConvention)
            {
                case DayCountConvention.ActActUS:
                    return (date2 - date1).Days / ((double)paymentFrequency * (date3 - date1).Days);
                default:
                    throw new Exception(string.Format(@"No DayCountFactor method for DayCountConvention = {0}", dayCountConvention));
            }
        }

        
        /// <summary>
        /// Present value for a bond.
        /// </summary>
        public static double PresentValue(double notional, double couponRate, PaymentFrequency paymentFrequency, DateTime evaluationDate, DateTime[] paymentDates, double yield)
        {
            double pv = 0;
            double df = 0;
            double couponPayment = paymentFrequency == PaymentFrequency.Zero ? 0 : notional * couponRate / (double) paymentFrequency;
            foreach(DateTime date in paymentDates)
            {
                double years = (date - evaluationDate).Days / CalendarDaysPerYear;
                df = Math.Pow(1.0 + yield, -years);
                pv += df;
            }
            pv *= couponPayment;
            pv += df * notional;
            return pv;
        }

        public static double Yield(double notional, double couponRate, PaymentFrequency paymentFrequency, DateTime evaluationDate, DateTime[] paymentDates, double price)
        {
            double T = (paymentDates.Last() - evaluationDate).Days / CalendarDaysPerYear;
            if(paymentFrequency == 0)
                return Math.Pow(notional / price, 1.0 / T) - 1.0;

            const double tolerance = 0.00001;
            const int maxLoops = 16;

            //Initial estimate for yield is based on first order Taylor expansion and annual frequency.
            double c = couponRate * notional;
            double y = (notional + T * c - price) / (price + (T - 1) * (notional + 0.5 * T * c)); 

            double impliedPrice = PresentValue(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, y);

            int nLoops = 0;
            while (Math.Abs(impliedPrice - price) > tolerance)
            {
                nLoops++;
                if (nLoops > maxLoops)
                    throw new Exception("Yield did not converge.");
                double dollarDur = DollarDuration(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, y);
                y = y - (impliedPrice - price) / dollarDur;
                impliedPrice = PresentValue(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, y);
            }
            return y;
        }

        /// <summary>
        /// The yield calculated using the "book" value of the bond in place of the current market price.
        /// </summary>
        /// <param name="notional"></param>
        /// <param name="couponRate"></param>
        /// <param name="paymentFrequency"></param>
        /// <param name="evaluationDate"></param>
        /// <param name="paymentDates"></param>
        /// <param name="purchasePrice"></param>
        /// <param name="purchaseDate"></param>
        /// <returns></returns>
        public static double BookYield(double notional, double couponRate, PaymentFrequency paymentFrequency, DateTime evaluationDate, DateTime[] paymentDates, double purchasePrice, DateTime purchaseDate)
        {
            DateTime maturityDate = paymentDates[paymentDates.Length - 1];
            double bv = BookValue(notional, purchasePrice, purchaseDate, evaluationDate, maturityDate);
            return Yield(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, bv);
        }

        /// <summary>
        /// Book value based on straight line amortization of the initial premium or discount.
        /// </summary>
        /// <param name="notional"></param>
        /// <param name="purchasePrice"></param>
        /// <param name="evaluationDate"></param>
        /// <param name="purchaseDate"></param>
        /// <param name="maturityDate"></param>
        /// <returns></returns>
        public static double BookValue(double notional, double purchasePrice, DateTime purchaseDate, DateTime evaluationDate, DateTime maturityDate)
        {
            double initialPremium = purchasePrice - notional;
            double fractionOfTimeRemaining = (double)(maturityDate - evaluationDate).Ticks / (maturityDate - purchaseDate).Ticks;
            double premiumRemaining = fractionOfTimeRemaining * initialPremium;
            return notional + premiumRemaining;
        }

        #region Duration
        /// <summary>
        /// -(d[present value] / d[yield]) x (1 + yield)/[present value]
        /// </summary>
        public static double Duration(double couponRate, PaymentFrequency paymentFrequency, DateTime evaluationDate, DateTime[] paymentDates, double yield)
        {
            //Note: Duration does not depend on the notional of the bond. Both dVdY and pv are proportional to the notional.
            //As long as we use the same notional in both formulas, the final result will be the same. 
            double dVdY = DollarDuration(100.00, couponRate, paymentFrequency, evaluationDate, paymentDates, yield);
            double pv = PresentValue(100.00, couponRate, paymentFrequency, evaluationDate, paymentDates, yield);
            return -dVdY * (1 + yield) / pv;
        }

        /// <summary>
        /// -(d[present value] / d[yield]) x (1/[present value])
        /// </summary>
        public static double ModifiedDuration(double couponRate, PaymentFrequency paymentFrequency, DateTime evaluationDate, DateTime[] paymentDates, double yield)
        {
            //Note: Modified duration does not depend on the notional of the bond. Both dVdY and pv are proportional to the notional.
            //As long as we use the same notional in both formulas, the final result will be the same. 
            double dVdY = DollarDuration(100.00, couponRate, paymentFrequency, evaluationDate, paymentDates, yield);
            double pv = PresentValue(100.00, couponRate, paymentFrequency, evaluationDate, paymentDates, yield);
            return -dVdY / pv;
        }

        /// <summary>
        /// d[present value] / d[yield]
        /// </summary>
        public static double DollarDuration(double notional, double couponRate, PaymentFrequency paymentFrequency, DateTime evaluationDate, DateTime[] paymentDates, double yield)
        {
            double dur = 0;
            double df = 0, years = 0;
            double couponPayment = paymentFrequency == PaymentFrequency.Zero ? 0 : notional * couponRate / (double)paymentFrequency;
            foreach (DateTime date in paymentDates)
            {
                years = (date - evaluationDate).Days / CalendarDaysPerYear;
                df = Math.Pow(1.0 + yield, -years - 1);
                dur -= years * df;
            }
            dur *= couponPayment;
            dur -= years * df * notional;
            return dur;
        }

        /// <summary>
        /// ([Price minus] - [Price plus]) / (2 x [Price zero] x dY).
        /// [Price zero] is the present value of the bond evaluate at the current yield.
        /// [Price minus] and [Price plus] are the present value calculate using (yield - dY) and (yield + dY), respectively.
        /// dY = 0.0001
        /// </summary>
        public static double EffectiveDuration(double notional, double couponRate, PaymentFrequency paymentFrequency, DateTime evaluationDate, DateTime[] paymentDates, double yield)
        {
            double pZero = PresentValue(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, yield);
            double pPlus = PresentValue(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, yield + 0.0001);
            double pMinu = PresentValue(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, yield - 0.0001);
            return 10000 * (pMinu - pPlus) / (2 * pZero);
        }
        #endregion

        #region Convexity
        /// <summary>
        /// d^2[present value] / d^2[yield] x (1/[present value])
        /// </summary>
        public static double Convexity(double couponRate, PaymentFrequency paymentFrequency, DateTime evaluationDate, DateTime[] paymentDates, double yield)
        {
            //Note: Convexity does not depend on the notional of the bond. Both cnvx and pv are proportional to the notional.
            //As long as we use the same notional in both formulas, the final result will be the same. 
            const double notional = 100;

            double cnvx = 0;
            double df = 0, years = 0;
            double couponPayment =  paymentFrequency == PaymentFrequency.Zero ? 0 : notional * couponRate / (double)paymentFrequency;
            foreach (DateTime date in paymentDates)
            {
                years = (date - evaluationDate).Days / CalendarDaysPerYear;
                df = Math.Pow(1.0 + yield, -years - 2);
                cnvx += years * (1 + years) * df;
            }
            cnvx *= couponPayment;
            cnvx += years * df * notional;
            double pv = PresentValue(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, yield);
            return cnvx / pv;
        }

        /// <summary>
        /// ([Price minus] + [Price plus]  - 2 x [Price zero]) / ([Price zero] x dY^2).
        /// [Price zero] is the present value of the bond evaluate at the current yield.
        /// [Price minus] and [Price plus] are the present value calculate using (yield - dY) and (yield + dY), respectively.
        /// dY = 0.0001
        /// </summary>
        public static double EffectiveConvexity(double notional, double couponRate, PaymentFrequency paymentFrequency, DateTime evaluationDate, DateTime[] paymentDates, double yield)
        {
            double pZero = PresentValue(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, yield);
            double pPlus = PresentValue(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, yield + 0.0001);
            double pMinu = PresentValue(notional, couponRate, paymentFrequency, evaluationDate, paymentDates, yield - 0.0001);
            return 100000000 * (pMinu + pPlus - 2 * pZero) / (pZero);
        }
        #endregion

        /// <summary>
        /// Converts a series of yields to one period forward rates. Assumes all yields are equally spaced (e.g. 1yr, 2yr, 3yr, ...).
        /// (1 + f0)(1 + f1)...(1 + fn) = (1 + yn) for all input yields, y0, y1, ..., yN.
        /// </summary>
        public static double[] YieldsToForwardRates(double[] yields)
        {
            double[] times = new double[yields.Length];
            for(int i = 0; i < yields.Length; i++)
                times[i] = i + 1;

            return YieldsToForwardRates(yields, times);
        }

        /// <summary>
        /// Converts a series of yields to one period forward rates.
        /// (1 + f0)(1 + f1)...(1 + fn) = (1 + yn) for all input yields, y0, y1, ..., yN.
        /// </summary>
        /// <param name="yields"></param>
        /// <param name="times">Time to maturity for each yield in yields. Typically in years (e.g. {1, 2, 3} for 1, 2 and 3 year yields).</param>
        /// <returns></returns>
        public static double[] YieldsToForwardRates(double[] yields, double[] times)
        {
            if(yields.Length != times.Length)
                throw new ArgumentException("Length of yields and times must be equal.");

            int n = yields.Length;
            double[] fwdRates = new double[n];

            for (int i = 0; i < n; i++)
            {
                fwdRates[i] = Math.Pow(1.0 + yields[i], times[i]);
                double d = 1;
                for (int j = i - 1; j >= 0; j--)
                    d *= fwdRates[j];
                fwdRates[i] /= d;
            }

            for (int i = 0; i < n; i++)
                fwdRates[i] -= 1.0;

            return fwdRates;
        }
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.