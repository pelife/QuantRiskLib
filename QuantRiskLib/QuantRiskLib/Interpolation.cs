using System;

namespace QuantRiskLib
{
    public class Interpolation
    {
        /// <summary>
        /// Given two points (x1, y1) and (x2, y2), determine the y value at xStar.
        /// If xStar is not between x1 and x2, this method will still return a value (technically extrapolation, not interpolation).
        /// x1 and x2 must have different values.
        /// </summary>
        public static double LinearInterpolation(double x1, double x2, double y1, double y2, double xStar)
        {
            if(x1 == x2)
                throw new ArgumentException("x1 and x2 cannot be eqaul.");
            double wt1 = (x2 - xStar) / (x2 - x1);
            return wt1 * y1 + (1 - wt1) * y2;
        }

        /// <summary>
        /// Given two points (x1, y1) and (x2, y2), determine the y value at xStar.
        /// If xStar is not between x1 and x2, this method will still return a value (technically extrapolation, not interpolation).
        /// x1 and x2 must have different values.
        /// </summary>
        public static double LinearInterpolation(double x1, double x2, double y1, double y2, double xStar, out double wt1, out double wt2)
        {
            if (x1 == x2)
                throw new ArgumentException("x1 and x2 cannot be eqaul.");
            wt1 = (x2 - xStar) / (x2 - x1);
            wt2 = 1 - wt1;
            return wt1 * y1 + (1 - wt1) * y2;
        }

        /// <summary>
        /// Given four points -- p1, p2, p3, p4 -- determine the z value for a fifth point with (x, y) coordinates (xStar, yStar)
        /// Points should be arrayed as:
        /// 2 - 4
        /// |   |
        /// 1 - 3
        /// </summary>
        public static double BilinearInterpolation(ThreeDPoint p1, ThreeDPoint p2, ThreeDPoint p3, ThreeDPoint p4, double xStar, double yStar, out double wt1, out double wt2, out double wt3, out double wt4)
        {
            if(!(p1.X == p2.X && p1.Y == p3.Y && p2.Y == p4.Y && p3.X == p4.X))
                throw new ArgumentException("Points are not lined up correctly");

            double x1 = p1.X;
            double x2 = p3.X;
            double y1 = p1.Y;
            double y2 = p2.Y;

            double k = 1.0 / ((x2 - x1) * (y2 - y1));
            wt1 = k * (x2 - xStar) * (y2 - yStar);
            wt2 = k * (x2 - xStar) * (yStar - y1);
            wt3 = k * (xStar - x1) * (y2 - yStar);
            wt4 = k * (xStar - x1) * (yStar - y1);
            double zStar = wt1 * p1.Z + wt2 * p2.Z + wt3 * p3.Z + wt4 * p4.Z;
            return zStar;
        }
    }

    public class ThreeDPoint
    {
        public double X;
        public double Y;
        public double Z;

        public ThreeDPoint(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public bool Equals(ThreeDPoint point)
        {
            return point.X == X && point.Y == Y && point.Z == Z;
        }
    }
}
