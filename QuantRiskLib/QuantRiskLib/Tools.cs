using System;

namespace QuantRiskLib
{
    ///Source: www.risk256.com

    public class Tools
    {

        /// <summary>
        /// Returns the most recent points from an array (the points with the highest index values).
        /// For functions in QuantRiskLib with a decay factor, the weight on the last element in the array is the highest.
        /// In many applications this corresponds with the assumption that the last element in the array represents the most recent data point.
        /// </summary>
        /// <param name="inArray">Array from which the data points will be selected.</param>
        /// <param name="length">The number of data points to be returned.</param>
        /// <returns></returns>
        public static double[] MostRecentValues(double[] inArray, int length)
        {
            if(inArray.Length == length) return inArray;

            if (inArray.Length < length)
                throw new ArgumentException("length of inArray must be greater than or equal to length");

            double[] outArray = new double[length];
            for (int i = 0; i < length; i++)
                outArray[i] = inArray[inArray.Length - length + i];
            return outArray;
        }

        /// <summary>
        /// Returns the most recent points from an array (the points with the highest index values).
        /// For functions in QuantRiskLib with a decay factor, the weight on the last element in the array is the highest.
        /// In many applications this corresponds with the assumption that the last element in the array represents the most recent data point.
        /// </summary>
        /// <param name="inArray">Array from which the data points will be selected. Array is TxN, where T is number of time periods, and N is number of variables.</param>
        /// <param name="length">The number of data points to be returned.</param>
        /// <returns></returns>
        public static double[,] MostRecentValues(double[,] inArray, int length)
        {
            int T = inArray.GetLength(0);
            int N = inArray.GetLength(1);
            if (T < length)
                throw new ArgumentException("length of inArray must be greater than or equal to length");

            if(T == length) return inArray;

            double[,] outArray = new double[length, N];
            
            // TODO: Test if this works.
            // Array.Copy(inArray, outArray, T*N);

            for (int i = 0; i < length; i++)
                for (int j = 0; j < N; j++)
                    outArray[i, j] = inArray[T - length + i, j];
            return outArray;
        }

        /// <summary>
        /// Returns a two dimensional array, length x (nLags + 1), where the first column is the most recent values of inArray, 
        /// the second column is the first lag of inArray (i.e. the last value in the column is the 2nd most recent point),
        /// the third column is the second lage of inArray, and so on.
        /// </summary>
        /// <param name="inArray">Time series of one variable.</param>
        /// <param name="length"></param>
        /// <param name="nLags"></param>
        /// <returns></returns>
        public static double[,] MostRecentValuesPlusLagArray(double[] inArray, int length, int nLags)
        {
            if (inArray.GetLength(0) < length + nLags)
                throw new ArgumentException("lenght of inArray must be greater than or equal to length + nLags");

            double[,] outArray = new double[length, nLags + 1];
            for (int i = 0; i < length; i++)
                for (int j = 0; j < (nLags + 1); j++)
                    outArray[i, j] = inArray[inArray.Length - length + i - j];

            return outArray;
        }

        /// <summary>
        /// Returns the lagged values of the most recent values.
        /// If inArray is (t-3), (t-2), (t-1), (t), with (t) being the most recent point, and nLags = 1, then the output array will be (t-3), (t-2), (t-1).
        /// </summary>
        /// <param name="inArray">Time series of one variable.</param>
        /// <param name="length"></param>
        /// <param name="nLags"></param>
        /// <returns></returns>
        public static double[] LagArray(double[] inArray, int length, int nLags)
        {
            if (inArray.GetLength(0) < length + nLags)
                throw new ArgumentException("lenght of inArray must be greater than or equal to length + nLags");
            double[] outArray = new double[length];
            for (int i = 0; i < length; i++)
                outArray[i] = inArray[inArray.Length - length + i - nLags];
            return outArray;
        }

        /// <summary>
        /// Returns true if all of the elements in the array are equal; otherwise, false.
        /// </summary>
        public static bool ArrayAllEqual(double[] array)
        {
            for (int i = 1; i < array.Length; i++)
            {
                if (array[i] == array[0]) continue;
                return false;
            }
            return true;
        }

        /// <summary>
        /// Returns true if all of the corresponding elements in the arrays are equal.
        /// That is a[0] = b[0], a[1] = b[1], ... a[n-1] = b[n-1].
        /// Returns false if arrays are not of equal length.
        /// </summary>
        public static bool ArraysEqual(double[] a, double[] b)
        {
            if (a.Length != b.Length) return false;

            for (int i = 0; i < a.Length; i++)
            {
                if (a[i] == b[i]) continue;
                return false;
            }
            return true;
        }
    }
}

//Disclaimer
//This code is freeware. The methods are not proprietary. Feel free to use, modify and redistribute. That said, if you plan 
//to use or redistribute give credit where credit is due and provide a link back to Risk256.com (or don't remove the link 
//and references already in the code). The code is intended primarily as an educational tool. No warranty is made as to the 
//code's accuracy. Use at your own risk.