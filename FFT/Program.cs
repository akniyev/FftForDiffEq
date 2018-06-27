using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using static System.Math;

namespace FFT
{
    internal static class Program
    {
        private static Complex[] expArray;
        private static Complex[] expArrayHalf;
        
        private static Complex[] IFft(Complex[] array)
        {
            var n = array.Length;
            
            if (n % 2 != 0 && n != 1) throw new ArgumentException();
            
            if (n == 1)
            {
                return array;
            }
            
            var omegaN =  new Complex(Cos(2 * PI / n), Sin(2 * PI / n));
            
            var omega = new Complex(1, 0);
            
            var even = new Complex[n / 2];
            var odd = new Complex[n / 2];
            for (var i = 0; i < n / 2; i++)
            {
                even[i] = array[2 * i];
                odd[i] = array[2 * i + 1];
            }

            var y0 = IFft(even);
            var y1 = IFft(odd);

            var y = new Complex[n];
            
            for (var k = 0; k < n/2; k++)
            {
                y[k] = y0[k] + omega * y1[k];
                y[k + n / 2] = y0[k] - omega * y1[k];
                omega = omega * omegaN;
            }
            
            return y;
        }
        private static Complex[] RealInvFft(double[] a)
        {
            var n = a.Length;

            var imag = new Complex(0, 1);

            var b = new Complex[n / 2];
            var h = new Complex[n / 2];
            for (int i = 0; i < n / 2; i++)
            {
                h[i] = new Complex(a[2 * i], a[2 * i + 1]);
            }

            var a1 = IFft(h);

            b[0] = new Complex(a1[0].Real, -a1[0].Imaginary);
            for (int i = 1; i < n / 2; i++)
            {
                b[i] = new Complex(a1[n / 2 - i].Real, -a1[n / 2 - i].Imaginary);
            }
            var g1 = new Complex[n / 2];
            var g2 = new Complex[n / 2];
            var b1 = new Complex[n];
            for (int i = 0; i < n / 2; i++)
            {
                g1[i] = (a1[i] + b[i]) / 2;
                g2[i] = -imag * (a1[i] - b[i]) / 2;
            }
            for (int i = 0; i < n / 2; i++)
            {
                b1[i] = g1[i] + expArrayHalf[i] * g2[i];
                b1[i + n / 2] = g1[i] - expArrayHalf[i] * g2[i];
            }

            return b1;
        }
        private static Complex[] RealInvFft2(double[] a)
        {
            var n = a.Length;
            var f = new double[2 * n];
            f[0] = 0;
            for (var i = 1; i < n; i++)
            {
                f[i] = a[i] / i;
            }
            f[n] = 0;
            for (var i = 1; i < n; i++)
            {
                f[2 * n - i] = -a[i] / i;
            }
            var expArray = new Complex[n];

            for (var k = 0; k < n; k++)
            {
                expArray[k] = new Complex(Cos(PI * k / n), Sin(PI * k / n));
            }
            var imag = new Complex(0, 1);

            var b = new Complex[n];
            var h = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                h[i] = new Complex(f[2 * i], f[2 * i + 1]);
            }

            var a1 = IFft(h);

            b[0] = new Complex(a1[0].Real, -a1[0].Imaginary);
            for (int i = 1; i < n; i++)
            {
                b[i] = new Complex(a1[n - i].Real, -a1[n - i].Imaginary);
            }
            var g1 = new Complex[n];
            var g2 = new Complex[n];
            var b1 = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                g1[i] = (a1[i] + b[i]) / 2;
                g2[i] = -imag * (a1[i] - b[i]) / 2;
            }
            for (int i = 0; i < n; i++)
            {
                b1[i] = g1[i] + expArray[i] * g2[i];
                //b1[i + n] = g1[i] - expArray[i] * g2[i];
            }
            return b1;
        }
        private static Complex[] NewIFft(Complex[] a)
        {
            var n = a.Length;

//            var expArray = new Complex[n];

//            for (var k = 0; k < n; k++)
//            {
//                expArray[k] = new Complex(Cos(PI * k / n), Sin(PI * k / n));
//            }
            var imag = new Complex(0, 1);

            var b = new Complex[n];

            var a1 = IFft(a);

            b[0] = new Complex(a1[0].Real, -a1[0].Imaginary);
            for (int i = 1; i < n; i++)
            {
                b[i] = new Complex(a1[n - i].Real, -a1[n - i].Imaginary);
            }
            var g1 = new Complex[n];
            var g2 = new Complex[n];
            var b1 = new Complex[2 * n];
            for (int i = 0; i < n; i++)
            {
                g1[i] = (a1[i] + b[i]) / 2;
                g2[i] = -imag * (a1[i] - b[i]) / 2;
            }
            for (int i = 0; i < n; i++)
            {
                b1[i] = g1[i] + expArray[i] * g2[i];
                //b1[i + n] = g1[i] - expArray[i] * g2[i];
            }
            return b1;
        }
        private static bool IsPowerOfTwo(int x)
        {
            return (x & (x - 1)) == 0;
        }

        private static IEnumerable<Complex> CalculateOld(IReadOnlyList<double> a)
        {
            var n = a.Count;
            var b1 = new double[2 * n];
            b1[0] = 0;
            for (var i = 1; i < n; i++)
            {
                b1[i] = a[i] / i;
            }
            b1[n] = 0;
            for (var i = 1; i < n; i++)
            {
                b1[2 * n - i] = -a[i] / i;
            }
            var a2 = IFft(b1.Select(x => new Complex(x, 0)).ToArray());
            var c1 = new Complex[n];
            for (var i = 0; i < c1.Length; i++)
            {
                c1[i] = Sqrt(2.0) / PI * a2[i] / ImagNum;
            }

            return c1;
        }

        private static IEnumerable<Complex> CalculateNew(IReadOnlyList<double> a)
        {
            var n = a.Count;
            var b = new double[2 * n];         
            b[0] = 0;
            for (var i = 1; i < n; i++)
            {
                b[i] = a[i] / i;
            }
            b[n] = 0;
            for (var i = 1; i < n; i++)
            {
                b[2 * n - i] = -a[i] / i;
            }
            var h = new Complex[n];
            for (int i = 0; i < n; i++)
            {
                h[i] = new Complex(b[2 * i], b[2 * i + 1]);
            }
           
            var a1 = NewIFft(h);
            var c = new Complex[n];
            
            for (var i = 0; i < c.Length; i++)
            {
                c[i] = Sqrt(2.0) / PI * a1[i] / ImagNum;
            }

            return c;
        }
        
        private static readonly Complex ImagNum = new Complex(0, 2);

        private static double[] GenerateRandomArray(int n)
        {
            var r = new Random();
            var result = new double[n];

            for (var i = 0; i < n; i++)
            {
                result[i] = r.NextDouble() * 100 - 50;
            }

            return result;
        }

        private static void Main(string[] args)
        {
            for (var i = 5; i < 16; i++)
            {
                var n = (int)Pow(2, i);
                    
                int iterationCount          = 1000;

                Console.Write("N = {0}; \t\t\t\t", n);
            
                expArray = new Complex[n];

                for (var k = 0; k < n; k++)
                {
                    expArray[k] = new Complex(Cos(PI * k / n), Sin(PI * k / n));
                }
            
                expArrayHalf = new Complex[n / 2];

                for (var k = 0; k < n / 2; k++)
                {
                    expArrayHalf[k] = new Complex(Cos(2 * PI * k / n), Sin(2 * PI * k / n));
                }

                var r = new Random();

                //Для сетки 1

                var sw = Stopwatch.StartNew();

                for (var iteration = 0; iteration < iterationCount; iteration++)
                {
                    sw.Stop();
                    var a = GenerateRandomArray(n);
                    sw.Start();
                    var c = CalculateNew(a);
                }
            
                var g1f = sw.Elapsed.TotalSeconds;
                Console.WriteLine("G1F {0:0.0000000}", g1f / iterationCount);
            }
//            int iterationCount          = 1000;
//            int n                       = 8192 * 8;
//
//            Console.WriteLine("N = {0}", n);
//            
//            expArray = new Complex[n];
//
//            for (var k = 0; k < n; k++)
//            {
//                expArray[k] = new Complex(Cos(PI * k / n), Sin(PI * k / n));
//            }
//            
//            expArrayHalf = new Complex[n / 2];
//
//            for (var k = 0; k < n / 2; k++)
//            {
//                expArrayHalf[k] = new Complex(Cos(2 * PI * k / n), Sin(2 * PI * k / n));
//            }
//
//            var r = new Random();
//
//            //Для сетки 1
//
//            var sw = Stopwatch.StartNew();
//
//            for (var iteration = 0; iteration < iterationCount; iteration++)
//            {
//                sw.Stop();
//                var a = GenerateRandomArray(n);
//                sw.Start();
//                var c = CalculateNew(a);
//            }
//            
//            var g1f = sw.Elapsed.TotalSeconds;
//            Console.WriteLine("G1F", g1f);


//            Console.WriteLine("Быстрое синус преобразование для первой сетки по новому");
//            Console.WriteLine("TIME:" + g1f / iterationCount);


//            sw.Restart();
//
//            for (var iteration = 0; iteration < iterationCount; iteration++)
//            {
//                sw.Stop();
//                var a = GenerateRandomArray(n);
//                sw.Start();
//                var c1 = CalculateOld(a);
//            }
//            
//            var g1ff = sw.Elapsed.TotalSeconds;


            
            
            
            
//            Console.WriteLine("Быстрое синус преобразование для первой сетки по старому");
//            Console.WriteLine("TIME:" + g1ff / iterationCount);

//            Console.WriteLine("Соотношение: {0}", g1ff / g1f);
//            foreach (var t in c1)
//            {
//                Console.WriteLine(t.Real);
//            }

            //sw.Restart();
            //var d = new double[n];
            //Console.WriteLine("Синус преобразование по формуле для первой сетки");
            //for (var j = 0; j < n; j++)
            //{
            //    double s = 0;
            //    for (var k = 1; k < n; k++)
            //    {
            //        s = s + 1.0 * a[k] * Sin(PI * k * j / n) / k;
            //    }
            //    d[j] = Sqrt(2.0) / PI * s;
            //}
            //var g1s = sw.Elapsed.TotalSeconds;
            //Console.WriteLine("TIME:" + g1s);

            //foreach (var t in d)
            //{
            //    Console.WriteLine(t);
            //}

            //Для сетки 2
//            sw.Restart();
//            var p = new Complex[2 * n];
//            var expArray = new Complex[2 * n];
//            sw.Restart();
//            for (var k = 0; k < 2 * n; k++)
//            {
//                expArray[k] = new Complex(Cos(PI * k / (2 * n)), Sin(PI * k / (2 * n)));
//            }
//            p[0] = 0;
//            for (var i = 1; i < n; i++)
//            {
//                p[i] = a[i] * expArray[i] / i;
//            }
//            p[n] = 0;
//            for (var i = n + 1; i < 2 * n; i++)
//            {
//                //p[2 * n - i] = a[i] * expArray[2 * n - i] / i;
//                p[i] = a[2 * n - i] * expArray[i] / (2 * n - i);
//            }
//
//            var a3 = IFft(p);
//            var q = new Complex[n];
//            for (var i = 0; i < n; i++)
//            {
//                q[i] = Sqrt(2.0) / PI * a3[i] / ImagNum;
//            }
//            Console.WriteLine("Быстрое синус преобразование для второй сетки");
//
//            var g2f = sw.Elapsed.TotalSeconds;
//            Console.WriteLine("TIME:" + g2f);
//
//            //for (var i = 0; i < q.Length; i++)
//            //{
//            //    Console.WriteLine(q[i].Real);
//            //}
//
//            sw.Restart();
//            var d1 = new double[n];
//
//            Console.WriteLine("Синус преобразование по формуле для второй сетки");
//
//            for (var j = 0; j < n; j++)
//            {
//                double s1 = 0;
//                for (var k = 1; k < n; k++)
//                {
//                    s1 = s1 + a[k] * Sin(PI * k * (2 * j + 1) / (2 * n)) / k;
//                }
//                d1[j] = Sqrt(2.0) / PI * s1;
//            }
//            var g2s = sw.Elapsed.TotalSeconds;
//            Console.WriteLine("TIME:" + g2s);
//
//            //foreach (var t in d1)
//            //{
//            //    Console.WriteLine(t);
//            //}
            Console.ReadLine();
        }
    }
}