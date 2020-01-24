using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace JOSTester
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine(Globals.qe.ToString("g4", System.Globalization.CultureInfo.InvariantCulture));
            JOS.Solver testowy = new JOS.Solver();
            testowy.dumpVersion();
            Console.WriteLine(testowy.f(0.1, 0.1, 0.1, 450e-7, 1.85, 9, 1e-20, 1e-20, 1e-20).ToString("g4", System.Globalization.CultureInfo.InvariantCulture));
            JOS.Multiplet exp = new JOS.Multiplet();
            //    
            exp.o2 = 1.0e-20;
            exp.o4 = 1.1e-20;
            exp.o6 = 1.2e-20;
            exp.LoadFromFile(@"r:\testszarp.txt");
            //Console.WriteLine( testowy.chi2(exp));
            //MathNet.Numerics.LinearAlgebra.Matrix<double> Hess ;
            //MathNet.Numerics.LinearAlgebra.Matrix<double> Grad;
            //testowy.CalculateHessian(exp, out Hess,out Grad);
            //Console.Write(Hess.ToString());
            string msg, lat;
            msg = "";
            testowy.FitLM(exp,out msg,out lat);
            //Console.WriteLine(msg);
            exp.CalculateRates();
            exp.ReportRates(out msg);
            //Console.WriteLine(exp.o2);
            Console.WriteLine("Report");
            Console.WriteLine(msg);
        }
    }
}
