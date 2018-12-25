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
            Console.WriteLine(JOS.Solver.qe.ToString("g4", System.Globalization.CultureInfo.InvariantCulture));
            JOS.Solver testowy = new JOS.Solver();
            testowy.dumpVersion();
            Console.WriteLine(testowy.f(0.1, 0.1, 0.1, 450e-7, 1.85, 9, 1e-20, 1e-20, 1e-20).ToString("g4", System.Globalization.CultureInfo.InvariantCulture));
            JOS.Experiment exp = new JOS.Experiment();
            /*
            exp.o2 = 1e-20;
            exp.o4 = 1e-20;
            exp.o6 = 1e-20;
            exp.n = 1.85;
            exp.TwoJPlusOne = 9;
            exp.fexp.Add(1e-6);
            exp.fexp.Add(2e-6);
            exp.u2.Add(0.02);
            exp.u2.Add(0.05);
            exp.u4.Add(0.0223);
            exp.u4.Add(0.055);
            exp.u6.Add(0.032);
            exp.u6.Add(0.015);
            exp.lambda0.Add(450e-7);
            exp.lambda0.Add(512e-7);
            Console.WriteLine(testowy.chi2(exp));
            */
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
            Console.WriteLine(msg);
        }
    }
}
