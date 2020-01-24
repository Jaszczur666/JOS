using MathNet.Numerics.LinearAlgebra;
using System;
using static Globals;
using static System.Math;
public static class Globals
{
    public const double pi = 3.141592653589793238462643383279502884;
    public const double m = 9.1093897e-28;//SI 9.10938291e-31;
    public const double c = 2.99792458e10;//SI 299792458;
    public const double h = 6.6260755e-27; //SI6.62606957e-34;
    public const double qe = 4.8032068e-10;
}
namespace JOS
{


    public class Solver
    {
        public void DumpVersion()
        {
            Console.WriteLine("0.0.nic");
        }

        //public Experiment exp;
        public Solver()
        {

        }
        public double F(double u2, double u4, double u6, double lambda0, double n, double TwoJPlusOne, double o2, double o4, double o6)
        {
            //n=1.469+3927.0/pow((1e7/lambda0),2);//test
            //Console.WriteLine(u2+" " + u4 + " " + u6 + " " + o2 + " " + o4 + " " + o6);
            return (((Pow(n * n + 2, 2.0) / (9 * n)) * (8 * pi * pi * m * c)) / (3 * h * TwoJPlusOne * lambda0)) * (u2 * o2 + u4 * o4 + u6 * o6);

        }
        public double Chi2(Multiplet exp)
        {
            int size;
            double tempchi2 = 0;
            //Console.WriteLine(exp.o2 + " " + exp.o4 + " " + exp.o6);
            size = exp.u2.Count;
            for (int i = 0; i < size; i++)
            {
                double tf = F(exp.u2[i], exp.u4[i], exp.u6[i], exp.lambda0[i], exp.n, exp.TwoJPlusOne, exp.o2, exp.o4, exp.o6);
                //Console.WriteLine(tf.ToString("G4")+" "+exp.fexp[i].ToString("G4"));
                //Console.WriteLine(tf/1.488E-07);
                tempchi2 += Pow(tf - exp.fexp[i], 2);
            };
            return tempchi2;
        }
        public double Chi2new(Multiplet exp, double no2, double no4, double no6)
        {
            double to2 = exp.o2;
            double to4 = exp.o4;
            double to6 = exp.o6;
            double rv;
            exp.o2 = no2;
            exp.o4 = no4;
            exp.o6 = no6;
            rv = Chi2(exp);
            exp.o2 = to2;
            exp.o4 = to4;
            exp.o6 = to6;
            return rv;

        }
        public double Residue(double fexp, double o2, double o4, double o6, double u2, double u4, double u6, double lambda, double n, double TwoJPlusOne)
        {
            double res = F(u2, u4, u6, lambda, n, TwoJPlusOne, o2, o4, o6) - fexp;
            return res;
        }
        public void CalculateHessian(Multiplet exp, out Matrix<double> Hess, out Matrix<double> Grad)
        {
            int i, size;
            //MathNet.Numerics.LinearAlgebra.Matrix<double> Hessian;
            double ro2, ro4, ro6, res, delta;
            size = exp.u2.Count;
            Matrix<double> Res = Matrix<double>.Build.Dense(size, 1);
            delta = 1e-26;
            Matrix<double> Jaco = Matrix<double>.Build.Dense(size, 3);
            //cout <<"Debug " <<o2<<" "<<o4<<" "<<o6<<std::endl;
            for (i = 0; i < size; i++)
            {
                res = Residue(exp.fexp[i], exp.o2, exp.o4, exp.o6, exp.u2[i], exp.u4[i], exp.u6[i], exp.lambda0[i], exp.n, exp.TwoJPlusOne);
                ro2 = Residue(exp.fexp[i], exp.o2 + delta, exp.o4, exp.o6, exp.u2[i], exp.u4[i], exp.u6[i], exp.lambda0[i], exp.n, exp.TwoJPlusOne);
                ro4 = Residue(exp.fexp[i], exp.o2, exp.o4 + delta, exp.o6, exp.u2[i], exp.u4[i], exp.u6[i], exp.lambda0[i], exp.n, exp.TwoJPlusOne);
                ro6 = Residue(exp.fexp[i], exp.o2, exp.o4, exp.o6 + delta, exp.u2[i], exp.u4[i], exp.u6[i], exp.lambda0[i], exp.n, exp.TwoJPlusOne);
                Jaco[i, 0] = (ro2 - res) / delta;
                Jaco[i, 1] = (ro4 - res) / delta;
                Jaco[i, 2] = (ro6 - res) / delta;
                Res[i, 0] = res;

            };
            //Console.WriteLine("Jacobian");
            //Console.WriteLine(Jaco.ToString());
            Hess = Jaco.Transpose() * Jaco;
            Grad = Jaco.Transpose() * Res;
            /*		std::cout <<"_________________________________"<<std::endl;
                    std::cout << Hess<<std::endl;
                    std::cout <<"_________________________________"<<std::endl;*/
        }
        double Exponent(double x)
        {
            return Floor(Log10(x));
        }
        double Mantissa(double x)
        {
            return x / Pow(10.0, Exponent(x));
        }
        string Lf(double x) // Latex format
        {
            return Mantissa(x).ToString("G3") + "$ \\cdot 10 ^{" + Exponent(x) + "}$";
        }
        public void FitLM(Multiplet exp, out string MSG, out string LATEX)
        {
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            //Matrix<double> Grad;
            Matrix<double> Hessiandiag;
            Matrix<double> parameters = Matrix<double>.Build.Dense(3, 1);
            Matrix<double> newparams;// = Matrix<double>.Build.Dense(3, 1);
            //Vector<double> error;// = Vector<double>.Build.Dense(3);
            LATEX = "";
            parameters[0, 0] = exp.o2;
            parameters[1, 0] = exp.o4;
            parameters[2, 0] = exp.o6;
            double no2, no4, no6, sumfexp, sumdfexp;
            double lambda, chi2s, chi2n;
            lambda = 1 / 1024.0;
            sumdfexp = 0;
            sumfexp = 0;
//            chi2s = 0;
            chi2n = 0;
            Hessiandiag = Matrix<double>.Build.Dense(3, 3, 1);
            MSG = "Num.\tChi2\tO2\tO4\tO6\r\n";
            Console.WriteLine("Beginnig fitting procedure.");
            for (int i = 1; i < 10; i++)
            {
                chi2s = Chi2(exp);
                CalculateHessian(exp, out Matrix<double> Hessian, out Matrix<double> Grad);
                /*Console.WriteLine("gradient");
                Console.WriteLine(Grad.ToString());
                Console.WriteLine("hessian");
                Console.WriteLine(Hessian.ToString());
                */
                Hessiandiag = Matrix<double>.Build.DiagonalOfDiagonalVector(Hessian.Diagonal());
                //Console.WriteLine(Hessiandiag.ToString());
                newparams = parameters - ((Hessian + lambda * Hessiandiag).Inverse() * Grad);
                //Console.WriteLine(newparams.ToString());
                no2 = (newparams[0, 0]);
                no4 = (newparams[1, 0]);
                no6 = (newparams[2, 0]);

                chi2n = Chi2new(exp, no2, no4, no6);
                Console.WriteLine(i + " " + chi2s.ToString("G6") + " " + no2.ToString("G6") + " " + no4.ToString("G6") + " " + no6.ToString("G6"));
                MSG += i.ToString() + "\t" + chi2s.ToString("G6") + "\t" + no2.ToString("G6") + "\t" + no4.ToString("G6") + "\t" + no6.ToString("G6") + "\r\n";
                if (chi2n < chi2s)
                {
                    parameters = newparams;
                    exp.o2 = no2;
                    exp.o4 = no4;
                    exp.o6 = no6;
                    lambda *= 1.1;
                }
                else
                {
                    lambda *= 1 / 1.1;
                }
            }
            Console.WriteLine("Fitting finished");
            MSG += "Fitting finished \r\n";
            DumpFitResult(exp, ref MSG, ref LATEX, sw, Hessiandiag, ref sumfexp, ref sumdfexp, chi2n);
        }

        private void DumpFitResult(Multiplet exp, ref string MSG, ref string LATEX, System.Diagnostics.Stopwatch sw, Matrix<double> Hessiandiag, ref double sumfexp, ref double sumdfexp, double chi2n)
        {
            Vector<double> error;
            int size = exp.u2.Count;
            error = (Hessiandiag.Inverse().Diagonal() * chi2n / (size - 3));
            for (int i = 0; i < 3; i++) error[i] = Abs(Sqrt(error[i]));
            Console.WriteLine("Errors " + error[0].ToString("F3") + " " + error[1].ToString("F3") + " " + error[2].ToString("F3"));
            MSG += "Parameters\t o2=" + exp.o2.ToString("G4") + "\t o4=" + exp.o4.ToString("G4") + "\t o6=" + exp.o6.ToString("G4") + "\r\n";
            MSG += "Errors \t do2=" + error[0].ToString("G2") + "\t do4=" + error[1].ToString("G2") + "\t do6=" + error[2].ToString("G2") + "\r\n";
            //LATEX+="$\\Omega_2$ ="+mantissa(o2).ToString("G4")+"$\\cdot$ 10$^{"+exponent(o2)+"}$ $\\Omega_4$ ="+mantissa(o4).ToString("G4")+"$\\cdot$ 10$^{"+exponent(o4)+"}$ $\\Omega_6$ ="+mantissa(o6).ToString("G4")+"$\\cdot$ 10$^{"+exponent(o6)+"}$ \\\\ \r\n";
            LATEX += "$\\Omega_2$ =" + Lf(exp.o2) + " $\\Omega_4$ =" + Lf(exp.o4) + " $\\Omega_6$ =" + Lf(exp.o6) + "\\\\ \r\n";
            LATEX += "$\\Delta\\Omega_2$ =" + Lf(error[0]) + " $\\Delta\\Omega_4$ =" + Lf(error[1]) + " $\\Delta\\Omega_6$ =" + Lf(error[2]) + " \\\\ \r\n";

            Console.WriteLine("Relative errors " + (100 * error[0] / exp.o2).ToString("G3") + "% " + (100 * error[1] / exp.o4).ToString("G3") + "% " + (100 * error[2] / exp.o6).ToString("G3") + "%");
            MSG += "Relative errors \t" + (100 * error[0] / exp.o2).ToString("G3") + "%\t" + (100 * error[1] / exp.o4).ToString("G3") + "%\t" + (100 * error[2] / exp.o6).ToString("G3") + "% \r\n";
            LATEX += "$\\frac{\\Delta\\Omega_2}{\\Omega_2}=$" + (100 * error[0] / exp.o2).ToString("G3") + "\\% $\\frac{\\Delta\\Omega_4}{\\Omega_4}=$" + (100 * error[1] / exp.o4).ToString("G3") + "\\% $\\frac{\\Delta\\Omega_6}{\\Omega_6}=$" + (100 * error[2] / exp.o6).ToString("G3") + " \\%  \\\\\r\n";
            Console.WriteLine("Effective relative error " + (100 * (error[0] / exp.o2 + error[1] / exp.o4 + error[2] / exp.o6)).ToString("G3") + "%");
            LATEX += "$\\frac{\\Delta f}{f}=$ " + (100 * (error[0] / exp.o2 + error[1] / exp.o4 + error[2] / exp.o6)).ToString("G3") + "\\% \r\n";
            MSG += "Wavenumber\tPexp\tPtheor\t(Pexp-Ptheor)/Pexp \r\n";
            for (int i = 0; i < size; i++)
            {
                double ftheor;
                ftheor = F(exp.u2[i], exp.u4[i], exp.u6[i], exp.lambda0[i], exp.n, exp.TwoJPlusOne, exp.o2, exp.o4, exp.o6);
                sumdfexp += Abs(Pow((exp.fexp[i] - ftheor), 2));
                sumfexp += exp.fexp[i] / size;
                Console.WriteLine(exp.fexp[i].ToString("G3") + " " + ftheor.ToString("G3") + "  " + ((exp.fexp[i] - ftheor) / exp.fexp[i]).ToString("G3"));
                MSG += (1.0 / exp.lambda0[i]).ToString() + "\t" + exp.fexp[i].ToString("G4") + "\t" + ftheor.ToString("G4") + "\t" + (Abs(100 * (exp.fexp[i] - ftheor) / exp.fexp[i])).ToString("G4") + "%\r\n";
            }
            MSG += "-----------------------------------------------\r\n";
            Console.WriteLine("-----------------------------------------------");
            Console.WriteLine("RMS = " + Sqrt(sumdfexp / (size - 3)) + " RMS/avg f " + 100 * Sqrt(sumdfexp / (size - 3)) / sumfexp + "%");
            LATEX += "RMS = " + Lf(Sqrt(sumdfexp / (size - 3))) + "$\\frac{RMS}{\\underline{f}}$= " + (100 * Sqrt(sumdfexp / (size - 3)) / sumfexp).ToString("G3") + "\\%\r\n\r\n";
            MSG += "RMS = " + Sqrt(sumdfexp / (size - 3)).ToString("G4") + " RMS/avg f = " + (100 * Sqrt(sumdfexp / (size - 3)) / sumfexp).ToString("G3") + "% \r\n";
            Console.WriteLine("Fitting procedure took " + sw.ElapsedMilliseconds / 1000.0 + " seconds");
            return;
        }
    }
}
