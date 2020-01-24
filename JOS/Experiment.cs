using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using static Globals;

namespace JOS
{
    public class Multiplet
    {

        public List<double> u2;
        public List<double> u4;
        public List<double> u6;
        public List<double> lambda0;
        public List<double> fexp;
        public List<double> Aii;
        public string absofilename;
        public double n;
        public double sumrate = 0;
        public double TwoJPlusOne;
        public Multiplet()
        {
            u2 = new List<double>();
            u4 = new List<double>();
            u6 = new List<double>();
            lambda0 = new List<double>();
            fexp = new List<double>();
            Aii = new List<double>();


        }
        public double o2, o4, o6;
        public void GetOmegas(Multiplet donor)
        {
            o2 = donor.o2;
            o4 = donor.o4;
            o6 = donor.o6;
        }
        public void LoadFromFile(string filename)
        {
            char[] charSeparators = new char[] { ' ', '\t' };
            string[] data = File.ReadAllLines(filename);
            string[] line = data[0].Split(charSeparators);
            TwoJPlusOne = double.Parse(line[0], CultureInfo.InvariantCulture);
            n = double.Parse(line[1], CultureInfo.InvariantCulture);
            absofilename = Path.GetFileNameWithoutExtension(filename);
            Console.WriteLine("n= " + n.ToString("g4", CultureInfo.InvariantCulture) + " 2J+1=" + TwoJPlusOne.ToString("g4", CultureInfo.InvariantCulture));
            for (int i = 1; i < data.Length; i++)
            {
                string[] pline = data[i].Split(charSeparators, StringSplitOptions.RemoveEmptyEntries);
                double tu2, tu4, tu6, tfosc, tlam;
                //foreach (string s in pline )Console.WriteLine(s);
                tfosc = double.Parse(pline[0], CultureInfo.InvariantCulture);
                tlam = 1.0 / double.Parse(pline[1], CultureInfo.InvariantCulture);

                tu2 = double.Parse(pline[2], CultureInfo.InvariantCulture);
                tu4 = double.Parse(pline[3], CultureInfo.InvariantCulture);
                tu6 = double.Parse(pline[4], CultureInfo.InvariantCulture);
                u2.Add(tu2);
                u4.Add(tu4);
                u6.Add(tu6);
                lambda0.Add(tlam);
                fexp.Add(tfosc);
                Console.WriteLine((1e7 * tlam).ToString("G4") + " " + tfosc.ToString("G4") + " " + tu2.ToString("G4") + " " + tu4.ToString("G4") + " " + tu6.ToString("G4"));
            }
        }
        public void LoadEmiFromFile(string filename)
        {
            char[] separator = { ' ', '\t' };
            string[] data = File.ReadAllLines(filename);
            string[] line = data[0].Split(separator);
            TwoJPlusOne = Double.Parse(line[0], CultureInfo.InvariantCulture);
            n = double.Parse(line[1], CultureInfo.InvariantCulture);
            for (int i = 1; i < data.Length; i++)
            {
                string[] lines = data[i].Split(separator);
                lambda0.Add(1.0 / double.Parse(lines[0], CultureInfo.InvariantCulture));
                u2.Add(double.Parse(lines[1], CultureInfo.InvariantCulture));
                u4.Add(double.Parse(lines[2], CultureInfo.InvariantCulture));
                u6.Add(double.Parse(lines[3], CultureInfo.InvariantCulture));

            }
        }
        public void CalculateRates()
        {

            for (int i = 0; i < u2.Count; i++)
            {
                double rate = ((64 * Math.Pow(pi, 4) * qe * qe) / (3 * h * Math.Pow(lambda0[i], 3) * TwoJPlusOne)) * (n * Math.Pow(n * n + 2, 2) / 9) * (u2[i] * o2 + u4[i] * o4 + u6[i] * o6);
                sumrate += rate;
                Aii.Add(rate);
                //  Console.WriteLine(rate);
            }
            //Console.WriteLine(sumrate + " " + 1e6 / sumrate);
        }
        public void ReportRates(out string msg)
        {
            int size = u2.Count;
            msg = "Wavenumber\t Transition rate\tBranching ratio\r\n";

            for (int i = 0; i < size; i++)
            {
                msg += (1.0 / lambda0[i]).ToString("G5") + "\t" + (Aii[i]).ToString("F1") + "\t" + (100 * Aii[i] / sumrate).ToString("F1") + "%\r\n";

            }
            msg += "Summary rate = " + (sumrate).ToString("G5") + " which amounts to lifetime " + (1000.0 / sumrate).ToString("G5") + " ms\r\n";

        }
        public void ReportRates(out string msg, out string latex)
        {
            int size = u2.Count;
            latex = "\r\n \\begin{tabular}{|l|l|l|l|}\r\n";
            latex += "\\hline\r\n";
            latex += "Multiplet $\\rightarrow$ & Wavelenumber [cm$^{-1}$]&Rate&Branching ratio\\cr \r\n";
            latex += "\\hline \r\n";
            msg = "Wavenumber\t Transition rate\tBranching ratio\r\n";

            for (int i = 0; i < size; i++)
            {
                msg += (1.0 / lambda0[i]).ToString("G5") + "\t" + (Aii[i]).ToString("F1") + "\t" + (100 * Aii[i] / sumrate).ToString("F1") + "%\r\n";
                latex += "& " + (1.0 / lambda0[i]).ToString("G5") + " & " + (Aii[i]).ToString("F1") + " & " + (100 * Aii[i] / sumrate).ToString("F1") + "\\% \\cr\r\n";
                latex += "\\hline \r\n";
            }
            msg += "Summary rate = " + (sumrate).ToString("G5") + " which amounts to lifetime " + (1000.0 / sumrate).ToString("G5") + " ms\r\n";
            latex += "\\multicolumn{4}{|c|}{$\\Sigma A$ = " + (sumrate).ToString("G5") + " which is equal to $\\tau$=" + (1000.0 / sumrate).ToString("G5") + " ms}\\cr\r\n";
            latex += "\\hline";
            latex += "\\end{tabular}\r\n";

        }
        public void DumpEmidata(out string msg)
        {
            msg = "";
            msg += "2J + 1 = " + TwoJPlusOne + " n = " + n + "\r\n";
            for (int i = 0; i < u2.Count; i++)
            {
                msg += 1 / lambda0[i] + "\t" + u2[i] + "\t" + u4[i] + "\t" + u6[i] + "\r\n";
            }
        }
    }

}
