using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Globalization;

namespace JOS
{
    public class Experiment
    {
        public List<double> u2;
        public List<double> u4;
        public List<double> u6;
        public List<double> lambda0;
        public List<double> fexp;
        public double n;
        public double TwoJPlusOne;
        public Experiment() {
            u2 = new List<double>();
            u4 = new List<double>();
            u6 = new List<double>();
            lambda0 = new List<double>();
            fexp = new List<double>();
        }
        public double o2, o4, o6;
        public void LoadFromFile(string filename)
        {
            char[] charSeparators = new char[] { ' ','\t' };
            string[] data = System.IO.File.ReadAllLines(filename);
           string[] line =data[0].Split(charSeparators);
            TwoJPlusOne = double.Parse(line[0], CultureInfo.InvariantCulture);
            n = double.Parse(line[1], CultureInfo.InvariantCulture);
            
            Console.WriteLine("n= "+n.ToString("g4", System.Globalization.CultureInfo.InvariantCulture) +" 2J+1="+TwoJPlusOne.ToString("g4", System.Globalization.CultureInfo.InvariantCulture));
            for (int i = 1; i < data.Length; i++)
            {
                string[] pline = data[i].Split(charSeparators, StringSplitOptions.RemoveEmptyEntries);
                double tu2, tu4, tu6, tfosc, tlam;
                //foreach (string s in pline )Console.WriteLine(s);
                tfosc = double.Parse(pline[0], CultureInfo.InvariantCulture);
                tlam = 1.0/double.Parse(pline[1], CultureInfo.InvariantCulture);
                
                tu2 = double.Parse(pline[2], CultureInfo.InvariantCulture);
                tu4 = double.Parse(pline[3], CultureInfo.InvariantCulture);
                tu6 = double.Parse(pline[4], CultureInfo.InvariantCulture);
                u2.Add(tu2);
                u4.Add(tu4);
                u6.Add(tu6);
                lambda0.Add(tlam);
                fexp.Add(tfosc);
                Console.WriteLine(tlam+" "+tfosc+" "+tu2 + " " + tu4 + " " + tu6);
            }
        }
    }
    
}
