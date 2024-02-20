using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OdeLibrary;

namespace HIV_MLMv1
{
    class Individual
    {
        public LambdaOde VirusDynamics { get; set; }
        int ID { get; set; }
        int externalTime;
        List<Tuple<double, double>> VirusDistribution { get; set; } // First variable is amount, second variable is beta.
        List<StateType> History { get; }
        List<List<double>> StateHistory;
        List<List<Tuple<double, double>>> VirusHistory { get; }
        public Individual(int tiem, int id, StateType Init, List<Tuple<double, double>> VD) 
        {
            ID = id;
            externalTime = tiem;
            History = new List<StateType>();
            StateHistory = new List<List<double>> { new List<double>(), new List<double>()};
            VirusHistory = new List<List<Tuple<double, double>>>();


            VirusDistribution = VD;

            VirusDynamics = new LambdaOde
            {
                InitialConditions = Init, // Currently arbitrary initial conditions
                OdeObserver = (x, t) => { History.Add(x); StateHistory[0].Add(x[0]); StateHistory[1].Add(x[1]); },
                OdeSystem = (x, dxdt, t) =>
                {
                    //Parameters currently taken from own R script (death rate from paper, replication too i believe. h1 determined)
                    //Paper supplied by rob has no (apparently) usefull parameters (?).
                    const double d1 = 0.01;
                    const double d2 = 0.01;
                    const double r = 0.111;
                    const double delta = 0.5;
                    const double a = 0.1;
                    const double h1 = 1.09901 * 1000000;
                    const double h2 = 100000;
                    double sum = 0;
                    double PerCapitaDeathToVirusTCells = 0;

                    // T cells are dxdt[0], Effector cells are dxdt[1], rest is Virus cells.
                    int c = 2;
                    foreach (Tuple<double, double> virus in VirusDistribution)
                    {
                        sum += x[c];
                        PerCapitaDeathToVirusTCells += virus.Item2 * x[c];
                        dxdt[c] = virus.Item2 * x[0] * x[c] - delta * x[1] * x[c];
                        c++;
                    }

                    dxdt[0] = r * x[0] - d1 * x[0] - (x[0] * sum) / h1 - x[0] * PerCapitaDeathToVirusTCells; // T cells
                    dxdt[1] = (a * x[1] * sum) / (h2 + x[1] + sum) - d2 * x[1];  // Effector cells
                }
            };
        }
        //Checks for both death and mutation.
        public bool ComputedOnce(double tcellCutoff, double mr, Random rGen, double jumplimit, double newVirusAmount)
        {
            return ComputedOnce(tcellCutoff, mr, rGen, 99999, jumplimit, newVirusAmount);
        }
        public bool ComputedOnce(double tcellCutoff, double mr, Random rGen, int VirusLimit, double jumplimit, double newVirusAmount)
        {
            //Checks to see if the individual still has enough T cells to continue living: returns true if not.
            
            StateType LS = History.Last();
            if (LS[0] <= tcellCutoff)
                return true;
            //Updates the VirusDistribution 
            //If a viral strain has more than 10% of the virus amount with which it usually infects, go to the next population.
            List<Tuple<double, double>> tempVirusDistribution = new List<Tuple<double, double>>();
            StateType NewLS = new StateType();
            double cutoff = 0.1 * newVirusAmount;
            for (int i = 0; i < LS.Count; i++)
            {
                if (i < 2)
                    NewLS.Add(LS[i]);
                //Console.WriteLine(LS[i]);
                else if (LS[i] > cutoff)
                {
                    tempVirusDistribution.Add(new Tuple<double, double>(LS[i], VirusDistribution[i - 2].Item2));
                    NewLS.Add(LS[i]);
                }   
            }
            VirusDistribution = tempVirusDistribution;
            LS = NewLS;
            //Weighted mutation chance, 
            if (VirusDistribution.Count < VirusLimit && rGen.NextDouble() <= mr)
            {
                double sum = VirusDistribution.Sum(x => x.Item1);
                double targetmutation = rGen.Next(0, (int)sum);
                double target = 0;
                for (int i = 0; i < VirusDistribution.Count; i++) {
                    if (VirusDistribution[i].Item1 <= targetmutation)
                        targetmutation -= VirusDistribution[i].Item1;
                    else
                    {
                        target = VirusDistribution[i].Item2;
                        //Remove virus which mutates, up to 'newVirusAmount'
                        if (newVirusAmount - VirusDistribution[i].Item1 < 0)
                            newVirusAmount = VirusDistribution[i].Item1;
                        VirusDistribution[i] = new Tuple<double, double>(VirusDistribution[i].Item1 - newVirusAmount, VirusDistribution[i].Item2);
                        break;
                    }
                }
                target += ((double)(rGen.Next(1, 100) - rGen.Next(1, 100)))/100*jumplimit;
                if (target < 0) target = 0; //dont have a negative beta
                VirusDistribution.Add(new Tuple<double, double>(newVirusAmount, target)); // Adds new virus to virusdistribution
                LS.Add(VirusDistribution.Last().Item1); //Adds the virus to the initial conditions
            }
            VirusHistory.Add(VirusDistribution);
            return false;
        }
        public double NewInfection(Random rGen)
        {
            double result = 0;
            double sum = VirusDistribution.Sum(x => x.Item1);
            double targetmutation = rGen.Next(0, (int)sum);
            for (int i = 0; i < VirusDistribution.Count; i++)
            {
                if (VirusDistribution[i].Item1 <= targetmutation)
                    targetmutation -= VirusDistribution[i].Item1;
                else
                {
                    result = VirusDistribution[i].Item2;
                    break;
                }
            }
            return result;
        }
    }
}
