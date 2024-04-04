using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OdeLibrary;
using MathNet.Numerics.Random;
using MathNet.Numerics.Distributions;

namespace HIV_MLMv1
{
    class Individual

        // POSSIBLE OPTIMIZATIONS:
        // VECTOR MULTIPLICATION INSTEAD OF LOOPING OVER LIST
        // PREMPTIVELY CALCULATING BINOMIAL DISTRIBUTIONS

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
                    //The following does not work: Somewhere, internally, a number of states is kept constant to what it was before, causing x to reset to its original size.
                    //THIS:
                    //if (!HasBeenUpdated)
                    //{ x = updatedLS; dxdt.AddRange(updatedRange); 
                    //    HasBeenUpdated = true; }
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
                    //Getting virus state from virusdistribution and putting it back in virusdistribution over and over is done to prevent n^2 time problems
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
        
        //VERY brute-force way of updating the ODE to allow for more states.
        public void ReconstructODE(StateType NewInit)
        {
            this.VirusDynamics.Dispose();
            this.VirusDynamics = new LambdaOde
            {
                InitialConditions = NewInit, // Currently arbitrary initial conditions
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
                    //Getting virus state from virusdistribution and putting it back in virusdistribution over and over is done to prevent n^2 time problems
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
        
        public bool ComputedOnce(double tcellCutoff, double mr, Random rGen, double jumplimit, double newVirusAmount)
        {
            return ComputedOnce(tcellCutoff, mr, rGen, 99999, jumplimit, newVirusAmount);
        }
        public bool ComputedOnce(double tcellCutoff, double mr, Random rGen, int VirusLimit, double jumplimit, double newVirusAmount)
        {
            //If a viral strain has more than 10% of the virus amount with which it usually infects, go to the next population.
            double cutoff = 0.1 * newVirusAmount;
            //SumV is used to calculate probability distributions
            double sumV = 0;
            //Checks to see if the individual still has enough T cells to continue living: returns true if not.
            StateType LS = History.Last();
            if (LS[0] <= tcellCutoff)
                return true;
            
            //Updates the VirusDistribution
            List<Tuple<double, double>> tempVirusDistribution = new List<Tuple<double, double>>();
            StateType NewLS = new StateType();
            NewLS.Add(LS[0]); NewLS.Add(LS[1]);

            int x = 2;
            foreach (Tuple<double, double> virus in this.VirusDistribution)
            {
                if (LS[x] > cutoff)
                {
                    tempVirusDistribution.Add(new Tuple<double, double>(LS[x], VirusDistribution[x - 2].Item2));
                    NewLS.Add(LS[x]);
                    sumV += LS[x];
                    x++;
                }
            }

            //Mutations using a binomial distribution
            int newSum = (int)sumV;
            var normalDist = new Binomial(mr, newSum);
            var sample = normalDist.Sample();

            if (VirusDistribution.Count < VirusLimit)
            {
                for(int trial = 0; trial < sample; trial++) { 
                double targetmutation = rGen.Next(0, newSum);
                double target = 0;
                for (int i = 0; i < VirusDistribution.Count; i++)
                {
                    if (VirusDistribution[i].Item1 <= targetmutation)
                        targetmutation -= VirusDistribution[i].Item1;
                    else
                    {
                        target = VirusDistribution[i].Item2;
                        //Remove virus which mutates, up to 'newVirusAmount'
                        if (newVirusAmount - VirusDistribution[i].Item1 > 0)
                            newVirusAmount = VirusDistribution[i].Item1;
                        tempVirusDistribution[i] = new Tuple<double, double>(VirusDistribution[i].Item1 - newVirusAmount, VirusDistribution[i].Item2);
                        break;
                    }
                }
                target += ((double)(rGen.Next(1, 100) - rGen.Next(1, 100))) / 100 * jumplimit;
                if (target < 0) target = 0; //dont have a negative beta
                tempVirusDistribution.Add(new Tuple<double, double>(newVirusAmount, target)); // Adds new virus to virusdistribution
                NewLS.Add(tempVirusDistribution.Last().Item1); //Adds the virus to the initial conditions
                }
            }
            if (VirusDistribution.Count < VirusLimit && sample > 0)
                ReconstructODE(NewLS);
            VirusDistribution = tempVirusDistribution;
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
