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
        public Individual(int tiem, int id, StateType Init) : this(tiem, id, Init, new List<Tuple<double, double>> { new Tuple<double, double>(1000, 0.001) })
        {
            ID = id;
            externalTime = tiem;
            History = new List<StateType>();

            VirusDynamics = new LambdaOde
            {
                
                InitialConditions = Init, // Currently arbitrary initial conditions
                OdeObserver = (x, t) => {History.Add(x);},
                OdeSystem = (x, dxdt, t) => 
                {
                    const double d1 = 0;
                    const double d2 = 0;
                    const double r = 0;
                    double beta1 = 0;
                    double beta2 = 0;
                    const double delta = 0;
                    const double a = 0;
                    const double h1 = 0;
                    const double h2 = 0;
                    double sum = x[2] + x[3];
                    // T cells are dxdt[0], Effector cells are dxdt[1], rest is Virus cells. EXPERIMENTAL STUFF HERE

                    /*int c = 2;
                    foreach (Tuple<double, double> virus in VirusDistribution)
                    {
                        sum += virus.Item1;
                        dxdt[c] = virus.Item2 * x[0] * x[c] - delta * x[1] * x[c];
                        c++;
                    }*/

                    dxdt[0] = r * x[0] - d1*x[0] - (x[0] * sum)/h1 - x[0]*(beta1 * x[2]+beta2 *x[3]); // T cells

                    dxdt[1] = (a * x[1] * sum) / (h2 + x[1] + sum) - d2 * x[1];  // Effector cells

                    dxdt[2] = beta1 * x[0]*x[1] - delta * x[1] * x[2]; // Infected Cells
                    dxdt[3] = beta2 * x[0] * x[1] - delta * x[1] * x[2];

                    
                }
            };
        }
        public Individual(int tiem, int id, StateType Init, List<Tuple<double, double>> VD) 
        {
            
        }
        //Checks for both death and mutation.
        public bool ComputedOnce(double tcellCutoff, double mr, Random rGen, int limit)
        {
            //Checks to see if the individual still has enough T cells to continue living: returns true if not.
            if (VirusDynamics.InitialConditions[0] <= tcellCutoff)
                return true;
            //Updates the VirusDistribution 
            //If a viral strain has more than 10 infected cells, go to the next simulation
            StateType LS = History.Last();
            for(int i = 2; i <= LS.Count; i++)
            {
                if (LS[i] > 10)
                {
                    VirusDistribution[i - 2] = new Tuple<double, double>(LS[i], VirusDistribution[i - 2].Item2);
                    VirusDynamics.InitialConditions.Add(LS[i]);
                }
            }
            //Weighted mutation chance, 
            if (VirusDistribution.Count < limit && rGen.NextDouble() <= mr)
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
                        break;
                    }
                }
                target += ((double)(rGen.Next(1, 100) - rGen.Next(1, 100))) / 1000;
                VirusDistribution.Add(new Tuple<double, double>(100, target)); // Adds new virus to virusdistribution
                VirusDynamics.InitialConditions.Add(VirusDistribution.Last().Item1); //Adds the virus to the initial conditions
            }
            return false;
        }
        public double NewInfection()
        {
           double result = 0;
           return result;
        }
    }
}
