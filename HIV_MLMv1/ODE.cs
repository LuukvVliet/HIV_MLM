using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OdeLibrary;

namespace HIV_MLMv1
{
    class ODE
    {
        public LambdaOde VirusDynamics;
        public List<StateType> History = new List<StateType>();
        public List<double> VirusBetas;
        public double virusGrowth;
        public ODE()
        {
            virusGrowth = 0;

            VirusDynamics = new LambdaOde
            {
                
                OdeObserver = (x, t) => { History.Add(x); },
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
                    foreach (double beta in VirusBetas)
                    {
                        sum += x[c];
                        double netGrowth = beta * x[0] * x[c];
                        virusGrowth += netGrowth;

                        PerCapitaDeathToVirusTCells += beta * x[c];
                        dxdt[c] = netGrowth - delta * x[1] * x[c];
                        c++;
                    }

                    dxdt[0] = r * x[0] - d1 * x[0] - (x[0] * sum) / h1 - x[0] * PerCapitaDeathToVirusTCells; // T cells
                    dxdt[1] = (a * x[1] * sum) / (h2 + x[1] + sum) - d2 * x[1];  // Effector cells
                }
            };
        }
        public void SetInit(StateType x, List<double> betas)
        {
            VirusDynamics.InitialConditions = x;
            VirusBetas = betas;
            virusGrowth = 0;
        }
        public int GetVirusGrowth()
        {
            int growth = 0;
            growth = (int)virusGrowth;
            virusGrowth = 0;
            return growth;
        }
    }
}
