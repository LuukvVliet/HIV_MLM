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
        public List<double> VirusBetas;
        
        public double virusGrowth;
        const int sourceBase = 20000;
        public ODE(int type)
        {
            virusGrowth = 0;
            switch (type)
            {
                case 1:
                    VirusDynamics = new LambdaOde
                    {

                        OdeObserver = (x, t) => 
                        { 
                           // History.Add(x); IntHist.Add(new List<double> { x[0], x[1], x[2] }); 
                        },
                        OdeSystem = (x, dxdt, t) =>
                        {

                            //Parameters currently taken from own R script (death rate from paper, replication too i believe. h1 determined)
                            //Paper supplied by rob has no (apparently) usefull parameters (UNTRUE).
                            const double r = 0.111;
                            const double d1 = 0.01;
                            const double delta = 0.5;
                            const double h1 = 1.09901 * 1000000;
                            const double a = 1.1;
                            const double d2 = 0.1;
                            const double h2 = 10000;
                            double sum = 0;
                            double PerCapitaDeathToTCells = 0;
                            // T cells are dxdt[0], Effector cells are dxdt[1], rest is Virus cells.
                            //Getting virus state from virusdistribution and putting it back in virusdistribution over and over is done to prevent n^2 time problems
                            int c = 2;
                            foreach (double beta in VirusBetas)
                            {
                                sum += x[c];
                                double netGrowth = beta * x[0] * x[c];
                                virusGrowth += netGrowth;

                                PerCapitaDeathToTCells += beta * x[c];
                                dxdt[c] = netGrowth - delta * x[1] * x[c];
                                c++;
                            }

                            dxdt[0] = r * x[0] - d1 * x[0] - r*x[0]*(x[0] + sum) / h1 - x[0] * PerCapitaDeathToTCells; // T cells
                            dxdt[1] = (a * x[1] * sum) / (h2 + x[1] + sum) - d2 * x[1];  // Effector cells
                        }
                    };
                    break;
                case 2:
                    VirusDynamics = new LambdaOde
                    {

                        OdeObserver = (x, t) =>
                        {
                            // History.Add(x); IntHist.Add(new List<double> { x[0], x[1], x[2] }); 
                        },
                        OdeSystem = (x, dxdt, t) =>
                        {

                            //Parameters currently taken from own R script (death rate from paper, replication too i believe. h1 determined)

                            const double r = 0.111;
                            const double d1 = 0.01;
                            const double delta = 0.5;
                            const double h1 = 1.09901 * 1000000;
                            const double a = 1.1;
                            const double d2 = 0.1;
                            const double h2 = 10000;
                            //   r * x[0] - d1 * x[0] - r * x[0] * (x[0] + x[2]) / h1 - VirusBetas[0] * x[0] * x[2]; // T cells
                            dxdt[0] = x[0] * (r - r * (x[0] + x[2]) / h1 - d1 - VirusBetas[0] * x[2]);
                            dxdt[2] = VirusBetas[0] * x[0] * x[2] - delta * x[1] * x[2]; // Infected cells
                            dxdt[1] = (a * x[1] * x[2]) / (h2 + x[1] + x[2]) - d2 * x[1];  // Effector cells
                        }
                    };
                    break;
                case 3:
                    VirusDynamics = new LambdaOde
                    {

                        OdeObserver = (x, t) =>
                        {
                        },
                        OdeSystem = (x, dxdt, t) =>
                        {

                            //Parameters currently taken from own R script (death rate from paper, replication too i believe. h1 determined)
                            //Paper supplied by rob has no (apparently) usefull parameters (UNTRUE).
                            double source = sourceBase;
                            const double KT = 1000000;
                            const double r = 0;
                            const double d1 = sourceBase/KT;
                            const double deltaI = 0.5;
                            const double k = 1;
                            const double a = 1.1;
                            const double deltaE = 0.1;
                            const double hE = 10000;
                            const double eE = 1;
                            const double eI = 1/KT;
                            const double fraqL = 0.04;
                            const double activL = 0.05;

                            double effectorGrowth = 0;
                            double thisVirusGrowth = 0;

                            // T cells are dxdt[0], Effector cells are dxdt[1], rest is Virus cells.
                            //Getting virus state from virusdistribution and putting it back in virusdistribution over and over is done to prevent n^2 time problems
                            int c = 2;
                            foreach (double beta in VirusBetas)
                            {
                                
                                double netGrowth = beta * x[0] * x[c]/(1+eI*x[0]+eI*x[c]);
                                double netEffectorDeath = x[1] * x[c] / (hE + eE * x[1] + eE * x[c]);
                                effectorGrowth += a*netEffectorDeath;
                                thisVirusGrowth += netGrowth;
                                

                                dxdt[c] = (1-fraqL)*netGrowth + activL*x[c+1] - deltaI*x[c] - k*netEffectorDeath;
                                dxdt[c + 1] = fraqL * netGrowth - d1 * x[c + 1] - activL * x[c + 1];
                                c += 2;
                            }

                            dxdt[0] = source + r * x[0] - d1 * x[0] - thisVirusGrowth; // T cells
                            dxdt[1] = effectorGrowth - deltaE * x[1];  // Effector cells

                            virusGrowth += thisVirusGrowth;
                        }
                    };
                    break;
            }
            
        }
        public void SetInit(StateType x, List<double> betas)
        {
            VirusDynamics.InitialConditions = x;
            VirusBetas = betas;
            virusGrowth = 0;
        }
        //Function used to get the Virus growth over the last period. Resets virus growth back to 0.
        public int GetVirusGrowth()
        {
            int growth = 0;
            growth = (int)virusGrowth;
            if (growth < 0)
                growth = int.MaxValue;
            virusGrowth = 0;
            return growth;
        }
    }
}
