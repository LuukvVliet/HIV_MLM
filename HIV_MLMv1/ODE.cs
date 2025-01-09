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
        //betasVector is the list of possible betas; betas Number is the position in the betasVector list which the species is located in.
        public List<double> betasVector;
        public LambdaOde VirusDynamics;
        public List<int> betasNumber;
        
        public double virusGrowth;
        public int sourceBase = 20000;
        public double FractieLatent = 0.05;
        public double FractieActivatie = 0.01;
        public double AttackRate = 0.5;
        public double ImmuneRecog = 10000;
        public ODE(int type, List<double> betas)
        {
            betasVector = betas;
            virusGrowth = 0;
            switch (type)
            {
                case 1:
                    //Model WITH replication rate for latent population.
                    VirusDynamics = new LambdaOde
                    {

                        OdeObserver = (x, t) =>
                        {
                        },
                        OdeSystem = (x, dxdt, t) =>
                        {

                            double source = sourceBase;
                            const double K = 1088240;
                            const double r = 0.111;
                            const double d1 = 0.01;
                            const double deltaI = 0.5;
                            double a = AttackRate;
                            const double deltaE = 0.1;
                            double hE = ImmuneRecog;
                            const double eE = 1;
                            const double eI = 0;
                            double fraqL = FractieLatent;
                            double activL = FractieActivatie;

                            double effectorGrowth = 0;
                            double thisVirusGrowth = 0;
                            double virusLoad = 0;
                            double latentLoad = 0;
                            double virusxbeta = 0;
                            for (int i = 2; i < x.Count; i += 2)
                            {
                                virusLoad += x[i];
                                latentLoad += x[i + 1];
                                virusxbeta += x[i] * betasVector[betasNumber[(i / 2 - 1)]];
                            }
                            // T cells are dxdt[0], Effector cells are dxdt[1], rest is Virus cells.
                            //Getting virus state from virusdistribution and putting it back in virusdistribution over and over is done to prevent n^2 time problems
                            int c = 2;
                            foreach (int beta in betasNumber)
                            {

                                double netGrowth = betasVector[beta] * x[0] * x[c] / (1 + eI * x[0] + eI * virusLoad);
                                double netEffectorDeath = x[1] * x[c] / (hE + eE * x[1] + eE * virusLoad);
                                effectorGrowth +=  netEffectorDeath;
                                thisVirusGrowth += netGrowth;


                                dxdt[c] = (1 - fraqL) * netGrowth + activL * x[c + 1] - deltaI * x[c] - a * netEffectorDeath;
                                dxdt[c + 1] = fraqL * netGrowth +r*x[c+1] -x[c+1]*(x[0] + latentLoad) / K - d1 * x[c + 1] - activL * x[c + 1];
                                c += 2;
                            }

                            dxdt[0] = source + r * x[0] - d1 * x[0] - x[0] * (x[0] + latentLoad) / K - thisVirusGrowth; // T cells
                            dxdt[1] = effectorGrowth - deltaE * x[1];  // Effector cells

                            virusGrowth += thisVirusGrowth;
                        }
                    };  
                    break;

                    //Model containing lower levels of source, with a proper replication rate (Replication set to 0.2 for proper 'hill' effect)
                    //Also adds the 'latent load' for density dependence.
                    //MAYBE add the replication rate to the latent cell population aswell, although this does have some problems
                case 2:
                    VirusDynamics = new LambdaOde 
                    {

                        OdeObserver = (x, t) =>
                        {
                        },
                        OdeSystem = (x, dxdt, t) =>
                        {

                            double source = sourceBase;
                            const double K = 1088240;
                            const double r = 0.4;
                            const double d1 = 0.01;
                            const double deltaI = 0.5;
                            double a = AttackRate;
                            const double deltaE = 0.1;
                            double hE = ImmuneRecog;
                            const double eE = 1;
                            const double eI = 0;
                            double fraqL = FractieLatent;
                            double activL = FractieActivatie;

                            double effectorGrowth = 0;
                            double thisVirusGrowth = 0;
                            double virusLoad = 0;
                            double latentLoad = 0;
                            double virusxbeta = 0;
                            for (int i = 2; i < x.Count; i += 2)
                            {
                                virusLoad += x[i];
                                latentLoad += x[i + 1];
                                virusxbeta += x[i] * betasVector[betasNumber[(i / 2 - 1)]];
                            }
                            // T cells are dxdt[0], Effector cells are dxdt[1], rest is Virus cells.
                            //Getting virus state from virusdistribution and putting it back in virusdistribution over and over is done to prevent n^2 time problems
                            int c = 2;
                            foreach (int beta in betasNumber)
                            {

                                double netGrowth = betasVector[beta] * x[0] * x[c] / (1 + eI * x[0] + eI * virusLoad);
                                double netEffectorDeath =  x[1] * x[c] / (hE + eE * x[1] + eE * virusLoad);
                                effectorGrowth +=  netEffectorDeath;
                                thisVirusGrowth += netGrowth;


                                dxdt[c] = (1 - fraqL) * netGrowth + activL * x[c + 1] - deltaI * x[c] - a * netEffectorDeath;
                                dxdt[c + 1] = fraqL * netGrowth - d1 * x[c + 1] - activL * x[c + 1];
                                c += 2;
                            }

                            dxdt[0] = source + r * x[0] - d1 * x[0]-x[0]*(x[0]+latentLoad)/K - thisVirusGrowth; // T cells
                            dxdt[1] = effectorGrowth - deltaE * x[1];  // Effector cells

                            virusGrowth += thisVirusGrowth;
                        }
                    };
                    break;
                case 5:
                    // The original source with no replication rate model
                    VirusDynamics = new LambdaOde
                    {

                        OdeObserver = (x, t) =>
                        {
                        },
                        OdeSystem = (x, dxdt, t) =>
                        {

                            //Parameters currently taken from own R script (death rate from paper, replication too i believe. h1 determined)
                            //Paper supplied by rob has no (apparently) usefull parameters (UNTRUE).
                            double source = 10000;
                            const double r = 0;
                            const double d1 = 0.01;
                            const double deltaI = 0.5;
                            double a = AttackRate;
                            const double deltaE = 0.1;
                            double hE = ImmuneRecog;
                            const double eE = 1;
                            const double eI = 0;
                            double fraqL = FractieLatent;
                            double activL = FractieActivatie;

                            double effectorGrowth = 0;
                            double thisVirusGrowth = 0;
                            double virusLoad = 0;
                            double latentLoad = 0;
                            double virusxbeta = 0;
                            for (int i = 2; i < x.Count; i += 2)
                            {
                                virusLoad += x[i];
                                virusxbeta += x[i] * betasVector[betasNumber[(i / 2 - 1)]];
                            }
                            // T cells are dxdt[0], Effector cells are dxdt[1], rest is Virus cells.
                            //Getting virus state from virusdistribution and putting it back in virusdistribution over and over is done to prevent n^2 time problems
                            int c = 2;
                            foreach (int beta in betasNumber)
                            {

                                double netGrowth = betasVector[beta] * x[0] * x[c] / (1 + eI * x[0] + eI * virusLoad);
                                double netEffectorDeath = x[1] * x[c] / (hE + eE * x[1] + eE * virusLoad);
                                effectorGrowth += netEffectorDeath;
                                thisVirusGrowth += netGrowth;


                                dxdt[c] = (1 - fraqL) * netGrowth + activL * x[c + 1] - deltaI * x[c] - a * netEffectorDeath;
                                dxdt[c + 1] = fraqL * netGrowth - d1 * x[c + 1] - activL * x[c + 1];
                                c += 2;
                            }

                            dxdt[0] = source + r * x[0] - d1 * x[0] - thisVirusGrowth; // T cells
                            dxdt[1] = effectorGrowth - deltaE * x[1];  // Effector cells

                            virusGrowth += thisVirusGrowth;
                        }
                    };
                    break;
                case 3:
                    // The original source with no replication rate model
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
                            double d1 = sourceBase/KT;
                            const double deltaI = 0.5;
                            double a = AttackRate;
                            const double deltaE = 0.1;
                            double hE = ImmuneRecog;
                            const double eE = 1;
                            const double eI = 1/KT;
                            double fraqL = FractieLatent;
                            double activL = FractieActivatie;

                            double effectorGrowth = 0;
                            double thisVirusGrowth = 0;
                            double virusLoad = 0;
                            double virusxbeta = 0;
                            for (int i = 2; i < x.Count; i+=2)
                                {
                                    virusLoad += x[i];
                                    virusxbeta += x[i] * betasVector[betasNumber[(i/2-1)]];
                                }
                            // T cells are dxdt[0], Effector cells are dxdt[1], rest is Virus cells.
                            //Getting virus state from virusdistribution and putting it back in virusdistribution over and over is done to prevent n^2 time problems
                            int c = 2;
                            foreach (int beta in betasNumber)
                            {
                                
                                double netGrowth = betasVector[beta] * x[0] * x[c]/(1+eI*x[0]+eI*virusLoad);
                                double netEffectorDeath = x[1] * x[c] / (hE + eE * x[1] + eE * virusLoad);
                                effectorGrowth += netEffectorDeath;
                                thisVirusGrowth += netGrowth;
                                

                                dxdt[c] = (1-fraqL)*netGrowth + activL*x[c+1] - deltaI*x[c] - a*netEffectorDeath;
                                dxdt[c + 1] = fraqL * netGrowth - d1 * x[c + 1] - activL * x[c + 1];
                                c += 2;
                            }

                            dxdt[0] = source + r * x[0] - d1 * x[0] - thisVirusGrowth; // T cells
                            dxdt[1] = effectorGrowth - deltaE * x[1];  // Effector cells

                            virusGrowth += thisVirusGrowth;
                        }
                    };
                    break;
                    //Model with replication rate, no latent replication, without saturation function.
                case 4:
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
                            const double K = 1088240;
                            const double r = 0.2;
                            const double d1 = 0.01;
                            const double deltaI = 0.5;
                            double a = AttackRate;
                            const double deltaE = 0.1;
                            double hE = ImmuneRecog;
                            const double eE = 1;
                            double fraqL = FractieLatent;
                            double activL = FractieActivatie;

                            double effectorGrowth = 0;
                            double thisVirusGrowth = 0;
                            double virusLoad = 0;
                            double latentLoad = 0;
                            double virusxbeta = 0;
                            for (int i = 2; i < x.Count; i += 2)
                            {
                                virusLoad += x[i];
                                latentLoad += x[i + 1];
                                virusxbeta += x[i] * betasVector[betasNumber[(i / 2 - 1)]];
                            }
                            // T cells are dxdt[0], Effector cells are dxdt[1], rest is Virus cells.
                            //Getting virus state from virusdistribution and putting it back in virusdistribution over and over is done to prevent n^2 time problems
                            int c = 2;
                            foreach (int beta in betasNumber)
                            {

                                double netGrowth = betasVector[beta] * x[0] * x[c];
                                double netEffectorDeath = x[1] * x[c] ;
                                effectorGrowth += netEffectorDeath / (hE + eE * x[1] + eE * virusLoad);
                                thisVirusGrowth += netGrowth;


                                dxdt[c] = (1 - fraqL) * netGrowth + activL * x[c + 1] - deltaI * x[c] - a * netEffectorDeath;
                                dxdt[c + 1] = fraqL * netGrowth - d1 * x[c + 1] - activL * x[c + 1];
                                c += 2;
                            }

                            dxdt[0] = source + r * x[0] - d1 * x[0] - x[0] * (x[0] + latentLoad) / K - thisVirusGrowth; // T cells
                            dxdt[1] = effectorGrowth - deltaE * x[1];  // Effector cells

                            virusGrowth += thisVirusGrowth;
                        }
                    };
                    break;
            }
            
        }
        public void SetInit(StateType x, List<int> betas)
        {
            VirusDynamics.InitialConditions = x;
            betasNumber = betas;
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
