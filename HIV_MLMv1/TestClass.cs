using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OdeLibrary;

namespace HIV_MLMv1
{
    class TestClass
    {
        
        //Simulation parameters, with N being the starting size of the population.
        static void Main(string[] args)
        {
            int timelimit = 10000;
            List<Tuple<double, double>> VirusDistribution = new List<Tuple<double, double>> {
                //new Tuple<double, double>(500,0.0001),
                new Tuple<double, double>(100,0.00011)
                };
            StateType StartingInfected = new StateType {100000, 10, 100 };
            Individual test = new Individual(0, -1, StartingInfected, VirusDistribution);
            Solver SolveTest = new Solver();
            
            SolveTest.StepperCode = StepperTypeCode.RungeKutta4;

            int t = 0;
            Random x = new Random();
            while(t < timelimit)
            {
                SolveTest.Solve(test.VirusDynamics, 0, 0.001, 1);

                if (test.ComputedOnce(50, 0.5, x, 0.0001, 1))
                    break;
                if (t % 100 == 0) { 
                    int pause = 0; 
                }
                t++;   
            }
            
            double i = Math.Pow(3, 2);
            i = i;
        }
    }
}
