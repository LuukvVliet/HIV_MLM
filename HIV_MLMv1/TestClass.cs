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
            List<double> VirusBetas = new List<double> {
                0.00011
            };
            StateType StartingInfected = new StateType { 100000, 10, 100 };
            Individual test = new Individual(0, -1, StartingInfected, VirusBetas);
            Solver SolveTest = new Solver();
            ODE testODE = new ODE();

            SolveTest.StepperCode = StepperTypeCode.RungeKutta4;

            int t = 0;
            Random x = new Random();
            while (t < timelimit)
            {

                testODE.SetInit(test.InternalState, test.VBetas);         //Set the internals for the ODE
                SolveTest.Solve(testODE.VirusDynamics, 0, 1, 0.01); //Solve the ODE system for one timestep
                test.InternalState = testODE.VirusDynamics.InitialConditions; 

                if (test.ComputedOnce(50, 0.000001, x, 0.0001, 1))
                    break;
                if (t % 100 == 0)
                {
                    int pause = 0;
                }
                t++;
            }
            double i = Math.Pow(3, 2);
            i = i;
        }
    }
}
