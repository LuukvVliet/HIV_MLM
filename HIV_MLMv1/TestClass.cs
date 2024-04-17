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
        //Heuristic Observations:
        // Somehow, there needs to be some capping method to the virus population. 
        // The amount of strains run out of control, consuming insane amounts of processor power.

        // Solvers dont produce the same results?!

        // Actually, contrary what we discussed before, new virus infection amount does matter (slightly, perhaps).
        // Effectively, increasing the amount of virus which mutates and the mutation rate can have a very similar effect (up to a certain degree)

        //Simulation parameters, with N being the starting size of the population.


        static void Main(string[] args)
        {
            List<List<double>> Debug = new List<List<double>>();
            int timelimit = 10000;
            List<double> VirusBetas = new List<double> {
                0.000005
            };
            StateType StartingInfected = new StateType { 1000000, 1, 25 };
            Individual test = new Individual(0, -1, StartingInfected, VirusBetas);
            Solver SolveTest = new Solver();
            ODE testODE = new ODE(2);

            SolveTest.StepperCode = StepperTypeCode.RungeKutta4;

            int t = 0;
            Random x = new Random();
            while (t < timelimit)
            {

                testODE.SetInit(test.InternalState, test.VBetas);         //Set the internals for the ODE
                SolveTest.Solve(testODE.VirusDynamics, 0, 0.1, 1, IntegrateFunctionTypeCode.Adaptive); //Solve the ODE system for one timestep
                test.InternalState = testODE.VirusDynamics.InitialConditions;
                Debug.Add(test.InternalState.ToList());

                if (test.ComputedOnce(200, 0, x, 0.000005, 1, testODE.GetVirusGrowth()))
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
