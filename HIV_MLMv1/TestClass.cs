using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using OdeLibrary;

namespace HIV_MLMv1
{
    class TestClass
    {
        //Heuristic Observations:
        // Somehow, there needs to be some capping method to the virus population. 
        // The amount of strains run out of control, consuming insane amounts of processor power.


        // Actually, contrary what we discussed before, new virus infection amount does matter (slightly, perhaps).
        // Effectively, increasing the amount of virus which mutates and the mutation rate can have a very similar effect (up to a certain degree)

        //Simulation parameters, with N being the starting size of the population.


        static void Main(string[] args)
        {
            List<List<double>> Debug = new List<List<double>>();
            int timelimit = 10000;
           
            Solver SolveTest = new Solver();
            ODE testODE = new ODE(3);

            SolveTest.StepperCode = StepperTypeCode.RungeKutta4;

            string DirOut = "";
            DirOut = File.ReadAllText("C:\\Users\\lukxi\\source\\repos\\HIV_MLMv1\\ToBeIgnored\\DirOut.txt");

            int t = 0;
            Random x = new Random();
            //TestLoop
            List<double> betasSweep = new List<double>(); //For s tor
            List<double> ParamSweep = new List<double>();
            for(double Param = 0.000001; Param <0.001; Param *= 1.5)
            {
                double LowestPoint = 10000000;
                List<double> VirusBetas = new List<double> {
                0.000005
                 };
                StateType StartingInfected = new StateType { 1000000, 1, 25 };
                Individual test = new Individual(0, -1, StartingInfected, VirusBetas);
                int PrevCount = 0;

                while (t < timelimit)
                {

                    testODE.SetInit(test.InternalState, test.VBetas);         //Set the internals for the ODE
                    SolveTest.Solve(testODE.VirusDynamics, 0, 0.1, 1, IntegrateFunctionTypeCode.Adaptive); //Solve the ODE system for one timestep
                    test.InternalState = testODE.VirusDynamics.InitialConditions;
                    Debug.Add(test.InternalState.ToList());
                    
                    //Triggers if individual goes extinct
                    if (test.ComputedOnce(200, 0, x, 0.000005, 1, testODE.GetVirusGrowth()))
                    {
                        ParamSweep.Add(testODE.IntHist[0].Min());
                        break;
                    }
                    
                    //Triggers if virus goes extinct by seeing if the sequence of virus history stops keeping track
                    if (PrevCount == testODE.IntHist[0].Count)
                    {
                        ParamSweep.Add(testODE.IntHist[0].Min());
                        break;
                    }
                
                    PrevCount = testODE.IntHist[0].Count;
                    

                    if (t % 100 == 0)
                    {
                        int pause = 0;
                    }
                    t++;
                }
                betasSweep.Add(Param);
            }
        }
    }
}
