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

        //This piece of code generates a List of Lists, of which the first list is the used betas in the parameter sweep.
        //Each subsequent list then has the used Source as their first value
        //Each other value is either the lowest point of the viral load, -1 if the virus went extinct due to oscillatory behaviour 
        //and -2 if the individual went extinct.

        static void Main(string[] args)
        {
            int timelimit = 20000;

            Solver SolveTest = new Solver();

            SolveTest.StepperCode = StepperTypeCode.RungeKutta4;
            string DirOut = "";
           

            List<string> cond_threshold = new List<string>();

            for (int threshold = 20000; threshold <= 20000; threshold += 1000)
            {
                DirOut = "C:\\Users\\lukxi\\source\\repos\\HIV_MLMv1\\ToBeIgnored\\";
                cond_threshold.Add(threshold.ToString());
                StringBuilder build = new StringBuilder();
                for (int testrun = 0; testrun < 50; testrun++)
                {
                    int t = 0;
                    Random x = new Random();
                    ODE testODE = new ODE(3);
                    List<double> VirusBetas = new List<double> {
                    0.00000400005
                 };
                    StateType StartingInfected = new StateType { 1000000, 1, 25, 0 };
                    Individual test = new Individual(0, -1, StartingInfected, VirusBetas);
                    while (t < timelimit)
                    {

                        testODE.SetInit(test.InternalState, test.VBetas);                                       //Set the internals for the ODE
                        SolveTest.Solve(testODE.VirusDynamics, 0, 0.05, 1, IntegrateFunctionTypeCode.Adaptive); //Solve the ODE system for one timestep
                        test.InternalState = testODE.VirusDynamics.InitialConditions;                           //Retrieve the now solved internals

                        //Triggers if individual goes extinct
                        if (test.ComputedOnce(threshold, 0.00001, x, 0.000004, 25, testODE.GetVirusGrowth()))
                        {
                            build.Append(t.ToString());
                            build.AppendLine(",");
                            break;
                        }
                        else if (t == 19999)
                        {
                            build.Append(t.ToString());
                            build.AppendLine(",");
                            break;
                        }
                        t++;
                    }
                }
                DirOut += "Threshold" + threshold.ToString() +".txt";
                build = build.Remove(build.Length - 3, 3);
                using (StreamWriter fs = new StreamWriter(DirOut))
                {
                    fs.Write(build);
                }
            }
            /*for (int testrun = 0; testrun < 100; testrun++)
            {
               
                int t = 0;
                Random x = new Random();
                ODE testODE = new ODE(3);
                List<double> VirusBetas = new List<double> {
                    0.00000400005
                 };
                StateType StartingInfected = new StateType { 1000000, 1, 25, 0 };
                Individual test = new Individual(0, -1, StartingInfected, VirusBetas);
                while (t < timelimit)
                {

                    testODE.SetInit(test.InternalState, test.VBetas);                                       //Set the internals for the ODE
                    SolveTest.Solve(testODE.VirusDynamics, 0, 0.05, 1, IntegrateFunctionTypeCode.Adaptive); //Solve the ODE system for one timestep
                    test.InternalState = testODE.VirusDynamics.InitialConditions;                           //Retrieve the now solved internals

                    //Triggers if individual goes extinct
                    if (test.ComputedOnce(10000, 0.00001, x, 0.000004, 25, testODE.GetVirusGrowth()))
                    {
                        build.Append(t.ToString());
                        build.AppendLine(",");
                        break;
                    }
                    t++;
                    if (t == 19999)
                    {
                        build.Append(t.ToString());
                        build.AppendLine(",");
                        break;
                    }
                }
            }
            int inspect = 0;

            build = build.Remove(build.Length-3, 3);

            DirOut
            using(StreamWriter fs = new StreamWriter(DirOut))
            {
                fs.Write(build);
            } */
            /*bool BetaAdded = false;
            List<double> BetaList = new List<double>();
            //TestLoop(s)
            List<List<double>> externalSweep = new List<List<double>>(); 
            for (double sourceP = 50; sourceP < 100000; sourceP *= 1.5)
            {
                List<double> internalSweep = new List<double>();
                for (double Beta = 0.000001; Beta < 0.001; Beta *= 1.5)
                {
                    BetaList.Add(Beta);
                    
                    List<double> VirusBetas = new List<double> {
                    Beta
                 };
                    StateType StartingInfected = new StateType { 1000000, 1, 25 };
                    Individual test = new Individual(0, -1, StartingInfected, VirusBetas);
                    int PrevCount = 0;

                    while (t < timelimit)
                    {

                        testODE.SetInit(test.InternalState, test.VBetas);                                       //Set the internals for the ODE
                        SolveTest.Solve(testODE.VirusDynamics, 0, 0.05, 1, IntegrateFunctionTypeCode.Adaptive); //Solve the ODE system for one timestep
                        test.InternalState = testODE.VirusDynamics.InitialConditions;

                        //Triggers if individual goes extinct
                        if (test.ComputedOnce(200, 0, x, 0.000005, 1, testODE.GetVirusGrowth()))
                        {
                            internalSweep.Add(-2);
                            break;
                        }

                        //Triggers if virus goes extinct by seeing if the sequence of virus history stops keeping track
                        if (PrevCount == testODE.IntHist[0].Count)
                        {

                            internalSweep.Add(-1);
                            break;
                        }

                        PrevCount = testODE.IntHist[0].Count;


                        if (t % 100 == 0)
                        {
                            int pause = 0;
                        }
                        t++;
                        
                        if(t == timelimit-1)
                            internalSweep.Add(testODE.IntHist[0].Min());
                    }
                    t = 0;
                }
                internalSweep.Prepend(sourceP);
                externalSweep.Add(internalSweep);
            }*/
        }
    }
}
