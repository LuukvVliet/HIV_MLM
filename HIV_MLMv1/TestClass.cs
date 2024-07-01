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
            const double jump = 0.0000025;
            const double startbeta = 0.00000400005;
            const double mr = 0.000001;

            List<double> seenbetas = new List<double>();
            List<List<double>> trials = new List<List<double>>();
            
            int timelimit = 1000;

            Solver SolveTest = new Solver();

            SolveTest.StepperCode = StepperTypeCode.RungeKutta4;
            string DirOut = "C:\\Users\\lukxi\\source\\repos\\HIV_MLMv1\\ToBeIgnored\\";

            List<string> cond_threshold = new List<string>();
            for (int threshold = 10000; threshold <= 10000; threshold += 1000)
            {
                cond_threshold.Add(threshold.ToString());
                for (int testrun = 0; testrun < 1; testrun++)
                {
                    int t = 0;
                    Random x = new Random();
                    ODE testODE = new ODE(3);
                    List<double> VirusBetas = new List<double> {
                        startbeta
                 };
                    StateType StartingInfected = new StateType { 1000000, 1, 25, 0 };
                    Individual test = new Individual(0, -1, StartingInfected, VirusBetas);
                    while (t < timelimit)
                    {

                        testODE.SetInit(test.InternalState, test.VBetas);                                       //Set the internals for the ODE
                        SolveTest.Solve(testODE.VirusDynamics, 0, 0.05, 1, IntegrateFunctionTypeCode.Adaptive); //Solve the ODE system for one timestep
                        test.InternalState = testODE.VirusDynamics.InitialConditions;                           //Retrieve the now solved internals

                        //Triggers if individual goes extinct
                        if (test.ComputedOnce(threshold, mr , x, jump, 25, testODE.GetVirusGrowth()))
                        {
                          //  build.Append(t.ToString());
                          //  build.AppendLine(",");
                            break;
                        }
                        else if (t == 19999)
                        {
                           // build.Append(t.ToString());
                           // build.AppendLine(",");
                            break;
                        }
                        if(t%100 == 0)
                        { 
                            int paus = 1; 
                        }
                        t++;
                    }
                    if (true)
                    {
                        trials.Add(new List<double>()); trials.Add(new List<double>()); trials.Add(new List<double>()); //Adding the first 3 columns
                        for (int timepoint = 0; timepoint < t; timepoint++)
                        {
                            trials[0].Add(test.StateHistory[timepoint][0]);
                            trials[1].Add(test.StateHistory[timepoint][1]);
                            trials[2].Add(test.StateHistory[timepoint].Skip(2).Where((k,i) => i%2 == 0).Sum());
                            for (int virusstrain = 0; virusstrain < test.BetasHistory[timepoint].Count; virusstrain++)
                            {
                                int xLoc = (int)((test.BetasHistory[timepoint][virusstrain] - startbeta + 5*jump) / jump); // 5 * jump to determine the location in the table. If one adds more columns to the left, add 1 to this 4 for each column to prevent overlapping.
                                while(trials.Count < xLoc + 1)
                                {
                                    List<double> vlist = new List<double>();
                                    vlist.AddRange(Enumerable.Repeat(0.0, t));
                                    trials.Add(vlist);
                                }
                                trials[xLoc][timepoint] = test.StateHistory[timepoint][2+2 * virusstrain];
                            }
                        }
                    }
                }
                using (StreamWriter fs = new StreamWriter(DirOut + "Traj"+"A"+".txt"))
                {
                    //Add metadata
                    fs.WriteLine("#Threshold:"+threshold.ToString() + "-jump:" + jump.ToString() + "-mr:" + mr.ToString());
                    fs.Flush();

                    //Add columnnames
                    StringBuilder build = new StringBuilder();
                    build.Append("T E I_tot");
                    for (int col = -1; col < trials.Count-4; col++)
                    {
                        build.Append(' ');
                        build.Append("I_"+(startbeta + col*jump).ToString());
                    }
                    fs.WriteLine(build);
                    fs.Flush();

                    //Write the file
                    for (int row = 0; row < timelimit; row++)
                    {
                        build = new StringBuilder();
                        foreach (List<double> column in trials)
                        {
                            /*for (int before = 0; before < row; before++)
                                build.Append("0 ");*/
                            build.Append(column[row].ToString());
                            build.Append(' ');
                        }
                        build = build.Remove(build.Length - 1, 1);
                        fs.WriteLine(build);
                        fs.Flush();
                    }
                }
                
            }
            /* 
             StreamWriter final = new StreamWriter(DirOut + "Summary.txt");
             final.WriteLine("Threshold Mean STD");
             final.Flush();
             string print = "";
             for(int x = 1000; x <= 20000; x+= 1000)
             {
                 StreamReader aa = new StreamReader(DirOut+ "Threshold" + x.ToString() + ".txt");
                 string input = aa.ReadToEnd();
                 string[] values = input.Split(',');
                 int[] v = values.Select(int.Parse).ToArray();
                 double mean = v.Sum() / v.Length;
                 double variance = v.Select(x => Math.Pow(((double)x - mean), 2)).Sum() / v.Length;
                 double STD = Math.Sqrt(variance);

                 final.WriteLine(x.ToString() + " " + mean.ToString() + " " + STD.ToString());
                 final.Flush();
             }*/




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
