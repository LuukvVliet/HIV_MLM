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
        //This piece of code generates a List of Lists, of which the first list is the used betas in the parameter sweep.
        //Each subsequent list then has the used Source as their first value
        //Each other value is either the lowest point of the viral load, -1 if the virus went extinct due to oscillatory behaviour 
        //and -2 if the individual went extinct.

        static void Main1(string[] args)
        {

            List<double> betasList = new List<double>();
            
            const double mr = 0.000001;

            for(int counters = -2; counters <98; counters++)
                betasList.Add(0.000002 + 0.0000005 *counters);


            int timelimit = 200;

            Solver SolveTest = new Solver();

            SolveTest.StepperCode = StepperTypeCode.RungeKutta4;

            StringBuilder metaFile = new StringBuilder();
            for (double xAxis = 1; xAxis <= 1; xAxis += 0.01)
            {
                xAxis = Math.Round(xAxis, 2);
                string DirOut = "C:\\Users\\lukxi\\source\\repos\\HIV_MLMv1\\ToBeIgnored\\ReplicationRate_Individual\\";
                DirOut += "Surviving" + "\\";
                Directory.CreateDirectory(DirOut);
                //Console.WriteLine(xAxis);
                using (StreamWriter fs = new StreamWriter(DirOut + "betasweep.txt"))
                {
                    fs.WriteLine("#mr: " + mr + " timelimit: "+ timelimit);
                    fs.Flush();
                    for (int yAxis = 0; yAxis <= 99; yAxis += 1)
                    {
                    int startbeta = yAxis;
                    //yAxis = Math.Round(yAxis, 2);

                    for (int testrun = 0; testrun < 1; testrun++)
                    {
                        List<List<double>> trials = new List<List<double>>();
                        List<double> seenbetas = new List<double>();


                        int t = 0;
                        Random x = new Random();
                        ODE testODE = new ODE(2, betasList);
                        testODE.sourceBase = 5000;
                        // testODE.AttackRate = xAxis;
                        // testODE.ImmuneRecog = yAxis;

                        List<int> VirusBetas = new List<int> {
                        startbeta
                        };
                        StateType StartingInfected = new StateType { 1000000, 1, 25, 0 };
                        Individual test = new Individual(0, -1, StartingInfected, VirusBetas);
                        
                            while (t < timelimit)
                            {

                                testODE.SetInit(test.InternalState, test.IntBetas);                                       //Set the internals for the ODE
                                SolveTest.Solve(testODE.VirusDynamics, 0, 0.05, 1, IntegrateFunctionTypeCode.Adaptive); //Solve the ODE system for one timestep
                                test.InternalState = testODE.VirusDynamics.InitialConditions;                           //Retrieve the now solved internals

                                //Triggers if individual goes extinct
                                if (test.ComputedOnce(50000, mr, x, 25, testODE.GetVirusGrowth()))
                                {

                                    break;
                                }
                                t++;
                            }
                            if (true)
                            {
                                double xx = 0;
                                for (int virusstrain = 0; virusstrain < test.BetasHistory[test.StateHistory.Count-1].Count; virusstrain++)
                                {
                                    xx += test.StateHistory[test.StateHistory.Count - 1][2 + 2 * virusstrain];
                                }
                                fs.WriteLine(betasList[yAxis] + " " + test.StateHistory[test.StateHistory.Count-1][0] + " " + xx + " " + test.expectedMutations[0] +" "+ test.expectedMutations[1]);
                                fs.Flush();
                                
                            }

                            if (false)
                            {

                                trials.Add(new List<double>()); trials.Add(new List<double>()); trials.Add(new List<double>()); trials.Add(new List<double>()); //Adding the first 4 columns
                                foreach (double b in betasList)  //Adding the other columns
                                {
                                    List<double> vlist = new List<double>();
                                    vlist.AddRange(Enumerable.Repeat(0.0, t));
                                    trials.Add(vlist);

                                    List<double> llist = new List<double>();
                                    llist.AddRange(Enumerable.Repeat(0.0, t));
                                    trials.Add(llist);
                                }
                                for (int timepoint = 0; timepoint < t; timepoint++)
                                {
                                    trials[0].Add(test.StateHistory[timepoint][0]);
                                    trials[1].Add(test.StateHistory[timepoint][1]);

                                    int odd = 0;
                                    double a = 0;
                                    double b = 0;
                                    //This foreach loop calculates the total viral load and latent load for one timepoint
                                    foreach(double xxx in test.StateHistory[timepoint].Skip(2))
                                    {

                                        if(odd%2==0)
                                            a += xxx;
                                        else
                                            b += xxx;
                                        odd++;
                                    }
                                    trials[2].Add(a);
                                    trials[3].Add(b);

                                    //trials[2].Add(test.StateHistory[timepoint].Skip(2).Where((k, i) => i % 2 == 0).Sum());
                                    //trials[3].Add(test.StateHistory[timepoint].Skip(2).Where((k, i) => i % 2 == 1).Sum());
                                    for (int virusstrain = 0; virusstrain < test.BetasHistory[timepoint].Count; virusstrain++)
                                    {
                                        int xLoc = (test.BetasHistory[timepoint][virusstrain]*2) + 4;
                                        trials[xLoc][timepoint] = test.StateHistory[timepoint][2 + 2 * virusstrain];
                                        trials[xLoc+1][timepoint] = test.StateHistory[timepoint][3 + 2 * virusstrain];
                                    }
                                }
                                //Add metadata onto files
                                fs.WriteLine("#Source: "+ yAxis.ToString() + "mr: " + mr.ToString());
                                fs.Flush();

                                //Add columnnames
                                StringBuilder build = new StringBuilder();
                                build.Append("time T E I_tot L_tot");
                                for (int col = 0; col < betasList.Count; col++)
                                {
                                    build.Append(' ');
                                    build.Append("I_" + Math.Round(betasList[col], 8).ToString());
                                    build.Append(' ');
                                    build.Append("L_" + Math.Round(betasList[col], 8).ToString());
                                }
                                fs.WriteLine(build);
                                fs.Flush();

                                //Write the file
                                for (int row = 0; row < t; row++)
                                {
                                    build = new StringBuilder();
                                    build.Append(row.ToString() + ' ');
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
