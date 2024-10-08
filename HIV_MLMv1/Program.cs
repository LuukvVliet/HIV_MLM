using OdeLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace HIV_MLMv1
{
    class Program
    {
        //Simulation parameters, feel free to change
        public static double vir50 = 13938; //Check the amount of target cells actually present in 1 ml of blood
        public static double virMax = 0.317;
        public static double virHill = 1.02;
        public static List<double> betasVector = new List<double>();

        public static double DeathProbability = 0.000075;          //Base Deathrate of Individuals
        public static double MutationProbability = 0.000005;    //Mutation chance of virus; should be scaled to amount of virus

        public static double TCellCutoff = 100000;             //At what threshold does an individual not have enough T cells to live anymore?
        public static int N = 20;                          //Starting amount of infected individuals.
        public static double IcellStart = 25;               //Starting amount of infected cells on new mutation. SHOULD BE equal to or higher than VirusAmountOnMutate due to cutoff calculation.

        public static int VirusAmountOnMutate = 25;         //How many infected cells 'transfer' between two strains when one mutates?



        //Recommended to leave these parameters be.
        public static int time = 0;
        public static int timeLimit = 100000; 
        public static Solver SolveSim = new Solver();
        public static StepperTypeCode UsedStepFunction = StepperTypeCode.RungeKutta4; // ODE step solver function type

        static void Main1(string[] args)
        {
            //Read the outputdirectory
            string DirOut = "C:\\Users\\lukxi\\source\\repos\\HIV_MLMv1\\ToBeIgnored\\TestRunSim";

            //How many replicates do we run?
            for (int trials = 1; trials <= 20; trials++)
            {
                time = 0;
                int nextID = 0; // ID tracking
                                //Initialize the population with a few infected individuals
                List<Individual> Population = new List<Individual> { };
                List<double> Graveyard = new List<double> { 0, 0, 0, 0, 0 };
                // GRAVEYARD FORMAT //
                // Deaths by T-Cell depletion
                // Deaths by chance
                // Average lifespan up till now

                for (int b = -2; b < 70; b++)
                    betasVector.Add(0.000004 * Math.Pow(1.03, b));

                for (int i = 0; i < N; i++)
                {
                    StateType StartingInfected = new StateType { 1000000, 1, 25, 0 };
                    List<int> StartVBetas = new List<int>();
                    StartVBetas.Add(2);
                    Population.Add(new Individual(time, nextID, StartingInfected, StartVBetas));
                    nextID++;
                }

                //Create the ODE:
                ODE Sim = new ODE(3, betasVector);

                // Simulation loop
                Random RInt = new Random();
                SolveSim.StepperCode = UsedStepFunction;
                using (StreamWriter fs = new StreamWriter(DirOut + "\\replicate_" + trials + ".txt"))
                {
                    fs.WriteLine("#HEADER: FraqLatent: " + Sim.FractieLatent + " Activatie: " + Sim.FractieActivatie + " mr: " + MutationProbability + " ");
                    fs.WriteLine("Time PopulationCount AvgAmountOfBetas AvgBeta TotalDeath EvoDeath NaturalDeath");

                    while (time < timeLimit && !(Population.Count == 0))
                    {
                        //Virusloads is the list containing every individual's viral load.
                        List<double> VirusLoads = new List<double>();
                        List<Individual> NextPopulation = new List<Individual> { };
                        foreach (Individual x in Population)
                        {
                            Sim.SetInit(x.InternalState, x.IntBetas);         //Set the internals for the ODE
                                                                              //Console.WriteLine(x.InternalState.Count);

                            SolveSim.Solve(Sim.VirusDynamics, 0, 0.05, 1, IntegrateFunctionTypeCode.Adaptive); //Solve the ODE system for one timestep
                            x.InternalState = Sim.VirusDynamics.InitialConditions; //Set the individual to the new internal state.

                            //After something is computed once, see if it mutates or if a virus is outcompeted.
                            if (x.ComputedOnce(TCellCutoff, MutationProbability, RInt, VirusAmountOnMutate, Sim.GetVirusGrowth()))
                            { 
                                Graveyard[0]++; 
                                Graveyard[2] = (Graveyard[2] * (Graveyard[0] + Graveyard[1]-1) + time - x.externalTime) / (Graveyard[0] + Graveyard[1]);
                            }

                            //If something doesnt die, add it to the next population
                            else if (RInt.NextDouble() > DeathProbability) NextPopulation.Add(x);
                            else { Graveyard[1]++; Graveyard[2] = (Graveyard[2] * (Graveyard[0] + Graveyard[1]-1) + time - x.externalTime) / (Graveyard[0] + Graveyard[1]); }

                            //New infection using distribution samples
                            double tempvload = 0;
                            for (int temp = 0; temp < x.VirusState.Count; temp += 2)
                                tempvload += x.VirusState[temp];
                            VirusLoads.Add(tempvload);

                        }



                        int totalVirusLoad = (int)VirusLoads.Sum();
                        double NewInfections = 0;
                        foreach (double vprob in VirusLoads)
                        {

                            NewInfections += CalcVirulence(vprob);
                        }

                        if (NewInfections >= 1)
                            NewInfections = (int)Math.Round(NewInfections);
                        else if (RInt.NextDouble() < NewInfections)
                            NewInfections = 1;
                        else
                            NewInfections = 0;

                        //LIMIT FOR POPULATION SIZE
                        if (Population.Count >= 1000)
                            NewInfections = 0;

                        while (NewInfections > 0)
                        {

                            //These are the states of the next infection
                            StateType NextInfection = new StateType { 1000000, 1, IcellStart, 0 };
                            List<int> NextBetas = new List<int>();
                            // By weighted draw, determine which individual and what beta jumps to the next individual.
                            int targetIndv = RInt.Next(0, totalVirusLoad);
                            double TI = (double)targetIndv;
                            int tracker = 0;
                            while (tracker < Population.Count)
                            {
                                //Weighted drawing of which individual to have infect
                                if (TI > VirusLoads[tracker])
                                {
                                    TI -= VirusLoads[tracker];
                                }
                                else
                                {
                                    //Weighted drawing of which virus to use for infection
                                    int targetVir = RInt.Next(0, (int)VirusLoads[tracker]);
                                    double TV = (double)targetVir - 1;
                                    int straincounter = 0;



                                    foreach (double strain in Population[tracker].VirusState)
                                    {
                                        if (straincounter % 2 == 1) { } //litte bit of a stinky hack, this prevents the straincounter to measure every other input; every other input is the latent population of the state

                                        else if (TV > strain)
                                        {
                                            TV -= strain;
                                        }
                                        else
                                        {
                                            NextBetas.Add(Population[tracker].IntBetas[straincounter / 2]);

                                            break;
                                        }
                                        straincounter++;
                                    }
                                    break;
                                }
                                tracker++;
                            }
                            //Population and Virusloads are not always equal in length
                            NextPopulation.Add(new Individual(time, nextID++, NextInfection, NextBetas));
                            NewInfections--;
                        };


                        if (time % 10 == 0 && time != 0)
                        {
                            int cnt = 0; int avgBeta = 0;
                            foreach (Individual alive in Population)
                            {
                                cnt += alive.IntBetas.Count;
                                avgBeta += alive.IntBetas.Sum() / alive.IntBetas.Count;
                            }
                            double cnt2 = (double)cnt / (double)Population.Count;
                            double avgBeta2 = (double)avgBeta / (double)Population.Count;
                            Console.WriteLine(cnt2 + " / " + avgBeta2);

                            fs.WriteLine(time + " " + Population.Count + " " + cnt2 + " " + avgBeta2 + " " + Graveyard[0] + " " + Graveyard[1] + " " + Graveyard[2]);



                        }
                        time++;
                        Population = NextPopulation;
                    }
                    int FINISH = 0;
                }
            }
            // ADD post-processing of lists
            // Print to .csv type file for further analysis of model function.
        }
    
    public static double CalcVirulence(double V)
    {
        double result;

            V = V * 100; // Modifier for translating virus particles to infective cells 
            result = virMax * Math.Pow(V, virHill) / (Math.Pow(V, virHill) + Math.Pow(vir50, virHill));
            result = 1-Math.Pow(1-result, 0.0027397260274); // The long number shown here is 1/365
        return result;
    }
  }
}
