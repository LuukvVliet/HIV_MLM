using OdeLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace HIV_MLMv1
{
    class Program
    {
        //Simulation parameters, feel free to change
        public static double vir50 = 13938; 
        public static double virMax = 0.85; //0.317
        public static double virHill = 1.02;
        public static List<double> betasVector = new List<double>();

        public static double DeathProbability = 0.00005;          //Base Deathrate of Individuals
        public static double MutationProbability = 0.00000005;    //Mutation chance of virus; should be scaled to amount of virus

        public static double TCellCutoff = 20000;             //At what threshold does an individual not have enough T cells to live anymore?
        public static int N = 50;                          //Starting amount of infected individuals.
        public static double IcellStart = 25;               //Starting amount of infected cells on new mutation. SHOULD BE equal to or higher than VirusAmountOnMutate due to cutoff calculation.

        public static int VirusAmountOnMutate = 25;         //How many infected cells 'transfer' between two strains when one mutates?
        public static int PopulationLimit = 100;            //Maximum limit of hosts in the simulation. Hosts are replaced upon infection.
        public static List<double> MutDistribution = new List<double> {0.85, 0.1, 0.03, 0.02, 0}; // NEEDS TO SUM TO 1!

        //Recommended to leave these parameters be.
        public static int time = 0;
        public static int timeLimit = 400000; 
        public static Solver SolveSim = new Solver();
        public static StepperTypeCode UsedStepFunction = StepperTypeCode.RungeKutta4; // ODE step solver function type

        static void Main(string[] args)
        {
            //Read the outputdirectory
            string DirOut = "C:\\Users\\lukxi\\source\\repos\\HIV_MLMv1\\ToBeIgnored\\TestRunSim_hill";

            //How many replicates do we run?
            for (int trials = 1; trials <= 1; trials++)
            {
                time = 0;
                int nextID = 0; // ID tracking
                                //Initialize the population with a few infected individuals
                List<Individual> Population = new List<Individual> { };
                List<double> Graveyard = new List<double> { 0, 0, 0, 0, 0 };
                List<double> notes = new List<double>();
                // GRAVEYARD FORMAT //
                // Deaths by T-Cell depletion
                // Deaths by chance
                // Average lifespan up till now

                /*for (int b = -2; b < 40; b++)
                    betasVector.Add(0.000002 * Math.Pow(1.1, b)+0.0000005*b);*/

                for (int b = -2; b < 70; b++)
                    betasVector.Add(0.000002 + 0.00000048 * 3 * b);
                for (int b = 0; b < 50; b++)
                    notes.Add(0);

                //for (int b = -2; b < 40; b++)
                //betasVector.Add(0.000004 + 0.0000015*b);
                //With these settings, the endpoint is reached in +/- 30 steps.

                for (int i = 0; i < N; i++)
                {
                    StateType StartingInfected = new StateType { 1000000, 1, 25, 0 };
                    List<int> StartVBetas = new List<int>();
                    StartVBetas.Add(3);
                    Population.Add(new Individual(time, nextID, StartingInfected, StartVBetas, MutDistribution));
                    nextID++;
                }
                //Create the ODE:
                ODE Sim = new ODE(5, betasVector);
                Sim.sourceBase = 5000;
                // Simulation loop
                Random RInt = new Random();
                SolveSim.StepperCode = UsedStepFunction;
                using (StreamWriter fs = new StreamWriter(DirOut + "\\replicate_" + trials + ".txt"))
                {
                    using (StreamWriter snd = new StreamWriter(DirOut + "\\sat+scale2" + trials + ".txt"))
                    {
                        fs.WriteLine("#HEADER: FraqLatent: " + Sim.FractieLatent + " Activatie: " + Sim.FractieActivatie + " mr: " + MutationProbability + " ");
                        fs.WriteLine("Time PopulationCount AvgAmountOfBetas AvgBeta TotalDeath EvoDeath AverageLifetime");

                        snd.WriteLine("#HEADER: mr: " + MutationProbability + " TCellCutoff: " + TCellCutoff + " PopulationLimit: " + PopulationLimit);
                        snd.Flush();
                        while (time < timeLimit && !(Population.Count == 0))
                        {
                            //Virusloads is the list containing every individual's viral load.
                            List<double> VirusLoads = new List<double>();
                            List<Individual> NextPopulation = new List<Individual> { };

                            /*if (time == 500)
                            {
                                NextPopulation.Add(new Individual(time, -9, new StateType { 1000000, 1, 25, 0 }, new List<int> { 3 }, MutDistribution));
                                NextPopulation.Add(new Individual(time, -9, new StateType { 1000000, 1, 25, 0 }, new List<int> { 3 }, MutDistribution));
                                NextPopulation.Add(new Individual(time, -9, new StateType { 1000000, 1, 25, 0 }, new List<int> { 3 }, MutDistribution));
                            }*/

                            foreach (Individual host in Population)
                            {
                                if (host.MutationTimer > 0)
                                {
                                    Sim.SetInit(host.InternalState, host.IntBetas);         //Set the internals for the ODE
                                    SolveSim.Solve(Sim.VirusDynamics, 0, 0.05, 1, IntegrateFunctionTypeCode.Adaptive); //Solve the ODE system for one timestep
                                    host.InternalState = Sim.VirusDynamics.InitialConditions; //Set the individual to the new internal state.
                                }
                                //After something is computed once, see if it mutates or if a virus is outcompeted.
                                if (host.ComputedOnce(TCellCutoff, MutationProbability, RInt, VirusAmountOnMutate, Sim.GetVirusGrowth()))
                                {
                                    Graveyard[0]++;
                                    Graveyard[2] = (Graveyard[2] * (Graveyard[0] + Graveyard[1] - 1) + time - host.externalTime) / (Graveyard[0] + Graveyard[1]);
                                }
                                //If something doesnt die, add it to the next population
                                else if (RInt.NextDouble() > DeathProbability) NextPopulation.Add(host);
                                else { Graveyard[1]++; Graveyard[2] = (Graveyard[2] * (Graveyard[0] + Graveyard[1] - 1) + time - host.externalTime) / (Graveyard[0] + Graveyard[1]); }

                                //New infection using distribution samples
                                double tempvload = 0;
                                for (int temp = 0; temp < host.VirusState.Count; temp += 2)
                                    tempvload += host.VirusState[temp];
                                VirusLoads.Add(tempvload);

                            }

                            

                            int totalVirusLoad = (int)VirusLoads.Sum();
                            double NewInfectionsExpected = 0;
                            List<double> InfectProbs = new List<double>();
                            foreach (double vprob in VirusLoads)
                            {
                                InfectProbs.Add(1000*CalcVirRevised(vprob));
                                
                                NewInfectionsExpected += CalcVirRevised(vprob);
                            }
                            double infectprobstotal = 1000*NewInfectionsExpected;

                            double NewInfections = 0;
                            for (double xuvu = 0; NewInfectionsExpected > 1; NewInfectionsExpected--)
                                NewInfections++;
                            if (RInt.NextDouble() < NewInfectionsExpected)
                                NewInfections++;

                            
                            /*if(time == 1 || time == 26384)
                            { NewInfections = 1; }*/

                            while (NewInfections > 0)
                            {
                                for (int kms = 0; kms < InfectProbs.Count; kms++)
                                {
                                    InfectProbs[kms] = InfectProbs[kms] / infectprobstotal;
                                }


                                //These are the states of the next infection
                                StateType NextInfection = new StateType { 1000000, 1, IcellStart, 0 };
                                List<int> NextBetas = new List<int>();
                                // By weighted draw, determine which individual and what beta jumps to the next individual.
                                double targetIndv = RInt.NextDouble();

                                int tracker = 0;
                                while (tracker < Population.Count)
                                {
                                    //Weighted drawing of which individual to have infect
                                    if (targetIndv > InfectProbs[tracker])
                                    {
                                        targetIndv -= InfectProbs[tracker];
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
                                                //Check for mutation on jump
                                                if(RInt.NextDouble() <= MutationProbability)
                                                {
                                                    int target = Population[tracker].IntBetas[straincounter / 2];
                                                    double whatjump = RInt.NextDouble();
                                                    whatjump = 0.99;
                                                    if (RInt.NextDouble() <= 0)  // Determine whether the next bin or previous bin should be mutated to
                                                    {
                                                        foreach (double pos in MutDistribution)
                                                        {
                                                            if (pos == MutDistribution[MutDistribution.Count - 1])
                                                            { target = RInt.Next(0, 20); break; }

                                                            target++;
                                                            if (pos > whatjump)
                                                                break;
                                                            whatjump -= pos;
                                                        }
                                                        NextBetas.Add(target);
                                                    }
                                                    else
                                                    {
                                                        foreach (double pos in MutDistribution)
                                                        {
                                                            if (pos == MutDistribution[MutDistribution.Count - 1])
                                                            { target = RInt.Next(7, 8); break; }
                                                            target--;
                                                            if (pos > whatjump)
                                                                break;
                                                            whatjump -= pos;
                                                        }
                                                        NextBetas.Add(target);
                                                    }
                                                    if (target < 0) target = 0;
                                                }
                                                else
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
                                double randomattack = 0.5; // + RInt.NextDouble() * 0.75;
                                //New infections override the older infections
                                if (NextPopulation.Count >= PopulationLimit && true)
                                {
                                    int tag = RInt.Next(0, NextPopulation.Count);
                                    NextPopulation[tag] = new Individual(time, nextID++, NextInfection, NextBetas, MutDistribution);
                                    NextPopulation[tag].AttackRate = randomattack;
                                }
                                else
                                {
                                    NextPopulation.Add(new Individual(time, nextID++, NextInfection, NextBetas, MutDistribution));
                                    NextPopulation[NextPopulation.Count-1].AttackRate = randomattack;
                                }
                                NewInfections--;
                            };


                            if (time % 365*4 == 0 && time != 0)
                            {

                                int cnt = 0; int avgBeta = 0;
                                List<int> betaCount = new List<int>();
                                for (int x = 0; x < 100; x++)
                                    betaCount.Add(0);
                                double averageKillingRate = 0;
                                foreach (Individual alive in Population)
                                {
                                    
                                    foreach (int beta in alive.IntBetas)
                                        betaCount[beta]++;
                                    cnt += alive.IntBetas.Count;
                                    avgBeta += alive.IntBetas.Sum() / alive.IntBetas.Count;
                                    averageKillingRate += alive.AttackRate;
                                }
                                averageKillingRate = averageKillingRate / Population.Count;
                                double cnt2 = (double)cnt / (double)Population.Count;
                                double avgBeta2 = (double)avgBeta / (double)Population.Count;
                                Console.WriteLine(avgBeta2 +" / " +averageKillingRate);

                                StringBuilder build2 = new StringBuilder();
                                foreach (double x in notes)
                                { build2.Append(x); build2.Append(" "); }
                                fs.WriteLine(build2);
                                fs.Flush();
                                
                                
                                //fs.WriteLine(time + " " + Population.Count + " " + cnt2 + " " + avgBeta2 + " " + Graveyard[0] + " " + Graveyard[1] + " " + Graveyard[2]);

                                StringBuilder build = new StringBuilder();
                                foreach (int x in betaCount)
                                { build.Append(x); build.Append(" "); }
                                build.Remove(build.Length-1, 1);
                                snd.WriteLine(build);
                                snd.Flush();
                            }
                            
                        
                        time++;
                        Population = NextPopulation;
                            
                        }
                    int FINISH = 0;
                    }
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
    public static double CalcVirRevised(double V)
        {
            double result;

            result = virMax * Math.Pow(V, 0.8) / (15000 + Math.Pow(V, 0.8));
            result = 1 - Math.Pow(1 - result, 0.0027397260274);

            return result;
        }
  }
}
