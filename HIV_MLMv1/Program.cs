using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OdeLibrary;

namespace HIV_MLMv1
{
    class Program
    {
        //Simulation parameters, feel free to change
        public static double InfectProbability = 0.02;      //Infection of new individuals
        public static double DeathProbability = 0;          //Base Deathrate of Individuals
        public static double MutationProbability = 0.05;    //Mutation chance of virus; should be scaled to amount of virus
        public static double MutationJumpLimit = 0.0001;    //The limit of how much a mutation can increase or decrease the beta parameter 
                                                            //NOTE: Increase and decrease of beta is done by a "normal" distribution (or an approximation) so getting the limit is highly unlikely. (as a matter of fact, it is one in ten thousand)
        public static double TCellCutoff = 100;             //At what threshold does an individual not have enough T cells to live anymore?
        public static int N = 100;                          //Starting amount of infected individuals.
        public static double IcellStart = 25;               //Starting amount of infected cells on new mutation.
        public static double IcellStartBeta = 0.0001;       //Starting beta of the first infective virus.
        public static int VirusCountLimit = 5;              //How many different strains of virus can exist in one individual?
        public static int VirusAmountOnMutate = 25;         //How many infected cells 'transfer' between two strains when one mutates?

        public static StateType StartingInfected = new StateType { 10000, 10, 100 }; //State of the first 5 infected individuals
        

        //Recommended to leave these parameters be.
        public static int time = 0;
        public static int timeLimit = 100000; //int.MaxValue; //timeLimit does (technically) count the timelimit, since the time starts on round 0.
        public static Solver SolveSim = new Solver();
        public static StepperTypeCode UsedStepFunction = StepperTypeCode.RungeKutta4; // ODE step solver function type

        static void Main2(string[] args)
        {
            int nextID = 0; // ID tracking
            int xxx = StartingInfected.Count;
            //Initialize the population with a few infected individuals
            List<Individual> Population = new List<Individual> { };
            List<Tuple<Individual, string>> Graveyard = new List<Tuple<Individual, string>> { };
            
            for(int i = 0; i < N; i++)
            {
                List<Tuple<double, double>> StartingVirusConditions = new List<Tuple<double, double>> { new Tuple<double, double>(IcellStart, IcellStartBeta) };
                Population.Add(new Individual(time, nextID, StartingInfected, StartingVirusConditions));
                nextID++;
            }

            // Simulation loop
            Random RInt = new Random();
            SolveSim.StepperCode = UsedStepFunction;
            while (time < timeLimit || !(Population.Count == 0))
            {
                List<Individual> NextPopulation = new List<Individual> { };
                foreach (Individual x in Population)
                {
                    SolveSim.Solve(x.VirusDynamics, 0 , 1, 0.01); //Solve the ODE system for one timestep
                    //If something doesnt die, add it to the next population
                    if (RInt.NextDouble() > DeathProbability) NextPopulation.Add(x);
                    else Graveyard.Add(new Tuple<Individual, string>(x, "Natural death"));
                    
                    //New infection
                    if (RInt.NextDouble() <= InfectProbability) 
                    {
                        StateType NextInfection = new StateType { StartingInfected[0], StartingInfected[1], IcellStart };
                        NextPopulation.Add(new Individual(time, nextID++, NextInfection, new List<Tuple<double, double>> { new Tuple<double, double>(IcellStart, x.NewInfection(RInt)) }));
                    }; 
                    
                    //After something is computed once, see if it mutates or if a virus is outcompeted.
                    if(x.ComputedOnce(TCellCutoff, MutationProbability, RInt, VirusCountLimit, MutationJumpLimit, VirusAmountOnMutate))
                        Graveyard.Add(new Tuple<Individual, string>(x, "T-cell depletion"));
                }
                time++;
                Population = NextPopulation;
            }
            int FINISH = 0;
        }
    }
}
