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
        //Simulation parameters, with N being the starting size of the population.
        public static double InfectProbability = 0.15; //Infection of new individuals
        public static double DeathProbability = 0.02;  //Base Deathrate of Individuals
        public static double MutationProbability = 0.05; //
        public static double TCellCutoff = 0.05; // At what threshold does an individual not have enough T cells to live anymore?
        public static int N = 5;
        public static int time = 0;
        public static int timeLimit = int.MaxValue; //timeLimit does (technically) count the timelimit, since the time starts on round 0.
        public static StateType StartingInfected = new StateType{10000, 10, 100};
        public static double IcellStart = 100;
        public static Solver SolveSim = new Solver();
        public static StepperTypeCode UsedStepFunction = StepperTypeCode.RungeKutta4; // ODE step solver function type
        public static int VirusCountLimit = 2; //How many different strains of virus can exist in one individual?

        static void Main(string[] args)
        {
            int nextID = 0; // ID tracking
            int xxx = StartingInfected.Count;
            //Initialize the population with a few infected individuals
            List<Individual> Population = new List<Individual> { };
            List<Individual> Graveyard = new List<Individual> { };
            for(int i = 0; i < N; i++)
            {
                Population.Add(new Individual(time, nextID, StartingInfected));
                nextID++;
            }

            // Simulation loop
            Random RInt = new Random();
            SolveSim.StepperCode = UsedStepFunction;
            while (time < timeLimit || Population.Count == 0)
            {
                List<Individual> NextPopulation = new List<Individual> { };
                foreach (Individual x in Population)
                {
                    SolveSim.Solve(x.VirusDynamics, 0 , 1, 0.01); //Solve the ODE system for one timestep
                    //If something doesnt die, add it to the next population
                    if (RInt.NextDouble() > DeathProbability) NextPopulation.Add(x);
                    else Graveyard.Add(x);
                    
                    //New infection
                    if (RInt.NextDouble() <= InfectProbability) 
                    {
                        StateType NextInfection = new StateType { StartingInfected[0], StartingInfected[1], IcellStart };
                        NextPopulation.Add(new Individual(time, nextID++, NextInfection, new List<Tuple<double, double>> { new Tuple<double, double>(IcellStart, x.NewInfection(RInt)) }));
                    }; 
                    
                    //After something is computed once, see if it mutates or if a virus is outcompeted.
                    x.ComputedOnce(TCellCutoff, MutationProbability, RInt, VirusCountLimit);
                }
                time++;
                Population = NextPopulation;
            }
            int FINISH = 0;
        }
    }
}
