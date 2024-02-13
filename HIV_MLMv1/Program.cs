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
        public static StateType StartingInfected = new StateType{1.0, 0.001, 0.1};
        public static Solver SolveSim = new Solver();
        public static StepperTypeCode UsedStepFunction = StepperTypeCode.RungeKutta4; // ODE step solver function type
        public static int VirusCountLimit = 2; //How many different strains of virus can exist in one individual?

        static void Main(string[] args)
        {
            List<StateType> History = new List<StateType>();
            List<double> xcoord = new List<double>();
            List<double> ycoord = new List<double>();

            LambdaOde VirusDynamics = new LambdaOde
            {
                InitialConditions = {100, 10}, // Currently arbitrary initial conditions
                OdeObserver = (x, t) => { History.Add(x); xcoord.Add(x[0]); ycoord.Add(x[1]);},
                OdeSystem = (x, dxdt, t) =>
                {

                    const double r = 1;
                    const double h = 1000;
                    const double d = 0.6;
                    const double cat = 0.01;

                    dxdt[0] = r * x[0] - r*(Math.Pow(x[0],2)) / h - d * x[0] - cat * x[0] * x[1];
                    dxdt[1] = cat * x[0] * x[1] - d * x[1];

                }
            };
            
            int nextID = 0; // ID tracking

            //Initialize the population with a few infected individuals
            List<Individual> Population = new List<Individual> { };
            
            for(int i = 0; i < N; i++)
            {
                Population.Add(new Individual(time, nextID, StartingInfected));
                nextID++;
            }

            // Simulation loop
            Random RInt = new Random();
            SolveSim.StepperCode = UsedStepFunction;
            while (time < timeLimit)
            {
                List<Individual> NextPopulation = new List<Individual> { };
                foreach (Individual x in Population)
                {
                    SolveSim.Solve(x.VirusDynamics, 0 , 1, 0.01); //Solve the ODE system for one timestep
                    if (RInt.NextDouble() > DeathProbability) NextPopulation.Add(x); 
                    
                    if (RInt.NextDouble() <= InfectProbability) 
                    {
                        RInt.Next();
                        Population.Add(new Individual(time, nextID++, StartingInfected));
                    }; // Make ODE transferrable
                    
                    
                    x.ComputedOnce(TCellCutoff, MutationProbability, RInt, VirusCountLimit);
                }
                time++;
                Population = NextPopulation;
            }
        }
    }
}
