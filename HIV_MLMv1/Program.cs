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
        public static double InfectProbability = 0.15;
        public static double DeathProbability = 0.02;
        public static double TCellCutoff = 0.05; // At what threshold does an individual not have enough T cells to live anymore?
        public static int N = 5;
        public static int time = 0;
        public static int timeLimit = int.MaxValue; //timeLimit does (technically) count the timelimit, since the time starts on round 0.
        public static StateType StartingInfected = new StateType{1.0, 0.001, 0.1};
        static void Main(string[] args)
        {
            List<StateType> History = new List<StateType>();
            Solver SolveTest = new Solver();
            LambdaOde VirusDynamics = new LambdaOde
            {
                InitialConditions = {100, 10}, // Currently arbitrary initial conditions
                OdeObserver = (x, t) => { History.Add(x); },
                OdeSystem = (x, dxdt, t) =>
                {

                    const double r = 1;
                    const double h = 1000;
                    const double d = 0.6;
                    const double cat = 0.01;

                    dxdt[0] = r * x[0] - x[0] / h - d * x[0] - cat * x[0] * x[1];
                    dxdt[1] = cat * x[0] * x[1] - d * x[1];

                }
            };
            SolveTest.ConvenienceSolve(VirusDynamics, 0, 0.001, 100);
            int i = 0;
            /*
            int nextID = 0; // ID tracking

            //Initialize the population with a few infected individuals
            List<Individual> Population = new List<Individual> { };
            
            for(int i = 0; i < N; i++)
            {
                Population.Add(new Individual(time, nextID, StartingInfected));
                nextID++;
            }


            Random RInt = new Random();
            // Simulation loop
            
            Solver SolveSim = new Solver();
            while (time < timeLimit)
            {
                List<Individual> NextPopulation = new List<Individual> { };
                foreach (Individual x in Population)
                {
                    // Solver.ConvenienceSolve(System, from timepoint, to timepoint, in steps)
                    SolveSim.Solve(x.VirusDynamics, 0 , 1, 0.01);
                    if (RInt.NextDouble() <= InfectProbability) 
                    {
                        
                        RInt.Next();
                        Population.Add(new Individual(time, nextID++, StartingInfected));
                    }; // Make ODE transferrable
                    
                    if(RInt.NextDouble() > DeathProbability) NextPopulation.Add(x); 
                    x.ComputedOnce(TCellCutoff);
                }
                time++;
                Population = NextPopulation;
            }*/
        }
    }
}
