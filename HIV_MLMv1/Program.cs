﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OdeLibrary;
using System.IO;
using System.Globalization;
using MathNet.Numerics.Distributions;

namespace HIV_MLMv1
{
    class Program
    {
        //Simulation parameters, feel free to change
        public static double vir50 = 13938; //Check the amount of target cells actually present in 1 ml of blood
        public static double virMax = 0.317;
        public static double virHill = 1.02;

        public static double DeathProbability = 0;          //Base Deathrate of Individuals
        public static double MutationProbability = 0.00001;    //Mutation chance of virus; should be scaled to amount of virus
        public static double MutationBinDistance = 0.00001;    //The limit of how much a mutation can increase or decrease the beta parameter 
                                                            //NOTE: Increase and decrease of beta is done by a "normal" distribution (or an approximation) so getting the limit is highly unlikely. (as a matter of fact, it is one in ten thousand)
        public static double TCellCutoff = 100;             //At what threshold does an individual not have enough T cells to live anymore?
        public static int N = 100;                          //Starting amount of infected individuals.
        public static double IcellStart = 25;               //Starting amount of infected cells on new mutation. SHOULD BE equal to or higher than VirusAmountOnMutate due to cutoff calculation.
        public static double IcellStartBeta = 0.0001;       //Starting beta of the first infective virus.
        public static int VirusCountLimit = 5;              //How many different strains of virus can exist in one individual?
        public static int VirusAmountOnMutate = 25;         //How many infected cells 'transfer' between two strains when one mutates?

        public static StateType StartingInfected = new StateType { 1000000, 5, 25, 0 }; //State of the first 5 infected individuals
        

        //Recommended to leave these parameters be.
        public static int time = 0;
        public static int timeLimit = 1000000; //int.MaxValue; //timeLimit does (technically) count the timelimit, since the time starts on round 0.
        public static Solver SolveSim = new Solver();
        public static StepperTypeCode UsedStepFunction = StepperTypeCode.RungeKutta4; // ODE step solver function type

        static void Main2(string[] args)
        {
            //Read the outputdirectory
            string DirOut = "";
            DirOut = File.ReadAllText("C:\\Users\\lukxi\\source\\repos\\HIV_MLMv1\\ToBeIgnored\\DirOut.txt");

            int nextID = 0; // ID tracking
            //Initialize the population with a few infected individuals
            List<Individual> Population = new List<Individual> { };
            List<Tuple<Individual, string>> Graveyard = new List<Tuple<Individual, string>> { };
            
            for(int i = 0; i < N; i++)
            {
                List<double> StartVBetas = new List<double>();
                StartVBetas.Add(IcellStartBeta);
                Population.Add(new Individual(time, nextID, StartingInfected, StartVBetas));
                nextID++;
            }
              
            //Create the ODE:
            ODE Sim = new ODE(3);

            // Simulation loop
            Random RInt = new Random();
            SolveSim.StepperCode = UsedStepFunction;
            while (time < timeLimit || !(Population.Count == 0))
            {
                List<double> VirusLoads = new List<double>();
                List<Individual> NextPopulation = new List<Individual> { };
                foreach (Individual x in Population)
                {
                    Sim.SetInit(x.InternalState, x.VBetas);         //Set the internals for the ODE
                    SolveSim.Solve(Sim.VirusDynamics, 0 , 1, 0.01); //Solve the ODE system for one timestep
                    x.InternalState = Sim.VirusDynamics.InitialConditions; //Set the individual to the new internal state.


                    //If something doesnt die, add it to the next population
                    if (RInt.NextDouble() > DeathProbability) NextPopulation.Add(x);
                    else { Graveyard.Add(new Tuple<Individual, string>(x, "Natural death")); continue; }
                    //After something is computed once, see if it mutates or if a virus is outcompeted.
                    if (x.ComputedOnce(TCellCutoff, MutationProbability, RInt, VirusCountLimit, MutationBinDistance, VirusAmountOnMutate, Sim.GetVirusGrowth()))
                    { Graveyard.Add(new Tuple<Individual, string>(x, "T-cell depletion")); continue; }

                    //New infection using distribution samples
                    VirusLoads.Add(x.VirusState.Sum());
                    
                }



                int totalVirusLoad = (int)VirusLoads.Sum();
                double NewInfections = totalVirusLoad;
                if (NewInfections > 1)
                    NewInfections = (int)Math.Round(NewInfections);
                else if (RInt.NextDouble() < NewInfections)
                    NewInfections = 1;

                while (NewInfections > 0)
                {
                    
                    //These are the states of the next infection
                    StateType NextInfection = new StateType { StartingInfected[0], StartingInfected[1], IcellStart };
                    List<double> NextBetas = new List<double>();
                    // By weighted draw, determine which individual and what beta jumps to the next individual.
                    int targetIndv = RInt.Next(0, totalVirusLoad);
                    double TI = (double)targetIndv;
                    int tracker = 0;
                    while (tracker < VirusLoads.Count)
                    {
                        //Weighted drawing of which individual to have infect
                        if(TI < VirusLoads[tracker])
                        {
                            TI -= VirusLoads[tracker];
                        }
                        else
                        {
                            //Weighted drawing of which virus to use for infection
                            int targetVir = RInt.Next(0, (int)VirusLoads[tracker]);
                            double TV = (double)targetVir;
                            int straincounter = 0;
                            foreach(double strain in Population[tracker].VirusState)
                            {
                                if (TV < strain)
                                {
                                    TV -= strain;
                                }
                                else
                                {
                                    NextBetas.Add(Population[tracker].VBetas[straincounter]);
                                }
                                straincounter++;
                            }
                        }
                        tracker++;
                    }

                    Population.Add(new Individual(time, nextID++, NextInfection, NextBetas));
                    NewInfections--;
                };
                time++;
                Population = NextPopulation;
                if(time%1000 == 0 && time != 0)
                {
                    foreach(Tuple<Individual, string> ind in Graveyard)
                    {

                    }
                }
            }
            int FINISH = 0;


            // ADD post-processing of lists
            // Print to .csv type file for further analysis of model function.
        }
    }
}
