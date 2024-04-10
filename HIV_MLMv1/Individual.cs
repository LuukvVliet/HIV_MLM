using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OdeLibrary;
using MathNet.Numerics.Random;
using MathNet.Numerics.Distributions;

namespace HIV_MLMv1
{
    class Individual

        // POSSIBLE OPTIMIZATIONS:
        // VECTOR MULTIPLICATION INSTEAD OF LOOPING OVER LIST
        // PREMPTIVELY CALCULATING BINOMIAL DISTRIBUTIONS
        // Aantal mutaties schaalt met de TOEVOER van virus, niet het totaal aantal virus
    {
        public LambdaOde VirusDynamics { get; set; }
        
        int ID { get; set; }
        int externalTime;
        public StateType InternalState { get; set; }
        public List<double> VirusState { get; set; }
        public List<double> VBetas { get; set; }

        public Individual(int tiem, int id, StateType Init, List<double> VB) 
        {
            ID = id;
            externalTime = tiem;
            InternalState = Init;
            VBetas = VB;
            VirusState = Init.Skip(2).ToList();
        }
        
        public bool ComputedOnce(double tcellCutoff, double mr, Random rGen, double jumplimit, double newVirusAmount, int VirusGrowth)
        {
            return ComputedOnce(tcellCutoff, mr, rGen, 99999, jumplimit, newVirusAmount, VirusGrowth);
        }
        public bool ComputedOnce(double tcellCutoff, double mr, Random rGen, int VirusLimit, double jumplimit, double newVirusAmount, int VirusGrowth)
        {
            //If a viral strain has more than 10% of the virus amount with which it usually infects, go to the next population.
            double cutoff = 0.1 * newVirusAmount;
            //SumV is used to calculate probability distributions
            double sumV = 0;
            //Checks to see if the individual still has enough T cells to continue living: returns true if not.
            StateType LS = InternalState;
            if (LS[0] <= tcellCutoff)
                return true;
            
            //Updates the VirusDistribution
            List<Tuple<double, double>> tempVirusDistribution = new List<Tuple<double, double>>();
            StateType NewLS = new StateType();
            List<double> NewBetas = new List<double>();
            NewLS.Add(LS[0]); NewLS.Add(LS[1]);

            int x = 2;
            foreach(double virusB in VBetas)
            {
                if (LS[x] > cutoff)
                {
                    
                    NewLS.Add(LS[x]);
                    NewBetas.Add(virusB);
                    sumV += LS[x];
                    x++;
                }
                else
                    Console.WriteLine(virusB);
            }


            //Mutations using a binomial distribution
            //Mutations according to the total net growth.
            int newSum = (int)sumV;
            var sample = Binomial.Sample(mr, VirusGrowth);

            this.VirusState = LS.Skip(2).ToList();
            if (VBetas.Count < VirusLimit)
            {
                
                for(int trial = 0; trial < sample; trial++) { 
                double targetmutation = rGen.Next(0, newSum);
                double target = 0;
                for (int i = 0; i < VBetas.Count; i++)
                {
                    if (VirusState[i] <= targetmutation)
                        targetmutation -= VirusState[i];
                    else
                    {
                        target = VBetas[i];
                        //Remove virus which mutates, up to 'newVirusAmount'
                        if (newVirusAmount - VirusState[i] > 0)
                            newVirusAmount = VirusState[i];
                        //Adds the new virus
                        NewLS.Add(VirusState[i] - newVirusAmount);
                        break;
                    }
                }
                target += ((double)(rGen.Next(1, 100) - rGen.Next(1, 100))) / 100 * jumplimit;
                if (target < 0) target = 0; //dont have a negative beta
                //Adds the new virus' new beta
                NewBetas.Add(target);
                }
            }

            InternalState = NewLS;
            VBetas = NewBetas;

            return false;
        }
        public double NewInfection(Random rGen)
        {
            double result = 0;
            double sum = InternalState.Skip(2).Sum();
            double targetmutation = rGen.Next(0, (int)sum);
            for (int i = 0; i < VirusState.Count; i++)
            {
                if (VirusState[i] <= targetmutation)
                    targetmutation -= VirusState[i];
                else
                {
                    result = VBetas[i];
                    break;
                }
            }
            return result;
        }
    }
}
