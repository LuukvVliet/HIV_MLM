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
        // Aantal mutaties schaalt met de TOEVOER van virus, niet het totaal aantal virus
    {
        int ID { get; set; }
        int externalTime;
        public StateType InternalState { get; set; }
        public List<double> VirusState { get; set; }
        public List<double> VBetas { get; set; }
        public List<StateType> StateHistory { get; set; }

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
        public bool ComputedOnce(double tcellCutoff, double mr, Random rGen, int VirusLimit, double NextBeta, double newVirusAmount, int VirusGrowth)
        {
            //If a viral strain has more than 10% of the virus amount with which it usually infects, go to the next population.
            double cutoff = 0.99 * newVirusAmount;
            //SumV is used to calculate probability distributions
            double sumV = 0;
            //Checks to see if the individual still has enough T cells to continue living: returns true if not.
            StateType LS = InternalState;
            StateHistory.Add(LS);
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


            int mutations = 0;
            //Mutations according to the total net growth.
            int newSum = (int)sumV;

            //Mutations according to the deterministic average:
            double expectedMut = mr * sumV;
            if (expectedMut > 1)
                mutations = (int)Math.Round(expectedMut);
            else if (rGen.NextDouble() < expectedMut)
                mutations = 1;
            //THIS one line of code is by far the most time consuming:
            //Mutations using a binomial distribution
           /* if (BinomialDatabase != null)
                sample = BinomialDatabase[VirusGrowth / 50].Sample();
            else sample = Binomial.Sample(mr, VirusGrowth);
           */
            this.VirusState = LS.Skip(2).ToList();
            if (VBetas.Count < VirusLimit)
            {
                
                for(int trial = 0; trial < mutations; trial++) { 
                double targetmutation = rGen.Next(0, newSum);
                double target = 0;
                for (int i = 0; i < VBetas.Count; i++)
                {
                    if (VirusState[i] <= targetmutation)
                        targetmutation -= VirusState[i];
                    //On else, this strain is chosen to mutate
                    else
                    {
                        target = VBetas[i]; // Saving the beta of the originally mutated strain.
                        //Remove virus which mutates, up to 'newVirusAmount'
                        if (newVirusAmount - VirusState[i] > 0)
                            newVirusAmount = VirusState[i];
                        VirusState[i] -= newVirusAmount;
                        // Determine whether the next bin or previous bin should be mutated to
                        if (rGen.NextBoolean())
                            target -= NextBeta;
                        else
                            target += NextBeta;
                        if (target < 0) target = 0;
                            // See if this bin already exists in the list
                        bool newV = true;
                        for (int why = 0; why < NewBetas.Count; why++)
                        {
                             if (NewBetas[why] == target)
                                {
                                    NewLS[why + 2] += newVirusAmount;
                                    newV = false;
                                    break;
                                }
                        }

                        // Bin does not exist, add new bin
                        if (newV)
                            {
                                NewLS.Add(newVirusAmount);
                                NewBetas.Add(target);
                            }
                        break;
                    }
                }
               
                
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
