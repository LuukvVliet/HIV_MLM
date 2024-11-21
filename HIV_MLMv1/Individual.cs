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
        
        public int ID { get; set; }
        public int externalTime { get; set; }
        public int DeathTimer { get; set; }
        public List<double> MutDistribution { get; set; }
        public StateType InternalState { get; set; }
        public List<double> VirusState { get; set; }
        public List<int> IntBetas { get; set; }
        public List<List<double>> StateHistory { get; set; }
        public List<List<int>> BetasHistory { get; set; }

        public List<double> expectedMutations { get; set; }

        public Individual(int tiem, int id, StateType Init, List<int> VB) 
        {
            ID = id;
            externalTime = tiem;
            InternalState = Init;
            IntBetas = VB;
            VirusState = Init.Skip(2).ToList();
            StateHistory = new List<List<double>>();
            BetasHistory = new List<List<int>>();
            MutDistribution = new List<double> {0.8, 0.15, 0.024, 0.02, 0.004, 0.001, 0.001};
            expectedMutations = new List<double> { 0.0, 0.0 };
        }
        

        public bool ComputedOnce(double tcellCutoff, double mr, Random rGen, double newVirusAmount, int VirusGrowth)
        {
            return ComputedOnce(tcellCutoff, mr, rGen, 99999, newVirusAmount, VirusGrowth);
        }
        public bool ComputedOnce(double tcellCutoff, double mr, Random rGen, int VirusLimit, double newVirusAmount, int VirusGrowth)
        {
            //If a viral strain has more than 10% of the virus amount with which it usually infects, go to the next population.
            double cutoff = 0.9 * newVirusAmount;
            //SumV is used to calculate probability distributions
            double sumV = 0;
            //Checks to see if the individual still has enough T cells to continue living: returns true if not.
            StateType LS = InternalState;

            // This piece of code keeps history of what happens in the individual; causes OutOfMemory errors by causing insane lengths of lists.
            StateHistory.Add(LS.ToList());
            BetasHistory.Add(IntBetas);
            /* */
            if(LS.Count <3)
            {
                return true;
            }
            if (LS[0] <= tcellCutoff)
            {
                DeathTimer++;
                if(DeathTimer > 20)
                    return true;
            }
            else
                DeathTimer = 0;
            
            

            //Updates the VirusDistribution
            StateType NewLS = new StateType();
            List<int> NewBetas = new List<int>();
            NewLS.Add(LS[0]); NewLS.Add(LS[1]);

            int x = 2;
            foreach(int virusB in IntBetas)
            {
                if ((LS.Count - 2) == 2 * IntBetas.Count) // Check whether the compartimental model is being ran
                    if (LS[x] + LS[x + 1] > cutoff)
                    {

                        NewLS.Add(LS[x]);
                        NewLS.Add(LS[x + 1]);
                        NewBetas.Add(virusB);
                        sumV += LS[x];
                        x += 2;
                    }
                    else
                    {
                        x += 2;
                    }
                else
                {
                    if (LS[x] > cutoff)
                    {

                        NewLS.Add(LS[x]);
                        NewBetas.Add(virusB);
                        sumV += LS[x];
                        x++;
                    }
                    else
                    {
                        x += 2;
                    }
                }
            }


            int mutations = 0;
            //Mutations according to the total net growth.
            int newSum = (int)sumV;

            //Mutations according to the deterministic average:
            double expectedMut = mr * sumV;
            expectedMutations[0] = expectedMut;
            expectedMutations[1] += expectedMut;
            if (expectedMut > 1)
                mutations = (int)Math.Round(expectedMut);
            else if (rGen.NextDouble() < expectedMut)
                mutations = 1;

            mutations = 0;
            this.VirusState = LS.Skip(2).ToList();
            if (IntBetas.Count < VirusLimit)
            {
                
                for(int trial = 0; trial < mutations; trial++) { 
                double targetmutation = rGen.Next(0, newSum);
                int target;
                //This code works for the latent cell population model.
                for (int i = 0; (i/2) < IntBetas.Count; i+=2)
                {
                    if (VirusState[i] <= targetmutation)
                        targetmutation -= VirusState[i];
                    //On else, this strain is chosen to mutate
                    else
                    {
                        //Remove virus which mutates, up to 'newVirusAmount'
                        if (newVirusAmount - VirusState[i] > 0)
                            newVirusAmount = VirusState[i];
                        VirusState[i] -= newVirusAmount;


                       target = IntBetas[i / 2]; // Saving the beta of the originally mutated strain.
                            double whatjump = rGen.NextDouble();
                            if (rGen.NextBoolean())  // Determine whether the next bin or previous bin should be mutated to
                            {
                                foreach(double pos in MutDistribution)
                                {
                                    if (pos == MutDistribution[MutDistribution.Count - 1])
                                    { target = rGen.Next(0, 35); break; }

                                    target++;
                                    if (pos > whatjump)
                                        break;
                                    whatjump -= pos;
                                }
                            }
                            else
                            {
                                foreach (double pos in MutDistribution)
                                {
                                    if (pos == MutDistribution[MutDistribution.Count - 1])
                                    { target = rGen.Next(0, 35); break; }
                                    target--;
                                    if (pos > whatjump)
                                        break;
                                    whatjump -= pos;
                                }
                            }
                        if (target < 0) target = 0;
                            // See if this bin already exists in the list


                        bool newV = true;
                        for (int why = 0; why < NewBetas.Count; why++)
                        {
                             if (NewBetas[why] == target)
                                {
                                    NewLS[why*2 + 2] += newVirusAmount;
                                    newV = false;
                                    break;
                                }
                        }

                        // Bin does not exist, add new bin
                        if (newV)
                            {
                                NewLS.Add(newVirusAmount);
                                NewLS.Add(0);
                                NewBetas.Add(target);
                            }
                        break;
                    }
                }
               
                
                }
            }
            //This resolves memory issues (somewhat)
            NewLS.Capacity = NewLS.Count;
            NewBetas.Capacity = NewBetas.Count;

            InternalState = NewLS;
            IntBetas = NewBetas;
            VirusState = NewLS.Skip(2).ToList();

            

            return false;
        }
    }
}
