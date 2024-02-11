using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OdeLibrary;

namespace HIV_MLMv1
{
    class Individual
    {
        public LambdaOde VirusDynamics { get; }
        int ID;
        int externalTime;
        StateType currentState;
        List<StateType> History;
        public Individual(int tiem, int id, StateType Init)
        {
            ID = id;
            externalTime = tiem;
            //RungaKutta STEPPER

            VirusDynamics = new LambdaOde
            {
                InitialConditions = Init, // Currently arbitrary initial conditions
                OdeObserver = (x, t) => {History.Add(x); currentState = x; },
                OdeSystem = (x, dxdt, t) => 
                {
                    const double d1 = 0;
                    const double d2 = 0;
                    const double r = 0;
                    const double beta = 0;
                    const double delta = 0;
                    const double a = 0;
                    const double h1 = 0;
                    const double h2 = 0;
                    
                    dxdt[0] = r * x[0] - d1*x[0] - (x[0] *x[1])/h1 - beta * x[0]* x[1]; // T cells
                    dxdt[1] = beta * x[0]*x[1] - delta * x[1] * x[2]; // Infected Cells
                    dxdt[2] = beta * x[0] * x[1] - delta * x[1] * x[2];
                    dxdt[3] = (a*x[1]*x[2])/(h2+x[1]+x[2]) -d2*x[2];  // Effector cells

                    // Vary H2 for individuals?
                    //  
                    // Function dxdt[1] should be a whole ensemble of different functions with different beta. 
                }
            };
        }
        public bool ComputedOnce(double tcellCutoff)
        {
            VirusDynamics.InitialConditions = currentState;
            //Checks to see if the individual still has enough T cells to continue living: returns true if not.
            if (currentState[0] <= tcellCutoff) return true;
            return false;
        }
        public double NewInfection()
        {
           double result = 0;
           return result;
        }
    }
}
