using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OdeLibrary;

namespace HIV_MLMv1
{
    class TestClass
    {

        public LambdaOde VirusDynamics { get; }
        int ID;
        int externalTime;
        StateType currentState;
        List<StateType> History;

        public TestClass(int tiem, int id, StateType Init)
        {
            ID = id;
            externalTime = tiem;
            //RungaKutta STEPPER

            VirusDynamics = new LambdaOde
            {
                InitialConditions = Init, // Currently arbitrary initial conditions
                OdeObserver = (x, t) => { History.Add(x);},
                OdeSystem = (x, dxdt, t) =>
                {
                    
                    const double r = 1;
                    const double h = 1000;
                    const double d = 0.6;
                    const double cat = 0.01;

                    dxdt[0] = r*x[0] - x[0]/h - d*x[0] -cat*x[0]*x[1];
                    dxdt[1] = cat*x[0]*x[1] - d*x[1];

                }
            };
        }
    }
}
