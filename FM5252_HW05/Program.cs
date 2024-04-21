using System;
using System.Xml;
using System.Linq;

/*
This program is a simple interface about Monte Carlo simulation.
User will enter the parameters including Stock Price, Strike Price, Risk Free Rate,
Volatility, Time to Maturity (Year), Number of Steps, and Number of Simulations.
Then the interface will show option price, standard error, and greeks for both call and put.

Refferences:
https://stackoverflow.com/questions/5336457/how-to-calculate-a-standard-deviation-array
https://stackoverflow.com/questions/13813166/read-user-input-of-double-type
https://stackoverflow.com/questions/42926301/how-to-return-multiple-values-in-c-sharp-7
https://learn.microsoft.com/en-us/dotnet/standard/base-types/standard-numeric-format-strings
*/

namespace HW5
{
    class program
    {
        static void Main(string[] args)
        {
            Run();
        }

        static void Run()
        {
            Console.WriteLine(" ");
            Console.WriteLine("======================================================================");
            Console.WriteLine("          Monte Carlo Simulator For Option Prices and Greeks          ");
            Console.WriteLine("======================================================================");
            Console.WriteLine(" ");

            // ask the user to enter the parameters
            Console.WriteLine("Please Enter");

            Console.Write("Stock Price: ");
            double S = Convert.ToDouble(Console.ReadLine());
            
            Console.Write("Strike Price: ");
            double K = Convert.ToDouble(Console.ReadLine());

            Console.Write("Risk Free Rate: ");
            double r = Convert.ToDouble(Console.ReadLine());

            Console.Write("Volatility: ");
            double sigma = Convert.ToDouble(Console.ReadLine());

            Console.Write("Time to Maturity (Year): ");
            double T = Convert.ToDouble(Console.ReadLine());

            Console.Write("Number of Steps: ");
            int steps = Convert.ToInt32(Console.ReadLine());

            Console.Write("Number of Simulations: ");
            int N = Convert.ToInt32(Console.ReadLine());
            
            // generate a N x steps random matrix
            NormalGenerator generator = new NormalGenerator();
            double[,] normal = generator.Generate(N, steps);

            MonteCarlo monteCarlo = new MonteCarlo();

            (double[] SimulateCall, double[] SimulatePut) = monteCarlo.Simulation(S, K, r, sigma, T, steps, N, normal);

            // calculate option prices
            (double CallPrice, double PutPrice) = monteCarlo.Price(SimulateCall, SimulatePut, r, T);
            
            // calculate standard error
            (double CallSE, double PutSE) = monteCarlo.StandardError(SimulateCall, SimulatePut);

            Greeks Greeks = new Greeks(S, K, r, sigma, T, steps, N, normal);
            // calculate greeks
            (double CallDelta, double PutDelta) = Greeks.Delta();
            (double CallGamma, double PutGamma) = Greeks.Gamma();
            (double CallVega, double PutVega) = Greeks.Vega();
            (double CallTheta, double PutTheta) = Greeks.Theta();
            (double CallRho, double PutRho) = Greeks.Rho();
            
            // print the results
            Console.WriteLine(" ");
            Console.WriteLine("======================================================================");
            Console.WriteLine("                          Simulation Results                          ");
            Console.WriteLine("======================================================================");
            Console.WriteLine(" ");
            Console.WriteLine("                        Call Option                   Put Option");
            Console.WriteLine("----------------------------------------------------------------------");
            Console.WriteLine("Price:          {0,20:F10}          {1,20:F10}", CallPrice, PutPrice);
            Console.WriteLine("Standard Error: {0,20:F10}          {1,20:F10}", CallSE, PutSE);
            Console.WriteLine("Delta:          {0,20:F10}          {1,20:F10}", CallDelta, PutDelta);
            Console.WriteLine("Gamma:          {0,20:F10}          {1,20:F10}", CallGamma, PutGamma);
            Console.WriteLine("Vega:           {0,20:F10}          {1,20:F10}", CallVega, PutVega);
            Console.WriteLine("Theta:          {0,20:F10}          {1,20:F10}", CallTheta, PutTheta);
            Console.WriteLine("Rho:            {0,20:F10}          {1,20:F10}", CallRho, PutRho);
            Console.WriteLine("----------------------------------------------------------------------");     
        }
    }

    // normally distributed random numbers using the Box-Muller
    class NormalGenerator
    {
        private Random random = new Random();

        public double[,] Generate(int N, int Steps)
        {
            double[,] normals = new double[N, Steps];

            // number of simulation
            for (int i = 0; i < N; i++)
            {   
                // number of steps
                for (int j = 0; j < Steps; j++)
                {
                    double x1 = random.NextDouble();
                    double x2 = random.NextDouble();

                    double z1 = Math.Sqrt(-2 * Math.Log(x1)) * Math.Cos(2 * Math.PI * x2);
                    normals[i, j] = z1;
                }
            }

            return normals;
        }
    }

    // payoff function for call and put
    class Payoff
    {
        public double func(double ST, double K, string OptionType)
        {
            if (OptionType == "call")
            {
                return Math.Max(ST - K, 0);
            }

            else
            {
                return Math.Max(K - ST, 0);
            }
        }
    }

    class MonteCarlo

    {
        public (double[], double[]) 
        Simulation(double S, double K, double r, double sigma, double T, int steps, int N,  double[,] normal)
        {
            Payoff payoff = new Payoff();

            double dt = (double)T / steps;
            
            // payoffs
            double[] CallPayoffs = new double[N];
            double[] PutPayoffs = new double[N];

            // simulation
            for (int i = 0; i < N; i++)
            {
                double ST = S;

                for (int j = 0; j < steps; j++)
                {
                    // price and se
                    ST *= Math.Exp((r - 0.5 * Math.Pow(sigma, 2)) * dt + sigma * Math.Sqrt(dt) * normal[i,j]);
                }

                // apply payoff function
                CallPayoffs[i] = payoff.func(ST, K, "call");
                PutPayoffs[i] = payoff.func(ST, K, "put");
            }

            return (CallPayoffs, PutPayoffs);
        }

        public (double, double) 
        Price(double[] CallPayoffs, double[] PutPayoffs, double r, double T)
        {
            // take the mean and discount back
            double CallPrice = CallPayoffs.Average() * Math.Exp(-r * T);
            double PutPrice = PutPayoffs.Average() * Math.Exp(-r * T);

            return (CallPrice, PutPrice);
        }

        public (double, double) 
        StandardError(double[] CallPayoffs, double[] PutPayoffs)
        {
            // Call Standard Error
            double CallMean = CallPayoffs.Average();
            // calculate std
            double CallSumSquares = CallPayoffs.Select(val => (val - CallMean) * (val - CallMean)).Sum();
            double CallStandardDev = Math.Sqrt(CallSumSquares / (CallPayoffs.Length - 1));
            // calculate standard error
            double CallStandardError = CallStandardDev / Math.Sqrt(CallPayoffs.Length);

            // Put Standard Error
            double PutMean = PutPayoffs.Average();
            double PutSumSquares = PutPayoffs.Select(val => (val - PutMean) * (val - PutMean)).Sum();
            double PutStandardDev = Math.Sqrt(PutSumSquares / (PutPayoffs.Length - 1));
            double PutStandardError = PutStandardDev / Math.Sqrt(PutPayoffs.Length);

            return(CallStandardError, PutStandardError);
        }
    }

    class Greeks
    {
        private double S;
        private double K;
        private double r;
        private double sigma;
        private double T;
        private int steps;
        private int N;
        private double[,] normal;
        MonteCarlo MonteCarlo = new MonteCarlo();

        public Greeks(double S, double K, double r, double sigma, double T, int steps, int N,  double[,] normal)
        {
            this.S = S;
            this.K = K;
            this.r = r;
            this.sigma = sigma;
            this.T = T;
            this.steps = steps;
            this.N = N;
            this.normal = normal;
        }
        
        public (double, double) Delta()
        {
            // ΔS
            double DeltaS = S * 0.01;

            // S + ΔS and S - ΔS for call and put
            (double[] CallPlusDelta, double[] PutPlusDelta) = MonteCarlo.Simulation(S * 1.01, K, r, sigma, T, steps, N, normal);
            (double[] CallMinusDelta, double[] PutMinusDelta) = MonteCarlo.Simulation(S * 0.99, K, r, sigma, T, steps, N, normal);

            (double CallPlusDelta_Price, double PutPlusDelta_Price) = MonteCarlo.Price(CallPlusDelta, PutPlusDelta, r, T);
            (double CallMinusDelta_Price, double PutMinusDelta_Price) = MonteCarlo.Price(CallMinusDelta, PutMinusDelta, r, T);

            // delta call and put
            double CallDelta = (CallPlusDelta_Price - CallMinusDelta_Price) / (2 * DeltaS);
            double PutDelta = (PutPlusDelta_Price - PutMinusDelta_Price) / (2 * DeltaS);

            return(CallDelta, PutDelta);
        }

        public (double, double) Gamma()
        {
            // ΔS
            double DeltaS = S * 0.01;

            // S + ΔS, S, and S - ΔS for call and put
            (double[] Call, double[] Put) = MonteCarlo.Simulation(S, K, r, sigma, T, steps, N, normal);
            (double[] CallPlusDelta, double[] PutPlusDelta) = MonteCarlo.Simulation(S * 1.01, K, r, sigma, T, steps, N, normal);
            (double[] CallMinusDelta, double[] PutMinusDelta) = MonteCarlo.Simulation(S * 0.99, K, r, sigma, T, steps, N, normal);

            (double Call_Price, double Put_Price) = MonteCarlo.Price(Call, Put, r, T);
            (double CallPlusDelta_Price, double PutPlusDelta_Price) = MonteCarlo.Price(CallPlusDelta, PutPlusDelta, r, T);
            (double CallMinusDelta_Price, double PutMinusDelta_Price) = MonteCarlo.Price(CallMinusDelta, PutMinusDelta, r, T);

            // gamma call and put
            double CallGamma = (CallPlusDelta_Price -(2 * Call_Price) + CallMinusDelta_Price) / Math.Pow(DeltaS, 2);
            double PutGamma = (PutPlusDelta_Price -(2 * Put_Price) + PutMinusDelta_Price) / Math.Pow(DeltaS, 2);

            return (CallGamma, PutGamma);
        }

        public (double, double) Vega()
        {
            // ΔSigma
            double DeltaSigma = sigma * 0.01;

            // Sigma + ΔSigma and Sigma - ΔSigma for call and put
            (double[] CallPlusVega, double[] PutPlusVega) = MonteCarlo.Simulation(S, K, r, sigma * 1.01, T, steps, N, normal);
            (double[] CallMinusVega, double[] PutMinusVega) = MonteCarlo.Simulation(S, K, r, sigma * 0.99, T, steps, N, normal);

            (double CallPlusVega_Price, double PutPlusVega_Price) = MonteCarlo.Price(CallPlusVega, PutPlusVega, r, T);
            (double CallMinusVega_Price, double PutMinusVega_Price) = MonteCarlo.Price(CallMinusVega, PutMinusVega, r, T);

            // vega call and put
            double CallVega = (CallPlusVega_Price - CallMinusVega_Price) / (2 * DeltaSigma);
            double PutVega = (PutPlusVega_Price - PutMinusVega_Price) / (2 * DeltaSigma);

            return(CallVega, PutVega);
            
        }

        public (double, double) Theta()
        {
            // ΔT
            double DeltaT = T * 0.01;

             // T + ΔT and T for call and put
            (double[] CallPlusTheta, double[] PutPlusTheta) = MonteCarlo.Simulation(S, K, r, sigma, T * 1.01, steps, N, normal);
            (double[] Call, double[] Put) = MonteCarlo.Simulation(S, K, r, sigma, T, steps, N, normal);

            (double CallPlusTheta_Price, double PutPlusTheta_Price) = MonteCarlo.Price(CallPlusTheta, PutPlusTheta, r, T * 1.01);
            (double Call_Price, double Put_Price) = MonteCarlo.Price(Call, Put, r, T);

            // theta call and put
            double CallTheta = (CallPlusTheta_Price - Call_Price) / DeltaT;
            double PutTheta = (PutPlusTheta_Price - Put_Price) / DeltaT;

            return(CallTheta, PutTheta);

        }
        public (double, double) Rho()
        {
            // Δr
            double DeltaR = r * 0.01;

            // r + Δr and r - Δr for call and put
            (double[] CallPlusRho, double[] PutPlusRho) = MonteCarlo.Simulation(S, K, r * 1.01, sigma, T, steps, N, normal);
            (double[] CallMinusRho, double[] PutMinusRho) = MonteCarlo.Simulation(S, K, r * 0.99, sigma, T, steps, N, normal);

            (double CallPlusRho_Price, double PutPlusRho_Price) = MonteCarlo.Price(CallPlusRho, PutPlusRho, r * 1.01, T);
            (double CallMinusRho_Price, double PutMinusRho_Price) = MonteCarlo.Price(CallMinusRho, PutMinusRho, r * 0.99, T);

            // rho call and put
            double CallRho = (CallPlusRho_Price - CallMinusRho_Price) / (2 * DeltaR);
            double PutRho = (PutPlusRho_Price - PutMinusRho_Price) / (2 * DeltaR);

            return(CallRho, PutRho);
        }

    }

}

