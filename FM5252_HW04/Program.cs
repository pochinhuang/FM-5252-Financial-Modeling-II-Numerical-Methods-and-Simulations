using System;

/*
This file contains three methods to generate Gaussian random numbers.
The user will select one method from Sum Twelve, Box-Muller, and Polar Rejection.
Then enter a ρ to generate correlated Gaussian random numbers.
*/

/*
References:
https://blog.csdn.net/weixin_34288121/article/details/92505865
https://learn.microsoft.com/zh-tw/dotnet/csharp/language-reference/keywords/access-modifiers
https://stackoverflow.com/questions/74426378/get-set-vs-get-private-set
https://www.w3schools.com/cs/cs_for_loop.php
*/


namespace Methods
{
    class program
    {
        static void Main(string[] args)
        {
            method_selection();
        }

        static void method_selection()
        {   
            // ask the user to enter a method
            Console.WriteLine("Please Enter a Method: \nsum_twelve, box_muller, polar_rejection");
            string Enter = Console.ReadLine();

            switch (Enter)
            {
                case "sum_twelve":
                   sum_twelve sum = new sum_twelve();
                   Console.WriteLine("Sum Twelve Method: {0}", sum.result_A);

                    double p_sum = correlation();
                    sum.sum_corr(p_sum);
                    Console.WriteLine("Sum Twelve Method Correlation: {0} {1}", sum.z_1, sum.z_2);                  

                   break;

                case "box_muller":
                    box_muller box = new box_muller();
                    Console.WriteLine("Box Muller Method: {0} {1}", box.z1, box.z2);

                    double p_box = correlation();
                    box.box_corr(p_box);
                    Console.WriteLine("Box Muller Method Correlation: {0} {1}", box.z_1, box.z_2);

                    break;
               
                case "polar_rejection":
                    polar_rejection polar = new polar_rejection();
                    Console.WriteLine("Polar Rejection Method: {0} {1}", polar.z_one, polar.z_two);

                    double p_pol = correlation();
                    polar.polar_corr(p_pol);
                    Console.WriteLine("Polar Rejection Method Correlation: {0} {1}", polar.z_1, polar.z_2);

                    break;
                
                // if the input is incorrect, re-ask the user to enter a method
                default:
                    Console.WriteLine("\n" + 
                                    "Error:\n" +
                                    "Incorrect Method");
                    method_selection();
                    
                    break;
            }
        }

        // ask the user to enter a ρ between -1 to 1
        // if the input is outside between -1 to 1, re-ask the user to enter a ρ between -1 to 1
        static double correlation()
        {
            double p;
            do
            {
                Console.WriteLine("Enter a correlation value between -1 to 1 :");
                p = double.Parse(Console.ReadLine());

                if (p < -1 || p > 1)
                {
                    Console.WriteLine("\n" +
                                      "Error:\n" +
                                      "Invalid correlation value.");
                }

            } while (p < -1 || p > 1);

            return p;
        }
    }
    class sum_twelve
    {
        private Random random = new Random();
        public double[] numbers = new double[12];
        public double result_A { get; private set; }
        public double result_B { get; private set; }
        public double z_1 { get; private set; }
        public double z_2 { get; private set; }
    public sum_twelve()
    {   
        // generate twelve uniform random numbers on the unit interval
        for (int i = 0; i < 12; i++)
        {
            numbers[i] = random.NextDouble();
        }
        double sum = 0;
        // add these twelve values together
        foreach (double num in numbers)
        {
            sum += num;
        }
        // subtract six from the total
        result_A = sum - 6;

        // result B is for sum_corr epsilon 2
        for (int i = 0; i < 12; i++)
        {
            numbers[i] = random.NextDouble();
        }
        sum = 0;
        foreach (double num in numbers)
        {
            sum += num;
        }
        result_B = sum - 6;
        }

        //In addition to Gaussian random numbers, we’d like to be able to generate correlated Gaussian random numbers
        public void sum_corr(double p)
        {
            z_1 = result_A;
            z_2 = (p * result_A) + Math.Sqrt(1 - Math.Pow(p, 2)) * result_B;
        }
            
        
    }

    class box_muller
    {
        private Random random = new Random();

        public double x1;
        public double x2;
        public double z1 { get; private set; }
        public double z2 { get; private set; }

        public double z_1 { get; private set; }
        public double z_2 { get; private set; }

        public box_muller()
        {
            // Generate two uniform random values on the unit interval
            x1 = random.NextDouble();
            x2 = random.NextDouble();

            // Input them to the following expressions
            // The resulting output is a pair of Gaussian random variables
            z1 = Math.Sqrt(-2 * Math.Log(x1)) * Math.Cos(2 * Math.PI * x2);
            z2 = Math.Sqrt(-2 * Math.Log(x1)) * Math.Sin(2 * Math.PI * x2);

        }
        //In addition to Gaussian random numbers, we’d like to be able to generate correlated Gaussian random numbers        
        public void box_corr(double p)
        {
            z_1 = z1;
            z_2 = (p * z1) + Math.Sqrt(1 - Math.Pow(p, 2)) * z2;
        }
            
    }


    class polar_rejection
    {
        private Random random = new Random();
        private double x_one;
        private double x_two;
        private double w;
        public double z_one { get; private set; }
        public double z_two { get; private set; }
        public double z_1 { get; private set; }
        public double z_2 { get; private set; }

        public polar_rejection()
        {   
            //If w > 1, repeat 
            do
            {
                // Generate two uniform random values, x1, x2
                x_one = random.NextDouble();
                x_two = random.NextDouble();
                // w = (x_1)^2 + (x_1)^2
                w = Math.Pow(x_one, 2) + Math.Pow(x_two, 2); 
            } while (w > 1);
            
            double c = Math.Sqrt(-2 * Math.Log(w) / w);
            z_one = c * x_one;
            z_two = c * x_two;
        }

        //In addition to Gaussian random numbers, we’d like to be able to generate correlated Gaussian random numbers
        public void polar_corr(double p)
        {
            z_1 = z_one;
            z_2 = (p * z_one) + Math.Sqrt(1 - Math.Pow(p, 2)) * z_two;
        }
    }

}


