using System;

namespace Test
{
    class Program
    {
        static void Main(string[] args)
        {
            int D = 2;
            int liczbaKukulek = 7;
            int[] lb = new int[2]{-2,-2};
            int[] ub = new int[2]{2,2};

            
            double[,] nest = new double[liczbaKukulek,D]; 

            double prawdopodobienstwoWykrycia = 0.25;

            int max_iteration = 600;

            for(int i=0;i<liczbaKukulek;i++)
            {
                for(int j=0;j<D;j++)
                {
                    Random rnd = new Random();
                    nest[i, j] = lb[j] + rnd.NextDouble() * (ub[j] - lb[j]);
                    //nest(i,j)=lb(:,j) + rand.*(ub(:,j)-lb(:,j));

                    //Console.Write(nest[i,j] + " |||| ");
                }
                //Console.Write("\n");
            } 
            
            double fx = fns(nest);

            double beta = 3/2;

            double sigma = Math.Pow(MathNet.Numerics.SpecialFunctions.Gamma(1+beta)*Math.Sin(Math.PI*beta/2) / 
                            MathNet.Numerics.SpecialFunctions.Gamma((1+beta)/2)*beta* Math.Pow(2,(beta-1)/2),1/beta);

            //sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta));
            //sigma nie dziala tak jak powinna :c
            //Console.Write(sigma);
            
            for(int iteration = 0;iteration < max_iteration; iteration++)
            {
                //[fnv, indf] = min(fx);
                //best = nest(indf,:);

                for(int j=0;j<liczbaKukulek;j++)
                {
                    //s=nest(j,:);
                    //u=randn(size(s))*sigma;
                    //v=randn(size(s));
                    //step=u./abs(v).^(1/beta);
                    //Xnew = s+rand(size(s)).*0.01.*step.*(s-best);



                }

            }
        }
        static double fns(double[,] array)
        {
            
            MathNet.Numerics.LinearAlgebra.Vector<double> vectX = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(7);
            MathNet.Numerics.LinearAlgebra.Vector<double> vectY = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(7);
            MathNet.Numerics.LinearAlgebra.Vector<double> vectZ = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(7);

            for(int i=0;i<7;i++)
            {
                vectX[i] = array[i , 1];
                vectY[i] = array[i , 2];
            }
            
            //vectZ = -(vectX.Multiply(vectX,vectX) + )

        }
    }
}
