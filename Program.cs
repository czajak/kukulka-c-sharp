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
            Random rnd = new Random();
            
            double[,] nest = new double[liczbaKukulek,D]; 

            double prawdopodobienstwoWykrycia = 0.25;

            int max_iteration = 600;

            for(int i=0;i<liczbaKukulek;i++)
            {
                for(int j=0;j<D;j++)
                {
                    nest[i, j] = lb[j] + rnd.NextDouble() * (ub[j] - lb[j]);
                    //nest(i,j)=lb(:,j) + rand.*(ub(:,j)-lb(:,j));

                    //Console.Write(nest[i,j] + " |||| ");
                }
                //Console.Write("\n");
            } 
            
            double fx = fns(nest);

            double beta = 3/2;

            double sigma = Math.Pow(MathNet.Numerics.SpecialFunctions.Gamma(1+beta)*Math.Sin(Math.PI*beta/2) / (MathNet.Numerics.SpecialFunctions.Gamma((1+beta)/2)*beta* Math.Pow(2,(beta-1)/2)),1/beta);

            //Console.Write(sigma);
            
            double[] s = new double[D];
            double[] u = new double[D];
            double[] v = new double[D];
            double[] step = new double[D];
            double[] Xnew = new double[D];
            double[] fnew = new double[liczbaKukulek];

            for(int iteration = 0;iteration < max_iteration; iteration++)
            {
                double[] minVect = new double[2];
                //minVect =

                //WYMAGA ZNALEZIENIA ODPOWIEDNIKA FUNKCJI min() MATLABA
                //[fnv, indf] = min(fx);

                double[] best = new double[D];
                for(int i=0;i<D;i++)
                {
                    best[i] = nest[indf,i];

                    //KOD NIE KOMPILUJE SIĘ PONIEWAŻ NIE MA ZMIENNEJ indf

                }
                for(int j=0;j<liczbaKukulek;j++)
                {
                    

                    for(int i=0;i<D;i++)
                    {
                        s[i] = nest[j,i];
                        u[i] = rnd.NextDouble()*sigma;
                        v[i] = rnd.NextDouble();
                        step[i] = u[i] / Math.Pow(Math.Abs(v[i]),(1/beta));
                        Xnew[i] = s[i] + (rnd.NextDouble() * 0.01 * step[i] * (s[i] - best[i]));

                        if(Xnew[i]>ub[i])
                        {
                            Xnew[i] = ub[i];
                        }  
                        else if(Xnew[i]<lb[i])
                        {
                            Xnew[i] = lb[i];
                        }
                    }
                    fnew = fns(Xnew);
                    /*      
                    
                    if fnew>fx(j,:)
                        nest(j,:)=Xnew;
                        fx(j,:)=fnew;
                        end         
                    end
                    [fmin, K1] = min(fx); 
                    best = nest(K1,:);    
                    */




                    /*
                            
                            
                            
                            K=rand(size(nest))<pa;
                            stepsizeK=rand*(nest(randperm(N),:)-nest(randperm(N),:));
                            new_nest=nest+stepsizeK.*K;
                            
                            for ii=1:size(nest,1)
                                s = new_nest(ii,:);
                                for j=1:size(s,2)
                                    if s(i)>ub(i)
                                        s(i)=ub(i);
                                    elseif s(i)<lb(i)
                                        s(i)=lb(i);
                                    end
                                end
                                new_nest(ii,:) = s;
                                
                                fnew = fns(s);
                                if fnew>fx(ii,:)
                                    nest(ii,:)=s;
                                    fx(ii,:)=fnew;
                                end
                            end
                            
                            [optival, optind] = min(fx);
                            Bestfx(iter) = optival;
                            BestX(iter,:)=nest(optind,:);
                            
                            disp(['Iteracja = ' num2str(iter) ' Najlepszy wynik = ' num2str(Bestfx(iter))]);
                        end
                        best
                        plot(Bestfx, 'linewidth', 2);
*/
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
