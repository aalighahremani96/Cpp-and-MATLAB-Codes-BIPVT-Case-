#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <ctime>
using namespace std;

/*inputs*/
const int loop=1;
const int forced=1;
const double nu=0.000015;
const double ro=1.2;
const double mu=0.0000151;
const double v_amb=2;
const double g=9.8;
const double A_pv=1;
const double b=1;
const double C_p=1005;
const double p_2=1.184;
const double k_c=0.036;
const double k_g=1;
const double k_i=0.035;
const double k_tg=0.033;
const double k_air=0.025;
const double L=1;
const double L_c=pow(300,-6);
const double L_g=0.003;
const double L_i=0.05;
const double L_tg=0.0005;
const double alpha_c=0.8;
const double beta_c=0.9;
const double alpha_tg=0.5;
const double theta=0;
const double etha_th=0.5;
const double t_g=0.95;
const double pr=0.72;
const double D_h=0.09;
const double ductdepth=0.06;
const double A_2=ductdepth*b;
const double I[10]={300, 400, 550, 650, 658, 600, 560, 420, 250, 110};
const double T_r[10]={32, 33, 35, 36, 42, 44, 42, 41, 41, 39,};     /*+273*/
const double T_amb[10]={29, 30, 32, 35, 38, 40, 41, 41, 40, 39};    /*+273*/
const double v[10]={2.67, 1.98, 1.42, 1.87, 1.73, 1.67, 1.77, 1.7, 1.83, 1.77};
const double u_tg=pow(L_tg/k_tg,-1);
const int gg=107000000;
const double h_aircare[10]={2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98, 2.98};
double Re[10]={0};
double gggg=1+pow((0.492/pr),(9/16));
double h_o=1.8+3.8*v_amb;   /*h_o(Ambient Side Natural Convection)*/
double T_sky[10]={0};
double etha_c=0.2;
double alpha;


int main()
{   /*Backsurface, Cell Temperature, Outlet air Temperature Insulation Temperature Initial Guesses*/
    double T_bs[loop][10]={38, 43, 47, 51, 57, 56, 55, 49, 46, 42}; /*+273*/
    double T_c[loop][10]={38, 41, 48, 50, 55, 55, 53, 51, 47, 42}; /*+273*/
    double T_2[loop][10]={33, 36, 39, 42, 46, 46, 46, 45, 43, 40}; /*+273*/
    double T_i[loop][10]={30, 30, 30, 30, 30, 30, 30, 30, 30, 30}; /*+273*/
    double T_airout[loop][10]={{0}};
    double dt[loop][10]={{0}};
    double Ra[loop][10]={{0}};
    double Raa[loop][10]={{0}};
    double Ra_2[loop][10]={{0}};
    double Gr[loop][10]={{0}};
    double f[loop][10]={{0}};
    double h_air[loop][10]={{0}};
    double h_r[loop][10]={{0}};
    double h_rad[loop][10]={{0}};
    double u_b[loop][10]={{0}};
    double u_T[loop][10]={{0}};
    double u_tair[loop][10]={{0}};
    double u_bb[loop][10]={{0}};
    double u_L[loop][10]={{0}};
    double mdot[loop][10]={{0}};
    double xxx[loop][10]={{0}};
    double zzz[loop][10]={{0}};
    double qdot[loop][10]={{0}};
    double electrical1[loop][10]={{0}};
    double thermal[loop][10]={{0}};
    double tinsulation[loop][10]={{0}};
    


    for (int t = 0; t < 10; t++)
    {
        T_sky[t]=T_amb[t]-6;
    }



    for (int w = 0; w < loop; w++)
    {
         for (int t=0; t <10; t++)
        {
            /*Ra, Gr, f*/    
            Ra[w][t]=abs(T_bs[w][t]-T_i[w][t])*gg*pow(ductdepth,3);
            Gr[w][t]=Ra[w][t]/pr;
            f[w][t]= 1.906 * pow((Gr[w][t]/pr),1/12);
            /*h_air(Duct Air)*/
            if (forced==0)
            {   
                dt[w][t]=T_bs[w][t]-T_i[w][t];
                Raa[w][t]=abs(dt[w][t])*(gg*nu*sin(theta))*ro*pow(L,3);
                h_air[w][t]=(k_air/D_h)*0.0965*pow(Raa[w][t],0.29); /*Natural Convection*/
            }
            else
            {   Re[t]=(p_2*D_h/mu)*v[t];

                if (Re[t]>2300)
                {
                    h_air[w][t]=(k_air/D_h)*0.023*pow(Re[t],0.8)*pow(pr,0.4);  /*Forced Convection(Turbulent)*/          
                }
                else
                {
                    h_air[w][t]= (k_air/D_h)*h_aircare[t];  /*Forced Convection(Laminar)*/
                } 
            }
        }

        for (int t=0; t <10; t++)
        {   /*h_r(Room Side Natural Convection)*/
            Ra_2[w][t]=gg*cos(theta)*abs((T_i[w][t]-T_r[t]))*pow(L,3);
            h_r[w][t]=(k_air/L)*(0.68+((0.67*pow(Ra_2[w][t],0.25))/(pow((gggg),(4/9)))));
        }
        for (int t=0; t <10; t++)
        {
        alpha=(t_g*((alpha_c*beta_c)+((alpha_tg)*(1-beta_c))))-(etha_c*beta_c*alpha_c);
        /*Radiative Heat Transfer Coefficient*/
        h_rad[w][t]=0.88*5.6704*pow(10,(-8))*(T_c[w][t]+T_sky[t]+546)*((pow(T_sky[t]+273,2))+(pow(T_c[w][t]+273,2)));
        /*Some More Heat Transfer Coefficients*/
        u_b[w][t]=pow(((L_tg/k_tg)+(1/h_air[w][t])),-1);
        u_T[w][t]=pow(((L_g/k_g) +(1/(h_o))),-1);
        u_tair[w][t]=pow(((1/u_b[w][t])+(1/u_T[w][t])),-1);
        u_bb[w][t]=pow(((1/h_air[w][t])+(1/h_r[w][t])+(L_i/k_i)),-1);
        u_L[w][t]=u_bb[w][t]+u_b[w][t];
        }

        for (int t=0; t <10; t++)
        {
            if (forced==1)
            {   /*Duct Mass Flow Rate is an input for Forced Convection*/
                mdot[w][t]=p_2*A_2*v[t];
            }
            else
            {   /*Duct Mass Flow Rate Should be Calculated for Natural Convection*/
                mdot[w][t]=pow(((2*(pow((A_2*p_2),2))*A_pv*I[t]*etha_th*g*L*sin(theta))/(T_2[w][t]*C_p*(2+((f[w][t]*L)/D_h)))),(1/3));
            }
        

            zzz[w][t]=((u_L[w][t]-((pow(u_b[w][t],2))/(u_T[w][t]+u_b[w][t]+h_rad[w][t])))*(b/(mdot[w][t]*C_p)));
            /*Outlet Air Temperature*/
            T_airout[w][t]=(((((u_bb[w][t])*T_r[t]+((u_b[w][t]*(u_T[w][t]+h_rad[w][t]))/(u_T[w][t]+u_b[w][t]+h_rad[w][t]))*T_amb[t]-6*((h_rad[w][t]*u_b[w][t])/(u_T[w][t]+u_b[w][t]+h_rad[w][t]))+((u_b[w][t])/(u_T[w][t]+u_b[w][t]+h_rad[w][t]))*I[t]*alpha))/((u_L[w][t]-((pow(u_b[w][t],2))/(u_T[w][t]+u_b[w][t]+h_rad[w][t])))))*(1-exp(-zzz[w][t]*L)))+T_r[t]*exp(-zzz[w][t]*L);


        }

    
        for (int t=0; t <10; t++)
        {   /*Backsurface Temperature*/
            T_bs[w][t]=(((h_air[w][t]+((u_tg*u_b[w][t])/(u_T[w][t]+u_b[w][t]+h_rad[w][t])))*((T_airout[w][t]+T_r[t])/2)+((u_tg*u_T[w][t])/(u_T[w][t]+u_b[w][t]+h_rad[w][t]))*T_amb[t]+((u_tg*h_rad[w][t])/(u_T[w][t]+u_b[w][t]+h_rad[w][t]))*T_sky[t]+I[t]*alpha*(u_tg/(u_T[w][t]+u_b[w][t]+h_rad[w][t])))/(u_tg+h_air[w][t]));

        }
        for (int t=0; t <10; t++)
        {   /*Solar Cell Temperature*/
            T_c[w][t]=(u_b[w][t]*((T_airout[w][t]+T_r[t])/2)+u_T[w][t]*T_amb[t]+h_rad[w][t]*T_sky[t]+I[t]*alpha)/(u_T[w][t]+u_b[w][t]+h_rad[w][t]);

        }
        for (int t=0; t <10; t++)
        {   /*heat gained*/
            qdot[w][t]=mdot[w][t]*C_p*(T_airout[w][t]-T_r[t]);
        }
        for (int t=0; t <10; t++)
        {
            electrical1[w][t]=etha_c*(1-(0.006*(T_c[w][t]+273-298)));   /*Electrical Efficiency*/
            thermal[w][t]=(((qdot[w][t]+electrical1[w][t]*b*L*I[t])/(I[t]*b*L)))*100;   /*Total Efficiency*/
        }
        for (int t=0; t <10; t++)
        {   /*Insulation Temperature*/
            tinsulation[w][t]=T_airout[w][t]-(((u_bb[w][t]*L_i)/(k_i))*(T_airout[w][t]-T_r[t]));
        }
        
        for (int j = 0; j < 10; j++)
        {
            T_2[w+1][j]=T_airout[w][j];
        }

    
    }
    /*outputs*/
    for (int i = 0; i < loop; i++)
    {
        cout<<"w is"<<i<<endl;
        for (int j = 0; j < 10; j++)
        {
            cout<<T_airout[i][j]<<endl; /*please change (T_airout[i][j]) with other variables for various output*/
        }
    }

}