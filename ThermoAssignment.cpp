#include<iostream>
#include<cmath>
using namespace std;
int main()
{

    float P1,T1,P2,T2,Tc,Pc,w;
    cout<<"Enter initial pressure :";
    cin>>P1;
    cout<<"Enter initial temperature :";
    cin>>T1;
    cout<<"Enter final pressure :";
    cin>>P2;
    cout<<"Enter final temperature :";
    cin>>T2;
    cout<<"Enter critical pressure :";
    cin>>Pc;
    cout<<"Enter critical temperature :";
    cin>>Tc;
    cout<<"Enter acentric factor :";
    cin>>w;
    float Pr1,Tr1,B0,B1,P1sat,Pr1_sat;
    float A,B,C,D;
    cout<<"Now, enter constants of Antoine Equation to calculate Psat (Kelvin and Pascal)"<<endl;
    cout<<"Enter A:";
    cin>>A;
    cout<<"Enter B:";
    cin>>B;
    cout<<"Enter C:";
    cin>>C;
    cout<<endl;
    float a,b,c,d,T;
    cout<<"Enter A of Cp:";
    cin>>a;
    cout<<"Enter B of Cp:";
    cin>>b;
    cout<<"Enter C of Cp:";
    cin>>c;
    cout<<"Enter D of Cp:";
    cin>>d;
    Pr1=P1/Pc;
    Tr1=T1/Tc;
    B0=0.083-(0.422/pow(Tr1,1.6));
    B1=0.139-(0.172/pow(Tr1,4.2));
    P1sat = exp(A-(B/(T1+C)));
    Pr1_sat= P1sat/Pc;
    bool check1=0;
    if(Tr1>=3) check1=1;
    else if(Tr1<=0.7)
    {
        if(Pr1<=Pr1_sat) check1=1;
        else
        {
            cout<<"Conditions of validity not satisfied at inlet since P1>P1sat!"<<endl;
        }
    }
    else
    {
        float Pr1UL,Pr1_satc,P1satc;
        P1satc = exp(A-(B/((Tc*0.7)+C)));
        Pr1_satc= P1satc/Pc;
        Pr1UL = (((Tr1-0.7)/2.3)*(3-Pr1_satc))+ Pr1_satc;
        if(Pr1<=Pr1UL) check1=1;
        else
        {
            cout<<"Conditions of validity not satisfied at inlet since P1>P1UL!"<<endl;
        }
    }

    cout<<endl;
    float Z1,Z0,Z;
    Z1=B1*Pr1/Tr1;
    Z0=1+(B0*Pr1/Tr1);
    Z=Z0+(w*Z1);
    float HR1,SR1;
    float R=8.314;
    float HR2,SR2;

    if(check1==1)
    {
        HR1=Pr1*R*Tc*(B0-(0.675/pow(Tr1,1.6))+ w*(B1-(0.722/pow(Tr1,4.2))));
        SR1=-1*Pr1*R*((0.675/pow(Tr1,2.6))+ (w*0.722/pow(Tr1,5.2)));
        cout<<"Residual Enthalpy for first real gas calculation is (+HR1): "<<(HR1)<<endl;
        cout<<"Residual Entropy for first real gas calculation is (+SR1) : "<<(SR1)<<endl;
    }
    cout<<endl;

    float Hig,Cp,Sig;
    T=T2/T1;
    Sig = R*((a*log(T))+ ((T-1)*((b*T1)+(((c*T1*T1)+(d/(T*T*T1*T1)))*((T+1)/2)))))-(R*log(P2/P1));
    Hig = R*(T2-T1)*(a+(b*T1*(T+1)/2)+(c*T1*T1*(pow(T,2)+T+1)/3)+(d/(T*T1*T1)));
    cout<<"Enthalpy change during ideal gas transition is (dHig):"<<Hig<<endl;
    cout<<"Entropy change during ideal gas transition is (dSig):"<<Sig<<endl;
    cout<<endl;

    float Pr2,Tr2,Pr2_sat,P2sat;
    Pr2=P2/Pc;
    Tr2=T2/Tc;
    B0=0.083-(0.422/pow(Tr2,1.6));
    B1=0.139-(0.172/pow(Tr2,4.2));
    P2sat = exp(A-(B/(T2+C)));
    Pr2_sat= P2sat/Pc;
    float check2=0;
    if(Tr2>=3) check2=1;
    else if(Tr2<=0.7)
    {
        if(Pr2<=Pr2_sat) check2=1;
        else
        {
            cout<<"Conditions of validity not satisfied at outlet since P2>P2sat!"<<endl;

        }
    }
    else
    {
        float Pr2UL,Pr2_satc,P2satc;
        P2satc = exp(A-(B/((Tc*0.7)+C)));
        Pr2_satc= P2satc/Pc;
        Pr2UL = (((Tr2-0.7)/2.3)*(3-Pr2_satc))+ Pr2_satc;
        if(Pr2<=Pr2UL) check2=1;
        else
        {
            cout<<"Conditions of validity not satisfied at outlet since P2>P2UL!"<<endl;
        }
    }

    Z1=B1*Pr2/Tr2;
    Z0=1+(B0*Pr2/Tr2);
    Z=Z0+(w*Z1);

    if(check2==1)
    {

        HR2=Pr2*R*Tc*(B0-(0.675/pow(Tr2,1.6))+ w*(B1-(0.722/pow(Tr2,4.2))));
        SR2=-1*Pr2*R*((0.675/pow(Tr2,2.6))+ (w*0.722/pow(Tr2,5.2)));
        cout<<"Residual Enthalpy for second real gas calculation is (+HR2) : "<<HR2<<endl;
        cout<<"Residual Entropy for second real gas calculation is (+SR2): "<<SR2<<endl;

    }

    float H,S;

    if(check1==1&&check2==1)
    {
        H=Hig+HR2-HR1;
        S=Sig+SR2-SR1;
        cout<<"Total enthalpy change is (DH): "<<H<<endl;
        cout<<"Total entropy change is (DS): "<<S<<endl;
    }
    else
    {
        cout<<endl;
        cout<<"FINAL RESULT: Total entropy and enthalpy change cannot be calculated due to lack of full information"<<endl;
    }

}


