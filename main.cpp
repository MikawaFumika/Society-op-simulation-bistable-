//
//  main.cpp
//  Bistable simulation
//
//  Created by Mikawa fumika on 15/5/22.
//  Copyright (c) 2015年 Ye yusong. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <time.h>


int main(int argc, const char * argv[])
{
    srand(time(0));
    
    int cell[50];
    
    double po[50];
    
    
    for(int i=0;i<25;i++)
        po[i]=(8.0f/8.0f);
    
    for(int i=25;i<50;i++)
        po[i]=(8.0f/8.0f);

    
    using namespace std;

            for(int i=0;i<50;i++)
            {double rdb=(double)(rand()%10000)/10000;
                if (rdb<0.5)
                    cell[i]=-1;
                else
                    cell[i]=1;}
 
    
    int loopCount=10000000;
    double T=0.2;
    double alpha=(double)1.0-(T/(((double)1.0)+T));
    double r1=1.1;
    double r2=1.1;
    
    ofstream file;
    file.open("result.txt");
    file<<"TotalLoopCount="<<loopCount<<endl;
    file<<"Alpha="<<alpha<<endl;
    file<<"T="<<T<<endl;
    file<<"R"<<"\t"<<r1<<"\t"<<r2<<endl;
    
    //设数组统计量
    double m1,m2;
   
    
    double M[51];
    for(int i=0;i<51;i++)
        M[i]=0;
   
    
    
    int state1=0;
    int state2=0;
    int transition=1;
    
    long int statetimesmall=0;
    long int statetimelarge=0;
    
    //主程序
    for(int i=0;i<loopCount;i++)
    {
        double rd=(double)(rand()%1000000)/1000000;
        
        int pos1=(rand()%50);
        
        m1=0;
        
        for(int j=0;j<50;j++)
            if(cell[j]==1)
                m1=m1+1;
        
        double m11=m1/50.0f;
        double m12=((double)1.0)-m11;
       
        
        if(rd<=alpha)
            
        {
            
            if(((double)(rand()%1000000)/1000000)<((double)po[pos1]))
                
            {
            if(cell[pos1]==1)
            {if( ((double)(rand()%1000000)/1000000)<(double) pow(m12,r2) )
                
                cell[pos1]=-1;
            else
                cell[pos1]=cell[pos1];}
                
            else if(cell[pos1]==-1)
                
            {if( ((double)(rand()%1000000)/1000000)<(double) pow(m11,r1) )
                
                cell[pos1]=1;
            else
                cell[pos1]=cell[pos1];}
            
            else
                cell[pos1]=cell[pos1];
                
            }
            
            
            else
                
            {
                cell[pos1]=cell[pos1];
            };
            
        }
        
        else
            
        cell[pos1]=cell[pos1]*(-1);
        
        m2=0;
        
        for(int j=0;j<50;j++)
            if(cell[j]==1)
                m2=m2+1;
        
        if(m2<15)
            state2=1;
        if(m2>35)
            state2=2;
        
        if(state2==1)
            statetimesmall=statetimesmall+1;
        else
            statetimelarge=statetimelarge+1;
        
        if(state1==state2)
            transition=transition;
        else
            transition=transition+1;
        
        if(m2<15)
            state1=1;
        if(m2>35)
            state1=2;
        
        
        for(int j=0;j<51;j++)
        {
            if(m2==j)
                M[j]=M[j]+1;
            else
                M[j]=M[j];
        }
        
    
    }
  
    
    
    
    file<<"lattice"<<endl;
    for(int i=1;i<51;i++)
        file<<i<<"\t";
   
    file<<endl;
    
    double Mt=0;
    for(int i=0;i<51;i++)
        Mt=Mt+M[i];
    
    for(int i=0;i<51;i++)
        M[i]=-log(M[i]/Mt);
    
    for(int i=0;i<51;i++)
        file<<M[i]<<"\t";
   
    file<<endl;
    
    //数值积分
    
    
    double R1[1000];
    
    for(int i=0;i<1000;i++)
        R1[i]=alpha*((double)1.0-(((double)i)+((double)0.5))/((double)1000.0))*((double) pow((((double)i)+((double)0.5))/((double)1000.0),r1))+((double)1.0-alpha)*((double)1.0-(((double)i)+((double)0.5))/((double)1000.0));
    
    double R2[1000];
    
    for(int i=0;i<1000;i++)
        R2[i]=alpha*((((double)i)+((double)0.5))/((double)1000.0))*((double) pow((double)1-(((double)i)+((double)0.5))/((double)1000.0),r2))+((double)1.0-alpha)*((((double)i)+((double)0.5))/((double)1000.0));
    

    double drift[1000];
    for(int i=0;i<1000;i++)
        drift[i]=(R1[i]-R2[i]);
    
    
    double var[1000];
    for(int i=0;i<1000;i++)
        var[i]=(R1[i]+R2[i]);
    
    double pvar=0;
    for(int i=0;i<1000;i++)
        pvar=pvar+var[i];
    pvar=pvar/5000000;
    
    double v[1001];
    for(int i=0;i<1001;i++)
        v[i]=(double)0.0;
    
    for(int i=0;i<1001;i++)
        for(int j=0;j<i;j++)
            v[i]=v[i]+((double)2.0)*((double)50.0)*(drift[j]/var[j])/(double)1000.0;
    
    for(int i=0;i<1001;i++)
        v[i]=-(v[i])+log(var[i]/((double)1.0-alpha));
    
    int a=0;
    int b=0;
    int c=0;
    
    for(int i=1;i<1000;i++)
    {if(a==0)
    
    {if(v[i-1]>v[i] && v[i+1]>v[i])
            a=i;
        else
            a=a;}
    else if(v[i-1]>v[i] && v[i+1]>v[i])
        c=i;
    else
        c=c;}
    
    for(int i=1;i<1000;i++)
    {if(v[i-1]<v[i] && v[i+1]<v[i])
            b=i;
        else
            b=b;}
    
    double dva;
    dva=((double)1000000.0)*((v[a+1])-((double)2.0)*(v[a])+(v[a-1]));
    double dvb;
    dvb=((double)1000000.0)*((v[b+1])-((double)2.0)*(v[b])+(v[b-1]));
    double dvc;
    dvc=((double)1000000.0)*((v[c+1])-((double)2.0)*(v[c])+(v[c-1]));
    
    
    
    double vp[51];
    for(int i=0;i<51;i++)
        vp[i]=v[20*i];
    for(int i=0;i<51;i++)
        file<<vp[i]<<"\t";
     file<<endl;
    
    double vt=0.00000000001;
    for(int i=0;i<51;i++)
        vt=vt+exp(-(vp[i]));
    file<<vt;
    file<<endl;
    for(int i=0;i<51;i++)
        vp[i]=(exp(-(vp[i])))/vt;
    file<<endl;
    for(int i=0;i<51;i++)
        vp[i]=-log((vp[i]));
    
    for(int i=0;i<51;i++)
        file<<vp[i]<<"\t";
    file<<endl;
    
    double pertime;
    pertime=((double)loopCount)/((double)transition);
    
    file<<"Pertime"<<"\t"<<pertime;
    file<<endl;
    file<<"Pvar"<<"\t"<<pvar<<endl;
    
    file<<"statetimesmall"<<"\t"<<statetimesmall<<endl;
    file<<"statetimelarge"<<"\t"<<statetimelarge<<endl;
    
    double Tpertimeac;
    Tpertimeac=(((double)1.0)/pvar)*((double)6.28)/(pow((-dva*dvb),0.5))*exp(v[b]-v[a]);
    double Tpertimeca;
    Tpertimeca=(((double)1.0)/pvar)*((double)6.28)/(pow((-dvc*dvb),0.5))*exp(v[b]-v[c]);
    file<<"TheoreticalTac"<<"\t"<<Tpertimeac<<endl;
    file<<"TheoreticalTca"<<"\t"<<Tpertimeca<<endl;
    
    
    file.close();
    
    return 0;
}




