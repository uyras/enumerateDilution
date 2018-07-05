/*
        kagome_enumerate.cpp

        Wang-Landau method for diluted Ising model (pyrochlore spin ice)

           programmed by YO          2008/04/07  (original for square lattice)
           programmed by YO          2016/07/18

        nla=nx*ny     : number of lattice sites

        g[E+3*N]        : log of DOS
               -3*N <= E <= N  (asymmetric because of frustration)
        g[ie]  0 <= ie <= nemax (= 4 * N)

        n      : number of iteration (0-nfinal)
        nfinal : final number of n   (default 24)
        f      : modification factor for g[ie] (1/2**n)
        factor : criterion of flat   (default 0.80)

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <bitset>
#include <string>

void period();
void spinset();
void bruteforce();
int energyTest();
void check();

#define nx 4
#define ny 4
#define nla (nx*ny)
#define hmax 6
#define dividor 100
#define nemax ((4*nla) + 2*nla*hmax)*dividor
#define neigh 6 //число ближайших соседей

int isp[nla];
int permute[nla];
int nn[neigh*nla];

int g[nemax+1];

int energy;

int iflag;

int nhole,nspins;

int h=0;

float hd = 0;

using namespace std;

int main(void)
{

    cout<<"N="<<nla<<endl;
    cout<<"input number of hole (zero) spins: ";
    cin>>nhole;
    //nhole=0;

    string hds;
    cout<<"input external field value H (float), which can be multiplied by "<<dividor<<" without residue: ";
    cin>>hds;
    hd = stof(hds);
    h=hd*dividor;

    printf("# DOS enumeration of Ising model ");
    printf("on the triangular lattice\n");


    printf("# linear size of triangular     L= %2d,  number of sites  N= %d,  field H= %f\n",nx,nla,(h/(double)dividor));
    printf("# concentration of holes  %5.3f,  number of holes  %d,  number of spins  %d\n",(double)nhole/nla,nhole,nla-nhole);

    nspins=nla-nhole;
    int total=1<<nspins; //total number of states

    for (int i=0; i<nla; ++i){
        if (i<nhole)
            permute[i]=0;
        else
            permute[i]=1;
    }

    period();

    /*ofstream of ("graph.dot");
    of<<"graph G {"<<endl;
    for (int i=0;i<N;i++){
        of<<i<<" -- "<<nn[i][0]<<";"<<endl;
        of<<i<<" -- "<<nn[i][1]<<";"<<endl;
        of<<i<<" -- "<<nn[i][2]<<";"<<endl;
        of<<i<<" -- "<<nn[i][3]<<";"<<endl;
        of<<i<<" -- "<<nn[i][4]<<";"<<endl;
        of<<i<<" -- "<<nn[i][5]<<";"<<endl;
    }
    of<<"}"<<endl;
    of.close();*/


    char fname[100];
    /*sprintf(fname,"g_%d.dat",nhole);
    ofstream f(fname);
    f<<"# number of spins N= "<<nspins<<endl;
    f<<"# E\tE/N\tg[E]\tg[E]/N"<<endl;*/

    int permNum=0;


    sprintf(fname,"g_%d_%s.dat",nhole,hds.c_str());
    ofstream f1(fname);
    f1<<"# number of spins N= "<<nspins<<endl;
    f1<<"# E\tE/N\tg[E]\tg[E]/N\tH"<<endl;

    do {
        spinset(); //устанавливаем спины в нужное значение
        /*cout<<"# permutation: ";
        for (int i=0; i<nla; ++i) cout<<permute[i];
        cout<<"; initial energy = "<<energy<<endl;*/
        bruteforce();
        //cout<<"next round"<<endl;

        for(int i=0; i<=nemax; i++){
            if (g[i] != 0) {
                /*f<<
                    -(i-3*N)<<"\t"<<
                    -(double)(i-3*N)/N<<"\t"<<
                    g[i]<<"\t"<<
                    g[i]/(double)(total)<<"\t"<<
                    endl;*/
                f1<<
                    (i/(double)dividor)-nla-hmax*nla<<"\t"<<
                    ((i/(double)dividor)-nla-hmax*nla)/nla<<"\t"<<
                    g[i]<<"\t"<<
                    g[i]/(double)(total)<<"\t"<<
                    hd<<"\t"<<
                    endl;
                break;
            }
        }

        permNum++;
    } while ( std::next_permutation(permute,permute+nla) );
    f1.close();

    //f.close();
    cout<<"# done with file"<<endl;
}

void period()
/*
       periodic boundary conditions for kagome lattice
                   system size = nx*ny
*/
{
    int vars[] = {
        13,14,4,2,5,6,
        14,15,1,3,6,7,
        15,16,2,4,7,8,
        16,13,3,1,8,5,
        12,9,8,6,1,4,
        9,10,5,7,2,1,
        10,11,6,8,3,2,
        11,12,7,5,4,3,
        13,14,12,10,6,5,
        14,15,9,11,7,6,
        15,16,10,12,8,7,
        16,13,11,9,5,8,
        4,1,16,14,9,12,
        1,2,13,15,10,9,
        2,3,14,16,11,10,
        3,4,15,13,12,11
    };
    for (int i=0;i<nla*neigh;i++){
        nn[i]=vars[i]-1;
    }
}

void check(){
    int res[nla];

    for (int i=0;i<nla;i++){
        res[i]=0;
    }


    for (int i=0;i<nla*neigh;i++){
        res[nn[i]]+=1;
    }


    for (int i=0;i<nla;i++){
        cout<<i+1<<"="<<res[i]<<endl;
    }
}

void spinset()
{

    for (int i=0; i<=nla-1; ++i){
        isp[i]=permute[i];
    }

    int la;

    energy=0;
    for (la=0; la < nla; la++){
        energy += isp[la]*(
                isp[nn[la*neigh+0]]+
                isp[nn[la*neigh+1]]+
                isp[nn[la*neigh+2]]+
                isp[nn[la*neigh+3]]+
                isp[nn[la*neigh+4]]+
                isp[nn[la*neigh+5]]
                )*dividor;
    }
    energy /= 2;

    for (la=0; la < nla; la++){
        energy += isp[la]*h;
    }

    for(int i=0; i<=nemax; i++){
        g[i]=0;
    }
    //cout<<energy<<endl;
}

void bruteforce()
{
    unsigned gray,gray_old=0,mask,temp,nsteps=1<<nspins;

    for (unsigned i=0;i<nsteps;++i){
        //cout<<energy<<"|"<<energy+(nla+hmax*nla)*dividor<<"|"<<energy+(nla+0*nla)*dividor<<"|"<<endl;

        gray = i ^ (i >> 1); //Gray's code. To enumerate all states, but with only 1 spins difference between next and previous state
        mask = gray ^ gray_old; //Evaluate which spin should be changed

        //std::cout << "gray = " << std::bitset<32>(gray)  << std::endl;
        //std::cout << "mask = " << std::bitset<32>(mask)  << std::endl;
        temp=mask;
        for (int j=0;j<nla;j++){
            if (isp[j]==0)
                continue;

            if (temp==1){
                isp[j]*=-1;
                energy += 2*isp[j]*(
                            isp[nn[j*neigh+0]]+
                            isp[nn[j*neigh+1]]+
                            isp[nn[j*neigh+2]]+
                            isp[nn[j*neigh+3]]+
                            isp[nn[j*neigh+4]]+
                            isp[nn[j*neigh+5]]
                            )*dividor;
                energy += 2*isp[j]*h;
                //cout<<energy<<"|"<<energyTest()<<endl;
                break;
            }

            temp = temp>>1;
        }

        //update the dos with energy
        g[energy+(nla+hmax*nla)*dividor]+=1;

        gray_old=gray;
    }
}

int energyTest(){
    int tenergy=0;
    for (int la=0; la < nla; la++){
        tenergy += isp[la]*(
                    isp[nn[la*neigh+0]]+
                    isp[nn[la*neigh+1]]+
                    isp[nn[la*neigh+2]]+
                    isp[nn[la*neigh+3]]+
                    isp[nn[la*neigh+4]]+
                    isp[nn[la*neigh+5]]
                    )*dividor;
    }
    tenergy /= 2;

    for (int la=0; la < nla; la++){
        tenergy += isp[la]*h;
    }
    return tenergy;
}
