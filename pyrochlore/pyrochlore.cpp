/*
        kagome_enumerate.cpp

        Wang-Landau method for diluted Ising model (pyrochlore spin ice)

           programmed by YO          2008/04/07  (original for square lattice)
           programmed by YO          2016/07/18

        nla = nx*ny*nz  : number of cubic unit cells
                16 sites in single cubic lattice
        N = 16*nla      : number of total spins

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

void period_sc();
void period();
void spinset();
void bruteforce();
int energyTest();
void check();

#define nx 1
#define ny nx
#define nz nx
#define nxy (nx*ny)
#define nla (nx*ny*nz)
#define N (16*nla)
#define hmax 6
#define dividor 150
#define nemax ((4*N) + 2*N*hmax)*dividor

int isp[N];
int permute[N];
int nn_sc[6*nla];
int nn[16*nla][6];

int g[nemax+1];

int energy;

int iflag;

int nhole,nspins;

int ha=0; // поле от спинов * 1/3
int hb=0; // полноценное поле от спинов

float hd = 0;

using namespace std;

int main(void)
{

    cout<<"N="<<N<<endl;
    cout<<"input number of hole (zero) spins: ";
    cin>>nhole;
    //nhole=0;

    string hds;
    cout<<"input external field value H (float), which can be multiplied by "<<dividor<<" without residue: ";
    cin>>hds;
    hd = stof(hds);
    ha=hd*(dividor/3);
    hb=hd*dividor;

    printf("# DOS enumeration of Ising model ");
    printf("on the pyrochlore lattice\n");


    printf("# linear size of cube     L= %2d,  number of sites  N= %d,  field H= %f\n",nx,N,hd);
    printf("# concentration of holes  %5.3f,  number of holes  %d,  number of spins  %d\n",(double)nhole/N,nhole,N-nhole);

    nspins=N-nhole;
    int total=1<<nspins; //total number of states

    for (int i=0; i<N; ++i){
        if (i<nhole)
            permute[i]=0;
        else
            permute[i]=1;
    }

    period_sc();
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
        for (int i=0; i<N; ++i) cout<<permute[i];
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
                    (i/(double)dividor)-N-hmax*N<<"\t"<<
                    ((i/(double)dividor)-N-hmax*N)/N<<"\t"<<
                    g[i]<<"\t"<<
                    g[i]/(double)(total)<<"\t"<<
                    hd<<"\t"<<
                    endl;
                break;
            }
        }

        permNum++;
    } while ( std::next_permutation(permute,permute+N) );
    f1.close();

    //f.close();
    cout<<"# done with file"<<endl;
}

void period_sc()
/*
       periodic boundary conditions for 3-d simple cubic lattice
                   system size = nx*ny*nz
*/
{
    int la, ix, ixy;

    for (la=0; la <= nla-1; la++){
        ix=((int)(la/nx))*nx;
        ixy=((int)(la/nxy))*nxy;
        nn_sc[la]       = (la+1)%nx   +ix;
        nn_sc[la+nla]   = (la-1+nx)%nx+ix;
        nn_sc[la+2*nla] = (la+nx) % nxy + ixy;
        nn_sc[la+3*nla] = (la-nx+nxy)% nxy + ixy;
        nn_sc[la+4*nla] = (la+nxy) % nla;
        nn_sc[la+5*nla] = (la-nxy+nla)% nla;
    }
}

void period()
/*
       periodic boundary conditions for pyrochlore lattice
       16 sites in a unit cubic cell
*/
{
    int la;

    for (la=0; la <= nla-1; la++){
        //=== 0
        nn[16*la][0] = la*16+1;
        nn[16*la][1] = la*16+4;
        nn[16*la][2] = la*16+7;
        nn[16*la][3] = nn_sc[nn_sc[la+3*nla]+5*nla]*16+3;
        nn[16*la][4] = nn_sc[nn_sc[la+nla]+5*nla]*16+6;
        nn[16*la][5] = nn_sc[nn_sc[la+nla]+3*nla]*16+9;

        //=== 1
        nn[16*la+1][3] = la*16;
        nn[16*la+1][5] = la*16+4;
        nn[16*la+1][4] = la*16+7;
        nn[16*la+1][0] = la*16+2;
        nn[16*la+1][1] = nn_sc[la+nla]*16+13;
        nn[16*la+1][2] = nn_sc[la+nla]*16+14;

        //=== 2
        nn[16*la+2][0] = la*16+3;
        nn[16*la+2][1] = la*16+10;
        nn[16*la+2][2] = la*16+11;
        nn[16*la+2][3] = la*16+1;
        nn[16*la+2][4] = nn_sc[la+nla]*16+13;
        nn[16*la+2][5] = nn_sc[la+nla]*16+14;

        //=== 3
        nn[16*la+3][3] = la*16+2;
        nn[16*la+3][5] = la*16+10;
        nn[16*la+3][4] = la*16+11;
        nn[16*la+3][0] = nn_sc[nn_sc[la+2*nla]+4*nla]*16;
        nn[16*la+3][1] = nn_sc[nn_sc[la+nla]+2*nla]*16+6;
        nn[16*la+3][2] = nn_sc[nn_sc[la+nla]+4*nla]*16+9;

        //=== 4
        nn[16*la+4][5] = la*16;
        nn[16*la+4][0] = la*16+1;
        nn[16*la+4][1] = la*16+7;
        nn[16*la+4][2] = la*16+5;
        nn[16*la+4][4] = nn_sc[la+3*nla]*16+10;
        nn[16*la+4][3] = nn_sc[la+3*nla]*16+15;

        //=== 5
        nn[16*la+5][0] = la*16+6;
        nn[16*la+5][1] = la*16+12;
        nn[16*la+5][2] = la*16+13;
        nn[16*la+5][3] = la*16+4;
        nn[16*la+5][5] = nn_sc[la+3*nla]*16+10;
        nn[16*la+5][4] = nn_sc[la+3*nla]*16+15;

        //=== 6
        nn[16*la+6][5] = la*16+5;
        nn[16*la+6][0] = la*16+12;
        nn[16*la+6][1] = la*16+13;
        nn[16*la+6][2] = nn_sc[nn_sc[la]+4*nla]*16;
        nn[16*la+6][3] = nn_sc[nn_sc[la+3*nla]]*16+3;
        nn[16*la+6][4] = nn_sc[nn_sc[la+3*nla]+4*nla]*16+9;

        //=== 7
        nn[16*la+7][4] = la*16;
        nn[16*la+7][0] = la*16+1;
        nn[16*la+7][5] = la*16+4;
        nn[16*la+7][1] = la*16+8;
        nn[16*la+7][2] = nn_sc[la+5*nla]*16+11;
        nn[16*la+7][3] = nn_sc[la+5*nla]*16+12;

        //=== 8
        nn[16*la+8][0] = la*16+9;
        nn[16*la+8][1] = la*16+14;
        nn[16*la+8][2] = la*16+15;
        nn[16*la+8][3] = la*16+7;
        nn[16*la+8][4] = nn_sc[la+5*nla]*16+11;
        nn[16*la+8][5] = nn_sc[la+5*nla]*16+12;

        //=== 9
        nn[16*la+9][4] = la*16+8;
        nn[16*la+9][5] = la*16+14;
        nn[16*la+9][0] = la*16+15;
        nn[16*la+9][1] = nn_sc[nn_sc[la]+2*nla]*16;
        nn[16*la+9][3] = nn_sc[nn_sc[la+5*nla]]*16+3;
        nn[16*la+9][2] = nn_sc[nn_sc[la+5*nla]+2*nla]*16+6;

        //=== 10
        nn[16*la+10][5] = la*16+2;
        nn[16*la+10][0] = la*16+3;
        nn[16*la+10][4] = la*16+11;
        nn[16*la+10][3] = la*16+15;
        nn[16*la+10][1] = nn_sc[la+2*nla]*16+4;
        nn[16*la+10][2] = nn_sc[la+2*nla]*16+5;

        //=== 11
        nn[16*la+11][5] = la*16+2;
        nn[16*la+11][0] = la*16+3;
        nn[16*la+11][1] = la*16+10;
        nn[16*la+11][3] = la*16+12;
        nn[16*la+11][4] = nn_sc[la+4*nla]*16+7;
        nn[16*la+11][2] = nn_sc[la+4*nla]*16+8;

        //=== 12
        nn[16*la+12][5] = la*16+5;
        nn[16*la+12][3] = la*16+6;
        nn[16*la+12][0] = la*16+11;
        nn[16*la+12][4] = la*16+13;
        nn[16*la+12][1] = nn_sc[la+4*nla]*16+7;
        nn[16*la+12][2] = nn_sc[la+4*nla]*16+8;

        //=== 13
        nn[16*la+13][5] = la*16+5;
        nn[16*la+13][4] = la*16+6;
        nn[16*la+13][0] = la*16+12;
        nn[16*la+13][1] = la*16+14;
        nn[16*la+13][3] = nn_sc[la]*16+1;
        nn[16*la+13][2] = nn_sc[la]*16+2;

        //=== 14
        nn[16*la+14][5] = la*16+8;
        nn[16*la+14][0] = la*16+9;
        nn[16*la+14][3] = la*16+13;
        nn[16*la+14][1] = la*16+15;
        nn[16*la+14][4] = nn_sc[la]*16+1;
        nn[16*la+14][2] = nn_sc[la]*16+2;

        //=== 15
        nn[16*la+15][5] = la*16+8;
        nn[16*la+15][3] = la*16+9;
        nn[16*la+15][4] = la*16+14;
        nn[16*la+15][0] = la*16+10;
        nn[16*la+15][1] = nn_sc[la+2*nla]*16+4;
        nn[16*la+15][2] = nn_sc[la+2*nla]*16+5;
    }

}

void check(){
    int res[N];

    for (int i=0;i<N;i++){
        res[i]=0;
    }


    for (int i=0;i<N;i++){
        cout<<i<<" - ";
        for (int j=0; j<6; ++j) {
            res[nn[i][j]]+=1;
            cout<<nn[i][j]<<",";
        }
        cout<<endl;
    }

    cout<<"----------------------------------"<<endl;


    for (int i=0;i<N;i++){
        cout<<i+1<<"="<<res[i]<<endl;
    }
}

void spinset()
{

    for (int i=0; i<=N-1; ++i){
        isp[i]=permute[i];
    }

    int la;

    energy=0;
    for (la=0; la < N; la++){
        energy += isp[la]*(isp[nn[la][0]]+isp[nn[la][1]]
                +isp[nn[la][2]])*dividor;
    }

    for (la=0; la < N; la++){
        if (la == 1 || la == 3 || la == 12 || la == 15)
            energy += isp[la]*hb;
        else
            energy += isp[la]*ha;
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
        //cout<<energy<<"|"<<energy+(N+hmax*N)*dividor<<"|"<<energy+(N+0*N)*dividor<<"|"<<endl;

        gray = i ^ (i >> 1); //Gray's code. To enumerate all states, but with only 1 spins difference between next and previous state
        mask = gray ^ gray_old; //Evaluate which spin should be changed

        //std::cout << "gray = " << std::bitset<32>(gray)  << std::endl;
        //std::cout << "mask = " << std::bitset<32>(mask)  << std::endl;
        temp=mask;
        for (int j=0;j<N;j++){
            if (isp[j]==0)
                continue;

            if (temp==1){
                isp[j]*=-1;
                energy += 2*isp[j]*(isp[nn[j][0]]+isp[nn[j][1]]
                        +isp[nn[j][2]]+isp[nn[j][3]]
                        +isp[nn[j][4]]+isp[nn[j][5]])*dividor;

                if (j == 1 || j == 3 || j == 12 || j == 15)
                    energy += 2*isp[j]*hb;
                else
                    energy += 2*isp[j]*ha;
                //cout<<energy<<"|"<<energyTest()<<endl;
                break;
            }

            temp = temp>>1;
        }

        //update the dos with energy
        g[energy+(N+hmax*N)*dividor]+=1;

        gray_old=gray;
    }
}

int energyTest(){
    int tenergy=0;
    for (int la=0; la < N; la++){
        tenergy += isp[la]*(isp[nn[la][0]]+isp[nn[la][1]]
                        +isp[nn[la][2]]+isp[nn[la][3]]
                        +isp[nn[la][4]]+isp[nn[la][5]])*dividor;
    }
    tenergy /= 2;

    for (int la=0; la < N; la++){
        if (la == 1 || la == 3 || la == 12 || la == 15)
            tenergy += isp[la]*hb;
        else
            tenergy += isp[la]*ha;
    }
    return tenergy;
}
