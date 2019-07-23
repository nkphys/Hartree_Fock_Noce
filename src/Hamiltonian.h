#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);
//zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);

class Hamiltonian {
public:

    Hamiltonian(Parameters& Parameters__, Coordinates&  Coordinates__, MFParams& MFParams__ )
        :Parameters_(Parameters__),Coordinates_(Coordinates__),MFParams_(MFParams__)

    {
        Initialize();
        Hoppings();
        HTBCreate();
    }


    void Initialize();    //::DONE
    void Hoppings();        //::DONE
    double GetCLEnergy();    //::DONE
    void InteractionsCreate();   //::DONE
    void Check_Hermiticity();  //::DONE
    void HTBCreate();   //::DONE
    double chemicalpotential(double muin,double Particles);    //::DONE

    double TotalDensity();   //::DONE
    double E_QM();   //::DONE

    void Diagonalize(char option);   //::DONE
    void copy_eigs(int i);  //::DONE


    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    int lx_, ly_, ns_, n_orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx,Ty,Tpxpy,Tpxmy;
    vector<double> eigs_,eigs_saved_,sx_,sy_,sz_;
    Mat_3_Complex_doub Pauli_;
    int X_,Y_,Z_;


    double HS_factor;

};



double Hamiltonian::chemicalpotential(double muin,double Particles){


    double mu_out;
    double n1,N;
    double dMubydN;
    double nstate = eigs_.size();
    dMubydN = 0.0005*(eigs_[nstate-1] - eigs_[0])/nstate;
    N=Particles;
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged=false;
    int final_i;


    if(1==2){
        for(int i=0;i<50000;i++){
            n1=0.0;
            for(int j=0;j<nstate;j++){
                n1+=double(1.0/( exp( (eigs_[j]-mu_out)*Parameters_.beta ) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if(abs(N-n1)<double(0.0001)){
                //cout<<abs(N-n1)<<endl;
                converged=true;
                final_i=i;
                break;
            }
            else {
                mu_out += (N-n1)*dMubydN;
                //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;

            }
        }

        if(!converged){
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else{
            //cout<<"mu converged, N = "<<n1<<" in "<<final_i<<" iters"<<endl;
        }

    }


    double mu1, mu2;
    double mu_temp = muin;
    //cout<<"mu_input = "<<mu_temp<<endl;
    if(1==1){
        mu1=eigs_[0];
        mu2=eigs_[nstate-1];
        for(int i=0;i<40000;i++){
            n1=0.0;
            for(int j=0;j<nstate;j++){
                n1+=double(1.0/( exp( (eigs_[j]-mu_temp)*Parameters_.beta ) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if(abs(N-n1)<double(0.0001)){
                //cout<<abs(N-n1)<<endl;
                converged=true;
                break;
            }
            else {
                if(n1 >N){
                    mu2=mu_temp;
                    mu_temp=0.5*(mu1 + mu_temp);
                }
                else{
                    mu1=mu_temp;
                    mu_temp=0.5*(mu2 + mu_temp);
                }

            }
            //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
        }

        if(!converged){
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else{
            //cout<<"mu converged, N = "<<n1<<endl;
        }

        mu_out = mu_temp;
    }

    return mu_out;
} // ----------


void Hamiltonian::Initialize(){


    ly_=Parameters_.ly;
    lx_=Parameters_.lx;
    ns_=Parameters_.ns;
    n_orbs_=Parameters_.n_orbs;
    X_=0;Y_=1;Z_=2;


    int space=Coordinates_.no_dof_ ;

    HTB_.resize(space,space);
    Ham_.resize(space,space);
    eigs_.resize(space);
    eigs_saved_.resize(space);


    Pauli_.resize(3);

    for(int comp=0;comp<3;comp++){
        Pauli_[comp].resize(2);
        for(int i=0;i<2;i++){
            Pauli_[comp][i].resize(2);
        }
    }

    //Pauli_[X_]
    Pauli_[X_][0][1]=complex<double>(1.0,0.0);
    Pauli_[X_][1][0]=complex<double>(1.0,0.0);

    //Pauli_[Y_]
    Pauli_[Y_][0][1]=complex<double>(0.0,-1.0);
    Pauli_[Y_][1][0]=complex<double>(0.0,1.0);

    //Pauli_[Z_]
    Pauli_[Z_][0][0]=complex<double>(1.0,0.0);
    Pauli_[Z_][1][1]=complex<double>(0.0,-1.0);



} // ----------

double Hamiltonian::TotalDensity(){

    double n1=0.0;
    /*
    for(int j=0;j<eigs_.size();j++){
        n1 +=  1.0f/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    */
    return n1;

} // ----------



double Hamiltonian::E_QM(){

    double E=0.0;
    for(int j=0;j<eigs_.size();j++){
        // E +=  (eigs_[j]-Parameters_.mus)/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
         E +=  (eigs_[j])/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    return E;

} // ----------



double Hamiltonian::GetCLEnergy(){

    double EClassical=0.0;

    for(int site=0;site<ns_;site++) {
        EClassical += (-1.0*(Parameters_.J_Hund - Parameters_.U_prime_onsite)*
                       (
                           (MFParams_.Local_Lz[site]*MFParams_.Local_Lz[site])+
                           (MFParams_.Local_Lx[site]*MFParams_.Local_Lx[site])+
                           (MFParams_.Local_Ly[site]*MFParams_.Local_Ly[site])
                           )
                       )
                +
                (1.0*(Parameters_.U_prime_onsite)*
                 (
                     (MFParams_.Local_Sz[site]*MFParams_.Local_Sz[site])+
                     (MFParams_.Local_Sx[site]*MFParams_.Local_Sx[site])+
                     (MFParams_.Local_Sy[site]*MFParams_.Local_Sy[site])
                     )
                 )
                +
                (-1.0*Parameters_.J_Hund*
                 (n_orbs_/2.0)*((n_orbs_/2.0)+1.0)
                 );

    }

    return EClassical;

} // ----------



void Hamiltonian::InteractionsCreate(){

    int space=Coordinates_.no_dof_;
    int row_, col_;

    for(int i=0;i<space;i++) {
        for(int j=0;j<space;j++) {
            Ham_(i,j)=HTB_(i,j);
        }
    }

    for(int site=0;site<ns_;site++){

        //TERM :: -2*U_prime (<S_vec>.S_vec)
        //(VALUE)*d^{dagger}_{site,orb,s1}d_{site,orb,s2}
        for(int orb=0;orb<n_orbs_;orb++){
            for(int s1=0;s1<2;s1++){
                for(int s2=0;s2<2;s2++){
                    col_=Coordinates_.Nc_dof(site,orb + s2*n_orbs_);
                    row_=Coordinates_.Nc_dof(site,orb + s1*n_orbs_);
                    Ham_(row_,col_) += (-2.0*Parameters_.U_prime_onsite)*
                                        (0.5*
                                        (
                                        (Pauli_[X_][s1][s2]*MFParams_.Local_Sx[site]) +
                                        (Pauli_[Y_][s1][s2]*MFParams_.Local_Sy[site]) +
                                        (Pauli_[Z_][s1][s2]*MFParams_.Local_Sz[site])
                                        )
                                        );

                }
            }
        }

        //TERM :: 2*(J_hund - U_prime) (<L_vec>.L_vec)
        //(VALUE)*d^{dagger}_{site,orb1,s}d_{site,orb2,s}
        for(int orb1=0;orb1<n_orbs_;orb1++){
           for(int orb2=0;orb2<n_orbs_;orb2++){
                for(int s=0;s<2;s++){
                    col_=Coordinates_.Nc_dof(site,orb2 + s*n_orbs_);
                    row_=Coordinates_.Nc_dof(site,orb1 + s*n_orbs_);
                    Ham_(row_,col_) += (2.0*(Parameters_.J_Hund - Parameters_.U_prime_onsite))*
                                        (0.5*
                                        (
                                        (Pauli_[X_][orb1][orb2]*MFParams_.Local_Lx[site]) +
                                        (Pauli_[Y_][orb1][orb2]*MFParams_.Local_Ly[site]) +
                                        (Pauli_[Z_][orb1][orb2]*MFParams_.Local_Lz[site])
                                        )
                                        );

                }
            }
        }

        //TERM :: 3/4*(J_hund + 2*U_prime) (n_{i})
        //(VALUE)*d^{dagger}_{site,orb1,s}d_{site,orb1,s}
        for(int orb=0;orb<n_orbs_;orb++){
                for(int s=0;s<2;s++){
                    col_=Coordinates_.Nc_dof(site,orb + s*n_orbs_);
                    row_=Coordinates_.Nc_dof(site,orb + s*n_orbs_);
                    Ham_(row_,col_) += ((3.0/4.0)*(Parameters_.J_Hund + 2.0*Parameters_.U_prime_onsite));

                }
        }




    }


} // ----------


void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0,0);
    complex<double>temp2;

    for(int i=0;i<Ham_.n_row();i++) {
        for(int j=0;j<Ham_.n_row();j++) {
            if(
                    abs(Ham_(i,j) - conj(Ham_(j,i)))>0.00001
                    ) {
                cout<<Ham_(i,j)<<endl;
                cout<<conj(Ham_(j,i))<<endl;

            }
            assert(
                        abs(Ham_(i,j) - conj(Ham_(j,i)))<0.00001
                        ); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}



void Hamiltonian::Diagonalize(char option){

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);


    char jobz=option;
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


}


void Hamiltonian::HTBCreate(){

    /*Convention
 index = orb + spin*n_orbs_ + site*2*n_orbs_;
orb0_up(site=0),orb1_up(site=0),orb2_up(site=0), orb0_dn(site=0), orb1_dn(site=0), orb2_dn(site=0)....site(n)...
*/
    int mx=Parameters_.TBC_mx;
    int my=Parameters_.TBC_my;
    complex<double> phasex, phasey;
    int l,m,a,b;
    complex<double> Boundary_val;

    if(Parameters_.PBC==true){
        Boundary_val=one_complex;
    }
    else{
        Boundary_val=zero_complex;
    }

    HTB_.fill(0.0);

    for(l=0;l<ns_;l++) {

        // * +x direction Neighbor
        if(Coordinates_.lx_>1){

            if(Coordinates_.indx(l)==(Coordinates_.lx_ -1)){
                phasex=Boundary_val*exp(iota_complex*2.0*(1.0*mx)*PI/(1.0*Parameters_.TBC_cellsX));
                phasey=one_complex;
            }
            else{
                phasex=one_complex;
                phasey=one_complex;
            }
            m = Coordinates_.neigh(l,0);
            for(int spin_=0;spin_<2;spin_++){
                for (int orb1=0;orb1<n_orbs_;orb1++){
                    b=Coordinates_.Nc_dof(l,orb1 + spin_*n_orbs_);
                    for(int orb2=0;orb2<n_orbs_;orb2++){
                        a=Coordinates_.Nc_dof(m,orb2 + spin_*n_orbs_);
                        //value*c^{\dagger}_{a}c_{b}
                        assert (a!=b);
                        if(a!=b){
                            HTB_(a,b)=Parameters_.Hopping_NN(orb2,orb1)*phasex;
                            HTB_(b,a)=conj(HTB_(a,b));
                        }
                    }
                }
            }
        }


        // * +y direction Neighbor
        if(Coordinates_.ly_>1){

            if(Coordinates_.indy(l)==(Coordinates_.ly_ -1)){
                phasex=one_complex;
                phasey=Boundary_val*exp(iota_complex*2.0*(1.0*my)*PI/(1.0*Parameters_.TBC_cellsY));
            }
            else{
                phasex=one_complex;
                phasey=one_complex;
            }
            m = Coordinates_.neigh(l,2);

            for(int spin_=0;spin_<2;spin_++){
                for (int orb1=0;orb1<n_orbs_;orb1++){
                    b=Coordinates_.Nc_dof(l,orb1 + spin_*n_orbs_);
                    for(int orb2=0;orb2<n_orbs_;orb2++){
                        a=Coordinates_.Nc_dof(m,orb2 + spin_*n_orbs_);
                        //value*c^{\dagger}_{a}c_{b}
                        assert (a!=b);
                        if(a!=b){
                            HTB_(a,b)=Parameters_.Hopping_NN(orb2,orb1)*phasey;
                            HTB_(b,a)=conj(HTB_(a,b));
                        }
                    }
                }
            }



        }

    }





    //Crystal_field
    for(int i=0;i<ns_;i++) {
        for(int spin_=0;spin_<2;spin_++){
            for(int orb=0;orb<n_orbs_;orb++) {
                a=Coordinates_.Nc_dof(i,orb + n_orbs_*spin_);
                HTB_(a,a)+=(Parameters_.Crystal_Field[orb]) + MFParams_.Disorder(Coordinates_.indx(i),Coordinates_.indy(i));
            }
        }
    }


    //HTB_.print();
} // ----------



void Hamiltonian::Hoppings(){
    //DOES SOMETHING EXACT i.e NOTHING :)

} // ----------

void Hamiltonian::copy_eigs(int i){

    int space=2*ns_;

    if (i == 0) {
        for(int j=0;j<space;j++) {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else {
        for(int j=0;j<space;j++) {
            eigs_saved_[j] = eigs_[j];
        }
    }

}


#endif
