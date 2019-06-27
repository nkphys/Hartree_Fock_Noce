#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "tensor_type.h"
#include "functions.h"

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#define PI acos(-1.0)

//n, a, lda, ipvt, work, lwork, info
extern "C" void   zgetri_(int *,std::complex<double> *, int *, int *, std::complex<double> *, int *, int *);

//zgetrf_ (&n, &n, &(_TEMP(0,0)), &n, &(ipvt[0]), &info);
extern "C" void   zgetrf_(int *,int *, std::complex<double> *, int *, int *, int *);


//zhetri (character UPLO, integer N, complex*16 dimension( lda, * ) A, integer LDA,
//integer IPIV, complex*16 dimension( * ) WORK, integer INFO)
extern "C" void   zhetri_(char *, int *, std::complex<double> *, int *, int *, std::complex<double> *, int *);


//zhetrf(uplo,n,a,lda,ipiv,work,lwork,info)
extern "C" void   zhetrf_(char *, int *, std::complex<double> *, int *, int *, std::complex<double> *, int *, int *);



class Observables{
public:

    Observables(Parameters& Parameters__, Coordinates& Coordinates__,
                MFParams& MFParams__, Hamiltonian& Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ns_(Parameters_.ns),n_orbs_(Parameters_.n_orbs)
    {
        Initialize();
    }

    void Initialize();

    void Calculate_Akw();
    void Calculate_Nw();
    void Calculate_IPR();
    void Calculate_Optical_Conductivity();
    void Get_Non_Interacting_dispersion();
    double Lorentzian(double x, double brd);
    complex<double> DOT_P(Mat_1_Complex_doub left, Mat_1_Complex_doub right);

    void Total_Energy_Average(double Curr_QuantE, double CurrE);

    void Calculate_Order_Params();
    void Calculate_Local_n_orb_resolved();
    void Get_OrderParameters_diffs();
    void Update_OrderParameters(int iter);
    void Update_OrderParameters_Modified_Broyden_in_ComplexSpace(int iter_count);
    void Invert_Beta();

    void Calculate_Single_Particle_Density_Matrix();
    complex<double> Two_particle_Den_Mat(int _alpha, int _beta, int _gamma, int _delta);
    void Calculate_two_point_correlations();


    Parameters& Parameters_;
    Coordinates& Coordinates_;
    MFParams& MFParams_;
    Hamiltonian& Hamiltonian_;
    int lx_,ly_,ns_,n_orbs_;
    double dosincr_,tpi_;
    vector<double> nia_,nib_,nic_;
    Matrix<double> SiSj_,dos;
    Mat_1_doub Local_Sz_, Local_Sx_, Local_Sy_;
    Mat_1_doub Local_Lz_, Local_Lx_, Local_Ly_;
    Mat_1_doub Local_n_orb_resolved;

    Mat_2_Complex_doub SP_Density_Matrix;


    // Declare Fields
    Matrix<double> Sz_obs, Sx_obs, Sy_obs;
    Matrix<double> Local_density_obs;
    double Error_OP_;
    double AVG_Total_Energy, AVG_Total_Energy_sqr;


    //Declare Broyden_Mixing vectors
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;


    //Declarations for Modified Broyden Method in Complex Space
    Mat_2_Complex_doub _Delta_F;
    Mat_2_Complex_doub _u;
    Mat_2_Complex_doub _A;
    Mat_1_Complex_doub _Fm, _Fm_minus1;
    Mat_1_Complex_doub _Delta_OPm, _Delta_OPm_minus1;
    Matrix<complex<double>> _Beta;
    Mat_1_Complex_doub _cm, _gammam;
    double w_minus1;
    Mat_1_doub w;




};
/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/



void Observables::Invert_Beta(){

    int n=_Beta.n_row();
    int lda=_Beta.n_col();
    int info;
    vector<int> ipvt;
    ipvt.resize(_Beta.n_col());
    vector<complex<double>> work(n);
    int lwork= n;


    char uplo='U';
    zhetrf_(&uplo, &n, &(_Beta(0,0)),&lda,&(ipvt[0]),&(work[0]),&lwork,&info);
    //cout<<"FACTORIZATION OF MATRIX:"<<endl;
    //_TEMP.print();

    zhetri_(&uplo, &n, &(_Beta(0,0)),&lda,&(ipvt[0]),&(work[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("Inverse: zgetri: failed with info!=0.\n");
    }

    //cout<<"INVERSE OF MATRIX:"<<endl;
    //_TEMP.print();

}



complex<double> Observables::DOT_P(Mat_1_Complex_doub left, Mat_1_Complex_doub right){
    complex<double> temp_;
    temp_=zero_complex;

    assert(left.size()==right.size());

    for(int i=0;i<left.size();i++){
        temp_ += conj(left[i])*right[i];
    }
    return temp_;

}

void Observables::Get_OrderParameters_diffs(){

    Error_OP_=0.0;

    for(int site=0;site<ns_;site++){
        Error_OP_ += abs((Local_Sz_[site] - MFParams_.Local_Sz[site])*
                         (Local_Sz_[site] - MFParams_.Local_Sz[site]));
        Error_OP_ += abs((Local_Sx_[site] - MFParams_.Local_Sx[site])*
                         (Local_Sx_[site] - MFParams_.Local_Sx[site]));
        Error_OP_ += abs((Local_Sy_[site] - MFParams_.Local_Sy[site])*
                         (Local_Sy_[site] - MFParams_.Local_Sy[site]));

        Error_OP_ += abs((Local_Lz_[site] - MFParams_.Local_Lz[site])*
                         (Local_Lz_[site] - MFParams_.Local_Lz[site]));
        Error_OP_ += abs((Local_Lx_[site] - MFParams_.Local_Lx[site])*
                         (Local_Lx_[site] - MFParams_.Local_Lx[site]));
        Error_OP_ += abs((Local_Ly_[site] - MFParams_.Local_Ly[site])*
                         (Local_Ly_[site] - MFParams_.Local_Ly[site]));
    }

    Error_OP_=sqrt(Error_OP_);

}

void Observables::Update_OrderParameters(int iter){


    //Simple mixing
    double alpha_OP=Parameters_.alpha_OP;

    if(Parameters_.Simple_Mixing==true){

        if(iter==0){
            cout<<"Using Simple Mixing to gain Self-Consistency"<<endl;
        }

        for(int site=0;site<ns_;site++){
            MFParams_.Local_Sz[site] = (1-alpha_OP)*MFParams_.Local_Sz[site]
                    + alpha_OP*Local_Sz_[site];

            MFParams_.Local_Sx[site] = (1-alpha_OP)*MFParams_.Local_Sx[site]
                    + alpha_OP*Local_Sx_[site];
            //MFParams_.Local_Sx[site] = 0.0;

            MFParams_.Local_Sy[site] = (1-alpha_OP)*MFParams_.Local_Sy[site]
                    + alpha_OP*Local_Sy_[site];
            // MFParams_.Local_Sy[site] = 0.0;

            MFParams_.Local_Lz[site] = (1-alpha_OP)*MFParams_.Local_Lz[site]
                    + alpha_OP*Local_Lz_[site];

            MFParams_.Local_Lx[site] = (1-alpha_OP)*MFParams_.Local_Lx[site]
                    + alpha_OP*Local_Lx_[site];
            // MFParams_.Local_Lx[site] = 0.0;

            MFParams_.Local_Ly[site] = (1-alpha_OP)*MFParams_.Local_Ly[site]
                    + alpha_OP*Local_Ly_[site];
            //MFParams_.Local_Ly[site] = 0.0;

        }

        //      Parameters_.mu_old = (1-alpha_OP)*Parameters_.mu_old + alpha_OP*Parameters_.mus;

    }



    /*
    //Declared in initialize function();
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;
     */

    vector<double> vec_V, vec_U;
    double Denominator_;
    vector<double> vec_L;
    int site;


    if(Parameters_.Broyden_Mixing==true){
        //assert(false);

        if(iter==0){
            cout<<"Using Broyden Mixing to gain Self-Consistency"<<endl;


            //Get Jinv_np1
            for(int i=0;i<(6*ns_);i++){
                for(int j=0;j<(6*ns_);j++){
                    if(i==j){
                        Jinv_np1[i][j]=-1.0*alpha_OP;
                    }
                    else{
                        Jinv_np1[i][j]=0.0;
                    }
                }}

            //Get F_n
            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site=Coordinates_.Nc(i,j);
                    F_n[site + ns_*0]=Local_Sz_[site] - MFParams_.Local_Sz[site];
                    F_n[site + ns_*1]=Local_Sx_[site] - MFParams_.Local_Sx[site];
                    F_n[site + ns_*2]=Local_Sy_[site] - MFParams_.Local_Sy[site];
                    F_n[site + ns_*3]=Local_Lz_[site] - MFParams_.Local_Lz[site];
                    F_n[site + ns_*4]=Local_Lx_[site] - MFParams_.Local_Lx[site];
                    F_n[site + ns_*5]=Local_Ly_[site] - MFParams_.Local_Ly[site];
                }
            }


            for(int i=0;i<(6*ns_);i++){
                Delta_x_n[i] =0.0;
                for(int j=0;j<(6*ns_);j++){
                    Delta_x_n[i] +=  -1.0*Jinv_np1[i][j]*F_n[j];
                }
            }



            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site=Coordinates_.Nc(i,j);
                    MFParams_.Local_Sz[site] +=Delta_x_n[site + ns_*0];
                    MFParams_.Local_Sx[site] +=Delta_x_n[site + ns_*1];
                    MFParams_.Local_Sy[site] +=Delta_x_n[site + ns_*2];
                    MFParams_.Local_Lz[site] +=Delta_x_n[site + ns_*3];
                    MFParams_.Local_Lx[site] +=Delta_x_n[site + ns_*4];
                    MFParams_.Local_Ly[site] +=Delta_x_n[site + ns_*5];
                }
            }

            //Copy Jinv_np1 to Jinv_n
            Jinv_n = Jinv_np1;

            //Copy F_n to F_nm1
            F_nm1=F_n;
        }

        else{
            //Get F_n
            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site=Coordinates_.Nc(i,j);
                    F_n[site + ns_*0]=Local_Sz_[site] - MFParams_.Local_Sz[site];
                    F_n[site + ns_*1]=Local_Sx_[site] - MFParams_.Local_Sx[site];
                    F_n[site + ns_*2]=Local_Sy_[site] - MFParams_.Local_Sy[site];
                    F_n[site + ns_*3]=Local_Lz_[site] - MFParams_.Local_Lz[site];
                    F_n[site + ns_*4]=Local_Lx_[site] - MFParams_.Local_Lx[site];
                    F_n[site + ns_*5]=Local_Ly_[site] - MFParams_.Local_Ly[site];
                }
            }

            //Get DeltaF_n
            for (int i=0;i<(6*ns_);i++){
                DeltaF_n[i] = F_n[i] - F_nm1[i];
            }

            //Get vec_V = Jinv_n*DeltaF_n
            vec_V.clear();
            vec_V.resize(6*ns_);

            for(int i=0;i<6*ns_;i++){
                vec_V[i] =0.0;
                for(int j=0;j<6*ns_;j++){
                    vec_V[i] += Jinv_n[i][j]*DeltaF_n[j];
                }
            }

            //Get vec_U = Delta_x_n^dagg*Jinv_n
            vec_U.clear();
            vec_U.resize(6*ns_);

            for(int i=0;i<6*ns_;i++){
                vec_U[i] =0.0;
                for(int j=0;j<6*ns_;j++){
                    vec_U[i] += Delta_x_n[j]*Jinv_n[j][i];
                }
            }

            // Get Denominator_=<Delta_x_n|vec_V>
            Denominator_=0.0;
            for(int i=0;i<6*ns_;i++){
                Denominator_ +=Delta_x_n[i]*vec_V[i];
            }


            //Get vec_L=  Delta_x_n - vec_V;
            vec_L.clear();
            vec_L.resize(6*ns_ );
            for(int i=0;i<6*ns_;i++){
                vec_L[i] = Delta_x_n[i] - vec_V[i];
            }


            //Get Mat_Temp [Remember to clear later on];
            Mat_2_doub Mat_Temp;
            Mat_Temp.resize(6*ns_ );
            for(int i=0;i<6*ns_ ;i++){
                Mat_Temp[i].resize(6*ns_);
                for(int j=0;j<6*ns_ ;j++){
                    Mat_Temp[i][j] = (vec_L[i]*vec_U[j])/(Denominator_);
                }
            }


            //Get Jinv_np1

            for(int i=0;i<6*ns_;i++){
                for(int j=0;j<6*ns_;j++){
                    Jinv_np1[i][j]  = Jinv_n[i][j]  + Mat_Temp[i][j];
                }
            }

            for(int i=0;i<6*ns_;i++){
                Delta_x_n[i] =0.0;
                for(int j=0;j<6*ns_;j++){
                    Delta_x_n[i] +=  -1.0*Jinv_np1[i][j]*F_n[j];
                }
            }



            for(int i=0;i<lx_;i++){
                for(int j=0;j<ly_;j++){
                    site=Coordinates_.Nc(i,j);
                    MFParams_.Local_Sz[site] +=Delta_x_n[site + ns_*0];
                    MFParams_.Local_Sx[site] +=Delta_x_n[site +  ns_*1];
                    MFParams_.Local_Sy[site] +=Delta_x_n[site +  ns_*2];
                    MFParams_.Local_Lz[site] +=Delta_x_n[site +  ns_*3];
                    MFParams_.Local_Lx[site] +=Delta_x_n[site +  ns_*4];
                    MFParams_.Local_Ly[site] +=Delta_x_n[site +  ns_*5];
                }
            }


            //Copy Jinv_np1 to Jinv_n
            Jinv_n = Jinv_np1;

            //Copy F_n to F_nm1
            F_nm1=F_n;


            //Clear Mat_Temp
            for(int i=0;i<6*ns_;i++){
                Mat_Temp[i].clear();
            }
            Mat_Temp.clear();

        }



        //Checking MFParams_.OParams symmetries

        //        for(int i=0;i<lx_;i++){
        //            for(int j=0;j<ly_;j++){
        //                site=Coordinates_.Nc(i,j);
        //                for(int state_i=0;state_i<6;state_i++){
        //                    for(int state_j=0;state_j<6;state_j++){
        //                    if(state_i==state_j){
        //                    if(MFParams_.OParams[site][state_i][state_j].imag()!=0.0){
        //                        cout<<"local densities order params have problem"<<endl;
        //                    }
        //                      assert(MFParams_.OParams[site][state_i][state_j].imag()==0.0);
        //                    }
        //                    else{
        //                    if(
        //                      abs(MFParams_.OParams[site][state_i][state_j] - conj(MFParams_.OParams[site][state_j][state_i]))>0.00001
        //                       ){
        //                        cout<<"Hermiticity of order params is an issue"<<endl;

        //                        cout<<site<<"   "<<state_i<<"    "<<state_j<<"    "<<MFParams_.OParams[site][state_i][state_j]<<endl;
        //                        cout<<site<<"   "<<state_j<<"    "<<state_i<<"    "<<MFParams_.OParams[site][state_j][state_i]<<endl;
        //                    }
        //                    assert(
        //                          abs(MFParams_.OParams[site][state_i][state_j] - conj(MFParams_.OParams[site][state_j][state_i]))<0.00001
        //                          );
        //                    }

        //                    }}}}



        //Maintain hermiticity of order params, because after many iterations they tend to loose it
        //        for(int i=0;i<lx_;i++){
        //            for(int j=0;j<ly_;j++){
        //                site=Coordinates_.Nc(i,j);
        //                for(int state_i=0;state_i<6;state_i++){
        //                    for(int state_j=state_i;state_j<6;state_j++){
        //                        if(state_j>state_i){
        //                            MFParams_.OParams[site][state_i][state_j]=conj(MFParams_.OParams[site][state_j][state_i]);
        //                        }
        //                        else{
        //                            assert(state_j==state_i);
        //                            MFParams_.OParams[site][state_i][state_j].imag(0.0);
        //                        }
        //                    }}}}


    }



}

void Observables::Update_OrderParameters_Modified_Broyden_in_ComplexSpace(int iter_count){

    /*
    //For details see your own notes at
    //"https://github.com/nkphys/3_ORB_SOC_Hartree_Fock/tree/master/Notes/Modified_Broyden_Method"

    //XXXXXXXXXXLiteratureXXXXXXXXXXXXX
    //Modified Broyden is used from "D. D. Johnson, Phys. Rev. B 38, 12807, 1988".
    //Look into this "https://arxiv.org/pdf/0805.4446.pdf" as well.
    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    double alpha_OP=Parameters_.alpha_OP;
    double normalization_;
    int site;
    int iter;


    iter=iter_count%Parameters_.ModBroydenCounter;


    if(iter==0){
        //****Getting ready for iters>0*********
        _Delta_F.clear();
        _u.clear();
        _A.clear();
        _cm.clear();
        _gammam.clear();
        w_minus1=Parameters_.w_minus1;
        w.clear();
        //************************************

        cout<<"Using Modified Broyden Mixing to gain Self-Consistency"<<endl;

        //Get Fm
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);

                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        _Fm[site*36 + state_i*6 + state_j]=OParams_[site][state_i][state_j]
                                - MFParams_.OParams[site][state_i][state_j];
                    }
                }
            }
        }


        for(int j=0;j<6*6*ns_;j++){
            _Delta_OPm[j] = alpha_OP*_Fm[j];
        }



        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);
                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        MFParams_.OParams[site][state_i][state_j] +=_Delta_OPm[site*36 + state_i*6 + state_j];

                    }
                }
            }
        }

        //Copy Jinv_np1 to Jinv_n
        _Delta_OPm_minus1 = _Delta_OPm;

        //Copy F_n to F_nm1
        _Fm_minus1=_Fm;

    }

    else{

        w.resize(iter);
        w[iter-1]=Parameters_.wn;

        //Get Fm******************
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);

                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        _Fm[site*36 + state_i*6 + state_j]=OParams_[site][state_i][state_j]
                                - MFParams_.OParams[site][state_i][state_j];
                    }
                }
            }
        }
        //******************************


        //Get DeltaFm/|DeltaFm|-------------------------//
        _Delta_F.resize(iter);
        _Delta_F[iter-1].resize(36*ns_);

        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);

                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        _Delta_F[iter-1][site*36 + state_i*6 + state_j]=_Fm[site*36 + state_i*6 + state_j]
                                - _Fm_minus1[site*36 + state_i*6 + state_j];
                    }
                }
            }
        }

        //Getting sqrt(<DeltaFm|DeltaFm>)
        normalization_=0.0;
        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);
                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        normalization_ += (conj(_Delta_F[iter-1][site*36 + state_i*6 + state_j])*
                                _Delta_F[iter-1][site*36 + state_i*6 + state_j]).real();
                    }
                }
            }
        }
        normalization_=sqrt(normalization_);

        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);
                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        _Delta_F[iter-1][site*36 + state_i*6 + state_j]=
                                _Delta_F[iter-1][site*36 + state_i*6 + state_j]*(1.0/normalization_);
                    }
                }
            }
        }
        //--------------------------------------------//


        //Getting Delta_n/|DeltaFm|-------------------//
        for(int i=0;i<36*ns_;i++){
            _Delta_OPm_minus1[i]=_Delta_OPm_minus1[i]*(1.0/normalization_);
        }
        //-----------------------------------------------//


        //Get u[iter-1]------------------------------//
        _u.resize(iter);
        _u[iter-1].resize(36*ns_);
        for(int i=0;i<36*ns_;i++){
            _u[iter-1][i] = (alpha_OP*_Delta_F[iter-1][i])  +  _Delta_OPm_minus1[i];
        }
        //-------------------------------------------------//


        //UPDATE _A----------------------------//

        Mat_1_Complex_doub temp_vec;
        temp_vec.resize(iter);
        _A.push_back(temp_vec);
        temp_vec.clear();

        for(int i=0;i<_A.size();i++){
            _A[i].resize(iter);
        }


        for(int i=0;i<iter;i++){
            if(i==(iter-1)){
                _A[i][i]=w[i]*w[i]*DOT_P(_Delta_F[i],_Delta_F[i]);
            }
            else{
                _A[iter-1][i]=w[iter-1]*w[i]*DOT_P(_Delta_F[i],_Delta_F[iter-1]);
                _A[i][iter-1]=w[iter-1]*w[i]*DOT_P(_Delta_F[iter-1],_Delta_F[i]);
            }
        }
        //---------------------------------------------//

        //Get Beta--------------------------------------//
        _Beta.resize(iter,iter);
        for(int i=0;i<iter;i++){
            for(int j=0;j<iter;j++){
                if(i==j){
                    _Beta(i,j) = (one_complex*w_minus1*w_minus1) + _A[i][j];
                }
                else{
                    _Beta(i,j) = _A[i][j];
                }
            }

        }

        //  cout<<"Before Inversion:"<<endl;
        //Matrix<complex<double>> Beta_before;
        //Beta_before=_Beta;
        //Beta_before.print();

        Invert_Beta();
        for(int i=0;i<_Beta.n_col();i++){
            for(int j=0;j<i;j++){
                _Beta(i,j)=_Beta(j,i);
            }
        }


        //cout<<"After Inversion:"<<endl;
        //Matrix<complex<double>> Beta_after;
        //Beta_after=_Beta;
        //Beta_after.print();

        //cout<<"Identity_check:"<<endl;
        //Matrix<complex<double>> Identity_check;
        //Identity_check=product(Beta_before, Beta_after);
        //Identity_check.print();


        //-----------------------------------------------//

        //Get _cm-------------------------------------------//
        _cm.clear();
        _cm.resize(iter);
        for(int i=0;i<iter;i++){
            _cm[i]=w[i]*DOT_P(_Delta_F[i],_Fm);
        }
        //---------------------------------------------------//

        //Get _gammam------------------------------------------//
        _gammam.clear();
        _gammam.resize(iter);
        for(int l=0;l<iter;l++){
            _gammam[l]=zero_complex;
            for(int k=0;k<iter;k++){
                _gammam[l] += _cm[k]*_Beta(k,l);
            }
        }
        //--------------------------------------------------//


        //Get _Delta_OPm-----------------------------------------//
        for(int i=0;i<36*ns_;i++){
            _Delta_OPm[i]=zero_complex;
            for(int n=0;n<iter;n++){
                _Delta_OPm[i] += (-1.0*one_complex)*w[n]*_gammam[n]*_u[iter-1][i];
            }
        }

        for(int i=0;i<36*ns_;i++){
            _Delta_OPm[i] += alpha_OP*_Fm[i];
        }
        //---------------------------------------------//




        for(int i=0;i<lx_;i++){
            for(int j=0;j<ly_;j++){
                site=Coordinates_.Nc(i,j);
                for(int state_i=0;state_i<6;state_i++){
                    for(int state_j=0;state_j<6;state_j++){
                        MFParams_.OParams[site][state_i][state_j] +=_Delta_OPm[site*36 + state_i*6 + state_j];

                    }
                }
            }
        }

        //Copy Jinv_np1 to Jinv_n
        _Delta_OPm_minus1 = _Delta_OPm;

        //Copy F_n to F_nm1
        _Fm_minus1=_Fm;

    }

    //Maintain hermiticity of order params, because after many iterations they tend to loose it
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);
            for(int state_i=0;state_i<6;state_i++){
                for(int state_j=state_i;state_j<6;state_j++){
                    if(state_j>state_i){
                        MFParams_.OParams[site][state_i][state_j]=conj(MFParams_.OParams[site][state_j][state_i]);
                    }
                    else{
                        assert(state_j==state_i);
                        MFParams_.OParams[site][state_i][state_j].imag(0.0);
                    }
                }}}}


*/
}


void Observables::Calculate_IPR(){
    /*
    double IPR;
    int c1, site;
    double eta = 0.001;
    int n_chosen=(Parameters_.ns*Parameters_.Fill*2.0) - 1;
    IPR=0.0;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;
                //  IPR += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                //         abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                //         Lorentzian( Parameters_.mus - Hamiltonian_.eigs_[n], eta);

                IPR += abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen))*
                        abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen));




            }
        }
    }

    cout<<"IPR for state no. "<<n_chosen<<" = "<<IPR<<endl;



    IPR=0.0;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){


                    IPR += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            Lorentzian( Parameters_.mus - Hamiltonian_.eigs_[n], eta);

                }
            }
        }
    }

    cout<<"IPR for near mu(with eta = "<<eta<<") = "<<IPR<<endl;




    string fileout_FermiState="Fermi_state_probability.txt";
    ofstream file_FermiState_out(fileout_FermiState.c_str());
    file_FermiState_out<<"#ix   iy   site   |Psi_{Fermi,up}(ix,iy)|^2 + |Psi_{Fermi,dn}(ix,iy)|^2"<<endl;

    double value;
    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);
            value = 0.0;

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;

                value += abs(Hamiltonian_.Ham_(c1,n_chosen))*abs(Hamiltonian_.Ham_(c1,n_chosen));
            }

            file_FermiState_out<<i<<"\t"<<j<<"\t"<<site<<"\t"<<value<<endl;

        }
        file_FermiState_out<<endl;
    }





    string fileout_Near_mu="Near_mu_probability.txt";
    ofstream file_Near_mu_out(fileout_Near_mu.c_str());
    file_Near_mu_out<<"#ix   iy   site   sum_{n}Lorentz(near_mu)*|Psi_{n,up}(ix,iy)|^2 + |Psi_{n,dn}(ix,iy)|^2"<<endl;


    for(int i=0;i<lx_;i++){
        for(int j=0;j<ly_;j++){
            site=Coordinates_.Nc(i,j);
            value = 0.0;

            for(int spin=0;spin<2;spin++){
                c1=site + spin*ns_;
                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    value += abs(Hamiltonian_.Ham_(c1,n))*abs(Hamiltonian_.Ham_(c1,n))*
                            Lorentzian(Parameters_.mus - Hamiltonian_.eigs_[n], eta);;
                }
            }

            file_Near_mu_out<<i<<"\t"<<j<<"\t"<<site<<"\t"<<value<<endl;

        }
        file_Near_mu_out<<endl;
    }



*/
}


void Observables::Calculate_Local_n_orb_resolved(){

    string fileout_local_den="Local_orb_spin_resolved_densities.txt";
    ofstream file_local_den_out(fileout_local_den.c_str());

    Local_n_orb_resolved.clear();
    Local_n_orb_resolved.resize(ns_*2*n_orbs_);

    int c1;
    int site;


    file_local_den_out<<"# site        ix         iy       n[orb=0,spin=0]      n[orb=0,spin=1]     n[orb=1,spin=0]     n[orb=1,spin=1] ........"<<endl;
    for(int iy=0;iy<ly_;iy++){
        for(int ix=0;ix<lx_;ix++){
            site = Coordinates_.Nc(ix,iy);
            file_local_den_out << site << "     " <<iy<<"     " <<ix<<"     ";


            for(int orb=0;orb<n_orbs_;orb++){
                for(int spin=0;spin<2;spin++){

                    c1=Coordinates_.Nc_dof(site,orb+n_orbs_*spin);
                    Local_n_orb_resolved[c1]=0;

                    for(int n=0;n<Hamiltonian_.eigs_.size();n++){
                        Local_n_orb_resolved[c1] += (
                                    ( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n))*
                                    (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                                    ).real();
                    }

                    file_local_den_out<<Local_n_orb_resolved[c1]<<"     ";

                }
            }
            file_local_den_out<<endl;
        }
        file_local_den_out<<endl;

    }

}

void Observables::Calculate_Order_Params(){


    int c1,c2;

    Local_Sz_.resize(ns_);
    Local_Sx_.resize(ns_);
    Local_Sy_.resize(ns_);

    Local_Lz_.resize(ns_);
    Local_Lx_.resize(ns_);
    Local_Ly_.resize(ns_);


    for(int site=0;site<ns_;site++){

        Local_Sz_[site]=0.0;
        Local_Sx_[site]=0.0;
        Local_Sy_[site]=0.0;
        Local_Lz_[site]=0.0;
        Local_Lx_[site]=0.0;
        Local_Ly_[site]=0.0;


        for(int n=0;n<Hamiltonian_.eigs_.size();n++){

            //Sz
            for(int orb=0;orb<n_orbs_;orb++){
                c1=Coordinates_.Nc_dof(site,orb+n_orbs_*0); //up
                c2=Coordinates_.Nc_dof(site,orb+n_orbs_*1); //down
                Local_Sz_[site] += (
                            ( 0.5*( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)
                                    - conj(Hamiltonian_.Ham_(c2,n))*Hamiltonian_.Ham_(c2,n)) )*
                            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                            ).real();
            }


            //Sx
            for(int orb=0;orb<n_orbs_;orb++){
                c1=Coordinates_.Nc_dof(site,orb+n_orbs_*0); //up
                c2=Coordinates_.Nc_dof(site,orb+n_orbs_*1); //down
                Local_Sx_[site] += (
                            ( 0.5*( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)
                                    + conj(Hamiltonian_.Ham_(c2,n))*Hamiltonian_.Ham_(c1,n)) )*
                            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                            ).real();
            }

            //Sy
            for(int orb=0;orb<n_orbs_;orb++){
                c1=Coordinates_.Nc_dof(site,orb+n_orbs_*0); //up
                c2=Coordinates_.Nc_dof(site,orb+n_orbs_*1); //down
                Local_Sy_[site] += (
                            (0.5*(-1.0*iota_complex)*( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)
                                                       - conj(Hamiltonian_.Ham_(c2,n))*Hamiltonian_.Ham_(c1,n)) )*
                            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                            ).real();
            }


            //Lz
            for(int spin_=0;spin_<2;spin_++){
                c1=Coordinates_.Nc_dof(site, 0 + n_orbs_*spin_);
                c2=Coordinates_.Nc_dof(site, 1 + n_orbs_*spin_);
                Local_Lz_[site] += (
                            ( 0.5*( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)
                                    - conj(Hamiltonian_.Ham_(c2,n))*Hamiltonian_.Ham_(c2,n)) )*
                            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                            ).real();
            }

            //Lx
            for(int spin_=0;spin_<2;spin_++){
                c1=Coordinates_.Nc_dof(site, 0 + n_orbs_*spin_);
                c2=Coordinates_.Nc_dof(site, 1 + n_orbs_*spin_);
                Local_Lx_[site] += (
                            ( 0.5*( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)
                                    + conj(Hamiltonian_.Ham_(c2,n))*Hamiltonian_.Ham_(c1,n)) )*
                            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                            ).real();
            }

            //Ly
            for(int spin_=0;spin_<2;spin_++){
                c1=Coordinates_.Nc_dof(site, 0 + n_orbs_*spin_);
                c2=Coordinates_.Nc_dof(site, 1 + n_orbs_*spin_);
                Local_Ly_[site] +=  (
                            (0.5*(-1.0*iota_complex)*( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)
                                                       - conj(Hamiltonian_.Ham_(c2,n))*Hamiltonian_.Ham_(c1,n)) )*
                            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                            ).real();
            }


        }//n

    }


    //checking Hermiticity of order parameters


}


void Observables::Calculate_Akw(){


    //---------Read from input file-----------------------//
    string fileout="Akw_orbital_resolved.txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.08;
    omega_min=Hamiltonian_.eigs_[0]-0.5-Parameters_.mus;omega_max=Hamiltonian_.eigs_[2*n_orbs_*ns_ -1]+0.5-Parameters_.mus;d_omega=0.0005;
    //---------------------------------------------------//


    int UP_=0;
    int DN_=1;

    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );

    ofstream file_Akw_out(fileout.c_str());
    file_Akw_out<<"# k_point   kx_i   ky_i    (ky_i*Parameters_.lx) + kx_i    omega_min + (d_omega*omega_ind)   omega_ind    Akw_orb0    Akw_orb1   ......"<<endl;

    int c1,c2;

    Mat_4_Complex_doub A_orb;
    A_orb.resize(n_orbs_);

    for(int orb=0;orb<n_orbs_;orb++){
        A_orb[orb].resize(Parameters_.ns);

        for (int i=0;i<Parameters_.ns;i++){
            A_orb[orb][i].resize(Parameters_.ns);

            for(int j=0;j<Parameters_.ns;j++){
                A_orb[orb][i][j].resize(omega_index_max);
            }
        }
    }


    //   complex<double> Nup_check(0,0);
    //   complex<double> Ndn_check(0,0);

    for (int j=0;j<Parameters_.ns;j++){
        cout<<"Akw for j="<<j<<" done"<<endl;
        for (int l=0;l<Parameters_.ns;l++){

            for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){

                for(int orb=0;orb<n_orbs_;orb++){

                    A_orb[orb][j][l][omega_ind]=zero_complex;

                    for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){

                        //c= l + or1*ns_ + ns_*orbs_*spin;
                        for(int spin=0;spin<2;spin++){
                            //Hamiltonian_.Ham_(c2,n) is nth eigenvector and c2th component [checked];
                            c1 = Coordinates_.Nc_dof(l,orb+ n_orbs_*spin);
                            c2 = Coordinates_.Nc_dof(j,orb + n_orbs_*spin);
                            A_orb[orb][j][l][omega_ind] +=  conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c2,n)*
                                    Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta);

                        }
                    }
                }
            }
        }
    }

    // cout << "Nup_check = "<<Nup_check<<endl;
    //  cout << "Ndn_check = "<<Ndn_check<<endl;

    complex<double> temp_orb;

    double kx,ky;
    int kx_i,ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;

    //--------For 1D in x direction-----
    /*
    ky_i=0;
    for(kx_i=-1;kx_i<=(Parameters_.lx);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    */
    //----------------------------------



    //--------\Gamma to X-----------------
    ky_i=0;
    for(kx_i=0;kx_i<=(Parameters_.lx/2);kx_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i=(Parameters_.lx/2);
    for(ky_i=1;ky_i<=(Parameters_.lx/2);ky_i++){
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------M to \Gamma[with one extra point,
    //                  because in gnuplor use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i=(Parameters_.lx/2) - 1;
    ky_i=(Parameters_.lx/2) - 1;
    for(kx_i=(Parameters_.lx/2) - 1;kx_i>=-1;kx_i--){
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    for(int k_point=0;k_point<k_path.size();k_point++){

        kx_i=k_path[k_point].first;
        ky_i=k_path[k_point].second;
        kx=(2.0*PI*kx_i)/(1.0*Parameters_.lx);
        ky=(2.0*PI*ky_i)/(1.0*Parameters_.ly);

        for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){

            file_Akw_out<< k_point<<"   "<<kx_i<<"   "<<ky_i<<"   "<<(ky_i*Parameters_.lx) + kx_i<<"    "<<
                           omega_min + (d_omega*omega_ind)<<"   "<<omega_ind<<"    ";

            for(int orb=0;orb<n_orbs_;orb++){
                temp_orb=zero_complex;

                for(int j=0;j<ns_;j++){
                    for(int l=0;l<ns_;l++){
                        temp_orb += one_complex*
                                exp(iota_complex*(kx*(Coordinates_.indx(j) - Coordinates_.indx(l)) +
                                                  ky*(Coordinates_.indy(j) - Coordinates_.indy(l))))*
                                A_orb[orb][j][l][omega_ind];
                    }
                }

                file_Akw_out<<temp_orb.real()<<"     "<<temp_orb.imag()<<"     ";
            }
            file_Akw_out<<endl;
        }
        file_Akw_out<<endl;
    }



}


void Observables::Calculate_Nw(){

    //---------Read from input file-----------------------//
    double omega_min, omega_max, d_omega;
    double eta = 0.05;
    omega_min=Hamiltonian_.eigs_[0]-0.5-Parameters_.mus;omega_max=Hamiltonian_.eigs_[2*n_orbs_*ns_ - 1]+0.5-Parameters_.mus;d_omega=0.005;
    //---------------------------------------------------//

    int c1;
    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );
    double temp_val ;


    //---------------------------------------------------------------------------------//
    //************************************Nw_jm****************************************//
    //---------------------------------------------------------------------------------//

    int c2;

    string fileout_orb="Nw_orbital_spin_resolved.txt";
    ofstream file_Nw_out_orb(fileout_orb.c_str());

    file_Nw_out_orb<<"#(w-mu)    orb0_up    orb1_up    .....     ";
    file_Nw_out_orb<<"orb0_dn    orb1_dn  ........."<<endl;

    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
        file_Nw_out_orb<<omega_min + (omega_ind*d_omega)<<"       ";


        for(int dof_type=0;dof_type<2*n_orbs_;dof_type++){
            temp_val=0.0;

            for(int site=0;site<Coordinates_.ns_;site++){
                c1=Coordinates_.Nc_dof_(site,dof_type);

                for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
                    temp_val += ( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n)*
                                  Lorentzian( omega_min + (omega_ind*d_omega) - (Hamiltonian_.eigs_[n] -Parameters_.mus), eta)).real();

                }

            }
            file_Nw_out_orb<<temp_val<<"      ";
        }
        file_Nw_out_orb<<endl;
    }

    file_Nw_out_orb<<"#actual mu = "<< Parameters_.mus<<", but shifted to 0"<<endl;

    //---------------------------------------------------------------------------------//
    //*********************************************************************************//






    int n_chosen=Parameters_.Total_Particles - 1;
    cout<<"Gap = "<<Hamiltonian_.eigs_[n_chosen+1] - Hamiltonian_.eigs_[n_chosen]<<endl;

    string fileout_Eigen="Eigen_spectrum.txt";
    ofstream file_Eigen_out(fileout_Eigen.c_str());

    file_Eigen_out<<"# n(eigen-index)   E(n)    Fermi(n)"<<endl;
    for(int n=0;n<Hamiltonian_.Ham_.n_row();n++){
        file_Eigen_out<<n<<"\t"<<Hamiltonian_.eigs_[n]<<"\t"<<(1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)) <<endl;
    }


}

void Observables::Calculate_Optical_Conductivity(){
    /*

    string fileout_sigma_w = "Optical_conductivity.txt";
    ofstream file_sigma_w_out(fileout_sigma_w.c_str());
    file_sigma_w_out<<"#omega   Sigma_xx  Sigma_yy"<<endl;

    //--------------------------------------------------//
    double omega_min, omega_max, d_omega;
    double eta = 0.05;
    omega_min=0.00001;omega_max=10.0;d_omega=0.001;
    //---------------------------------------------------//


    Mat_2_doub PSI_x, PSI_y;
    complex<double> value_x, value_y;
    int ipx, ipy;




    PSI_x.resize(2*ns_);
    PSI_y.resize(2*ns_);
    for(int n=0;n<2*ns_;n++){
        PSI_x[n].resize(2*ns_);
        PSI_y[n].resize(2*ns_);
    }



    for(int n=0;n<2*ns_;n++){
        for(int m=0;m<2*ns_;m++){

            value_x=zero_complex;
            value_y=zero_complex;

            for(int i=0;i<ns_;i++){
                ipx = Coordinates_.neigh(i,0);
                ipy = Coordinates_.neigh(i,2);



                for(int spin=0;spin<2;spin++){
                    value_x += ( conj(Hamiltonian_.Ham_(ipx + (ns_*spin),n))*Hamiltonian_.Ham_(i + (ns_*spin),m) )
                            -  ( conj(Hamiltonian_.Ham_(i + (ns_*spin),n))*Hamiltonian_.Ham_(ipx + (ns_*spin),m) );

                    value_y += ( conj(Hamiltonian_.Ham_(ipy + (ns_*spin),n))*Hamiltonian_.Ham_(i + (ns_*spin),m) )
                            -  ( conj(Hamiltonian_.Ham_(i + (ns_*spin),n))*Hamiltonian_.Ham_(ipy + (ns_*spin),m) );

                }


            }


            PSI_x[n][m] = abs(value_x)*abs(value_x);
            PSI_y[n][m] = abs(value_y)*abs(value_y);


        }
    }




    double sigma_x, sigma_y;
    double omega_val;


    int omega_index_max = int( (omega_max - omega_min)/(d_omega) );


    //cout<<"omega_index_max = "<<omega_index_max<<endl;
    for(int omega_ind=0;omega_ind<omega_index_max;omega_ind++){
        omega_val = omega_min + (omega_ind*d_omega);

        //cout<<omega_ind<<endl;

        sigma_x=0.0; sigma_y=0.0;
        for(int n=0;n<2*ns_;n++){
            for(int m=0;m<2*ns_;m++){

                if(n!=m){
                    sigma_x += (PSI_x[n][m])
                            *((1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)))
                            *((1.0/( exp((Parameters_.mus-Hamiltonian_.eigs_[m])*Parameters_.beta ) + 1.0)))
                            *Lorentzian( omega_min + (omega_ind*d_omega) + Hamiltonian_.eigs_[n] - Hamiltonian_.eigs_[m], eta);

                    sigma_y += (PSI_y[n][m])
                            *((1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)))
                            *((1.0/( exp((Parameters_.mus-Hamiltonian_.eigs_[m])*Parameters_.beta ) + 1.0)))
                            *Lorentzian( omega_min + (omega_ind*d_omega) + Hamiltonian_.eigs_[n] - Hamiltonian_.eigs_[m], eta);


                }
            }
        }

        sigma_x = sigma_x*PI*(1.0 - exp(-1.0*Parameters_.beta*(omega_val)))*(1.0/(omega_val*ns_));
        sigma_y = sigma_y*PI*(1.0 - exp(-1.0*Parameters_.beta*(omega_val)))*(1.0/(omega_val*ns_));

        file_sigma_w_out<<omega_val<<"     "<<sigma_x<<"     "<<sigma_y<<endl;

    }



*/
}

void Observables::Calculate_Single_Particle_Density_Matrix(){

    /*
      NOTE:
      SP_Density_Matrix[alpha][beta] = <c_{alpha^{daggger}} c_{beta}>
     */
    SP_Density_Matrix.resize(ns_*2*n_orbs_);
    for(int i=0;i<ns_*2*n_orbs_;i++){
        SP_Density_Matrix[i].resize(ns_*2*n_orbs_);
    }

    for(int alpha_=0;alpha_<ns_*2*n_orbs_;alpha_++){
        for(int beta_=0;beta_<ns_*2*n_orbs_;beta_++){
            SP_Density_Matrix[alpha_][beta_] = zero_complex;
            for(int n=0;n<ns_*2*n_orbs_;n++){
                SP_Density_Matrix[alpha_][beta_] += conj(Hamiltonian_.Ham_(alpha_,n))*Hamiltonian_.Ham_(beta_,n)*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));
            }
        }
    }

}

void Observables::Calculate_two_point_correlations(){


    int UP_, DOWN_;
    UP_=0;DOWN_=1;

    int _Sz, _Sx,_Sy, _Lz, _Lx,_Ly;

    Mat_1_string opr_type;
    opr_type.clear();
    opr_type.push_back("Sz");_Sz=0;
    opr_type.push_back("Sx");_Sx=1;
    opr_type.push_back("Sy");_Sy=2;
    opr_type.push_back("Lz");_Lz=3;
    opr_type.push_back("Lx");_Lx=4;
    opr_type.push_back("Ly");_Ly=5;


    Mat_3_Complex_doub Oprs_;
    Oprs_.resize(opr_type.size());
    for(int opr_no_=0;opr_no_<opr_type.size();opr_no_++){
        Oprs_[opr_no_].clear();
        Oprs_[opr_no_].resize(2*n_orbs_);
        for(int i=0;i<2*n_orbs_;i++){
            Oprs_[opr_no_][i].resize(2*n_orbs_);
        }
    }




    //Opr. =\sum_{orb1,orb2,spin1,spin2}Temp_Mat{orb1,spin1;orb2,spin2}c*_{orb1,spin1}c_{orb2,spin2}
    //    for(int orb1=0;orb1<3;orb1++){
    //        for(int spin1=0;spin1<2;spin1++){
    //            //orb1 + 3*spin1
    //            for(int orb2=0;orb2<3;orb2++){
    //                for(int spin2=0;spin2<2;spin2++){
    //                }
    //            }
    //        }
    //    }

    int opr_no;

    //Sz---
    assert(opr_type[_Sz]=="Sz");
    for(int orb1=0;orb1<n_orbs_;orb1++){
        for(int spin1=0;spin1<2;spin1++){
            Oprs_[_Sz][orb1 + n_orbs_*spin1][orb1 + n_orbs_*spin1]=one_complex*(0.5*(1.0-(2.0*spin1)));
        }}


    //Sx---
    assert(opr_type[_Sx]=="Sx");
    for(int orb1=0;orb1<n_orbs_;orb1++){
        Oprs_[_Sx][orb1 + n_orbs_*0][orb1 + n_orbs_*1]=one_complex*(0.5);
        Oprs_[_Sx][orb1 + n_orbs_*1][orb1 + n_orbs_*0]=one_complex*(0.5);
    }

    //Sy---
    assert(opr_type[_Sy]=="Sy");
    for(int orb1=0;orb1<n_orbs_;orb1++){
        Oprs_[_Sy][orb1 + n_orbs_*0][orb1 + n_orbs_*1]=(-1.0*iota_complex)*(0.5);
        Oprs_[_Sy][orb1 + n_orbs_*1][orb1 + n_orbs_*0]=(-1.0*iota_complex)*(-0.5);
    }

    //Lz---
    assert(opr_type[_Lz]=="Lz");
    for(int spin1=0;spin1<2;spin1++){
        Oprs_[_Lz][0 + n_orbs_*spin1][0 + n_orbs_*spin1]=one_complex*(0.5);
        Oprs_[_Lz][1 + n_orbs_*spin1][1 + n_orbs_*spin1]=one_complex*(-0.5);
    }

    //Lx---
    opr_no=_Lx;
    assert(opr_type[opr_no]=="Lx");
    for(int spin1=0;spin1<2;spin1++){
        Oprs_[_Lx][0 + n_orbs_*spin1][1 + n_orbs_*spin1]=one_complex*(0.5);
        Oprs_[_Lx][1 + n_orbs_*spin1][0 + n_orbs_*spin1]=one_complex*(0.5);
    }

    //Ly---
    opr_no=_Ly;
    assert(opr_type[opr_no]=="Ly");
    for(int spin1=0;spin1<2;spin1++){
        Oprs_[_Ly][0 + n_orbs_*spin1][1 + n_orbs_*spin1]=(-1.0*iota_complex)*(0.5);
        Oprs_[_Ly][1 + n_orbs_*spin1][0 + n_orbs_*spin1]=(-1.0*iota_complex)*(-0.5);
    }



    string corrs_out = "corrs.txt";
    ofstream file_corrs_out(corrs_out.c_str());
    file_corrs_out<<"#site_i   site_i(x)    site_i(y)    site_j   site_j(x)    site_j(y)     SS[site_i][site_j]     LL[site_i][site_j]"<<endl;

    int i1,i2,j1,j2;

    complex<double> temp_val_SS, temp_val_LL;
    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            temp_val_SS = zero_complex;
            temp_val_LL = zero_complex;

            for(int i_row=0;i_row<2*n_orbs_;i_row++){
                for(int i_col=0;i_col<2*n_orbs_;i_col++){

                    for(int j_row=0;j_row<2*n_orbs_;j_row++){
                        for(int j_col=0;j_col<2*n_orbs_;j_col++){


                            //S.S
                            for(int comp=_Sz;comp<_Sy+1;comp++){
                                if( (Oprs_[comp][i_row][i_col] !=zero_complex)
                                        &&
                                        (Oprs_[comp][j_row][j_col] !=zero_complex)   ){
                                    i1=Coordinates_.Nc_dof(i,i_row);
                                    i2=Coordinates_.Nc_dof(i,i_col);
                                    j1=Coordinates_.Nc_dof(j,j_row);
                                    j2=Coordinates_.Nc_dof(j,j_col);
                                    temp_val_SS += Oprs_[comp][i_row][i_col]*Oprs_[comp][j_row][j_col]*
                                            Two_particle_Den_Mat(i1,i2,j1,j2);
                                }
                            }

                            //L.L
                            for(int comp=_Lz;comp<_Ly+1;comp++){
                                if( (Oprs_[comp][i_row][i_col] !=zero_complex)
                                        &&
                                        (Oprs_[comp][j_row][j_col] !=zero_complex)   ){
                                    i1=Coordinates_.Nc_dof(i,i_row);
                                    i2=Coordinates_.Nc_dof(i,i_col);
                                    j1=Coordinates_.Nc_dof(j,j_row);
                                    j2=Coordinates_.Nc_dof(j,j_col);
                                    temp_val_LL += Oprs_[comp][i_row][i_col]*Oprs_[comp][j_row][j_col]*
                                            Two_particle_Den_Mat(i1,i2,j1,j2);
                                }
                            }

                        }
                    }
                }
            }
            file_corrs_out<<i<<setw(15)<<Coordinates_.indx(i)<<setw(15)<<Coordinates_.indy(i)<<setw(15)<<j<<setw(15)<<Coordinates_.indx(j)<<setw(15)<<Coordinates_.indy(j)<<
                            setw(15)<<temp_val_SS.real()<<//"\t"<<temp_val_SS.imag()<<
                            setw(15)<<temp_val_LL.real()<<//"\t"<<temp_val_LL.imag()<<
                            endl;

        }
    }



    //------------------------------//
}

complex<double> Observables::Two_particle_Den_Mat(int _alpha, int _beta, int _gamma, int _delta){

    complex<double> temp;
    complex<double> delta_gamma_beta;

    if(_gamma == _beta){
        delta_gamma_beta=one_complex;
    }
    else{
        assert(_gamma != _beta);
        delta_gamma_beta=zero_complex;
    }

    temp = (SP_Density_Matrix[_alpha][_beta]*SP_Density_Matrix[_gamma][_delta])
            +
            (SP_Density_Matrix[_alpha][_delta]*(delta_gamma_beta - SP_Density_Matrix[_gamma][_beta]));
    return temp;
}





void Observables::Get_Non_Interacting_dispersion(){

}


double Observables::Lorentzian(double x, double brd){
    double temp;

    temp = (1.0/PI)*( (brd/2.0)/ ( (x*x) + ((brd*brd)/4.0) ) );

    return temp;

}




void Observables::Total_Energy_Average(double Curr_QuantE, double CurrE){

    AVG_Total_Energy += Curr_QuantE + CurrE;
    AVG_Total_Energy_sqr += (Curr_QuantE + CurrE)*(Curr_QuantE + CurrE);
}





void Observables::Initialize(){


    F_n.resize(ns_*6); //F_n=x_n_out - x_n_in []
    F_nm1.resize(ns_*6);
    DeltaF_n.resize(ns_*6); //DeltaF_n=F_n - F-nm1;
    Delta_x_n.resize(ns_*6); //Delta_x_n= x_n_in - x_nm1_in;
    Jinv_n.resize(ns_*6);
    Jinv_np1.resize(ns_*6);


    _Fm.resize((ns_*6));
    _Delta_OPm.resize((ns_*6));
    _Fm_minus1.resize((ns_*6));
    _Delta_OPm_minus1.resize((ns_*6));



    for(int i=0;i<ns_*6;i++){
        Jinv_n[i]. resize(ns_*6);
        Jinv_np1[i]. resize(ns_*6);
    }



} // ----------












#endif // OBSERVABLES_H
