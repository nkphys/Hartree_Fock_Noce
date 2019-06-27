#ifndef Parameters_class
#define Parameters_class
#include "tensor_type.h"

class Parameters{

public:
    int lx, ly, ns, IterMax, RandomSeed;
    int n_orbs;
    double Convergence_Error;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    int lx_cluster, ly_cluster;
    double mus,Total_Particles,pi, mu_old;
    double J_Hund;
    double U_onsite;
    double U_prime_onsite;
    double Disorder_Strength, RandomDisorderSeed;
    bool PBC;

    bool Read_OPs;
    string File_OPs_in, File_OPs_out;

    Matrix<double> Hopping_NN;
    Mat_1_doub Crystal_Field;

    bool Simple_Mixing;
    bool Broyden_Mixing;
    bool BroydenSecondMethodMixing;
    double w_minus1,wn;
    int BroydenSecondMethodCounter;
    double alpha_OP;

    double Temperature,beta,Eav,maxmoment;

    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file,string match);
    string matchstring2(string file,string match);

};


void Parameters::Initialize(string inputfile_){



    double Simple_Mixing_double, Broyden_Mixing_double, BroydenSecondMethodMixing_double;
    double Read_OPs_double;
    string PBC_string;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: " << inputfile_ << endl;
    cout << "____________________________________" << endl;



    lx = int(matchstring(inputfile_,"Xsite"));
    ly = int(matchstring(inputfile_,"Ysite"));
    n_orbs = int(matchstring(inputfile_,"n_orbs"));
    if(n_orbs !=2){
        cout<<"Firstly figure out Eq-9 of Phys. Status Solidi B 251, No. 4, 907–911 (2014)"<<endl;
        cout<<"Right now this only works for n_orbs=2"<<endl;
        cout<<"For higher n_orbs add hiigher dimensional pauli matrices"<<endl;
    }
    assert(n_orbs==2);

    TBC_mx = int(matchstring(inputfile_,"TwistedBoundaryCond_mx"));
    TBC_my = int(matchstring(inputfile_,"TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_,"TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_,"TBC_cellsY"));
    BroydenSecondMethodCounter = int(matchstring(inputfile_,"BroydenSecondMethodCounter"));

    ns = lx*ly;
    cout << "TotalNumberOfSites = "<< ns << endl;

    Total_Particles = matchstring(inputfile_,"Total No. of particles");
    cout << "TotalNumberOfParticles = "<< Total_Particles << endl;

    IterMax = int(matchstring(inputfile_,"No_of_SelfConsistency_iters"));
    Convergence_Error=matchstring(inputfile_,"Convergence_Error");
    RandomSeed = matchstring(inputfile_,"RandomSeed");
    RandomDisorderSeed = matchstring(inputfile_,"RandomDisorderSeed");
    Disorder_Strength = matchstring(inputfile_,"Disorder_Strength");
    J_Hund = matchstring(inputfile_,"J_HUND");
    U_onsite = matchstring(inputfile_,"U_Onsite");
    U_prime_onsite = matchstring(inputfile_,"U_prime_Onsite");
    if( U_prime_onsite != (U_onsite - (2.0*J_Hund)) ){
        cout << " U_prime_onsite = (U_onsite - (2.0*J_Hund)) have to be used, ";
        cout << "because this code only works for rotational invariant cases"<<endl;
    }
    assert(U_prime_onsite == (U_onsite - (2.0*J_Hund)));

    alpha_OP = matchstring(inputfile_,"alpha_OP");
    w_minus1 = matchstring(inputfile_,"w_minus1");
    wn = matchstring(inputfile_,"wn");



    Dflag = 'N';

    Simple_Mixing_double=double(matchstring(inputfile_,"Simple_Mixing"));
    Broyden_Mixing_double=double(matchstring(inputfile_,"Broyden_Mixing"));
    BroydenSecondMethodMixing_double=double(matchstring(inputfile_,"Broyden_Second_Method_Mixing"));


    if(BroydenSecondMethodMixing_double==1.0){
         BroydenSecondMethodMixing=true;
        Broyden_Mixing=false;
        Simple_Mixing=false;
    }
    else{
         BroydenSecondMethodMixing=false;
        if(Broyden_Mixing_double==1.0){
            Broyden_Mixing=true;
            Simple_Mixing=false;

        }
        else if(Broyden_Mixing_double==0.0){
            Broyden_Mixing=false;
            Simple_Mixing=true;
            cout<<"Broyden_Mixing and  BroydenSecondMethodMixing, both are 0(false). So Simple mixing is used"<<endl;

        }

    }



    Read_OPs_double=double(matchstring(inputfile_,"Read_initial_OPvalues"));
    if(Read_OPs_double==1.0){
        Read_OPs=true;
    }
    else{
        Read_OPs=false;
    }


    PBC_string=matchstring2(inputfile_,"PBC");
    if(PBC_string=="true"){
        PBC=true;
    }
    else{
        PBC=false;
    }

    File_OPs_in=matchstring2(inputfile_,"Read_initial_OPvalues_file");
    File_OPs_out=matchstring2(inputfile_,"Write_Final_OPvalues_file");


    Hopping_NN.resize(n_orbs,n_orbs);
    string Nearest_Neigh_Hopping;
    Nearest_Neigh_Hopping=matchstring2(inputfile_, "Nearest_Neigh_Hopping_matrix");
    stringstream Hopping_stream(Nearest_Neigh_Hopping);

    for(int orb0=0;orb0<n_orbs;orb0++){
        for(int orb1=0;orb1<n_orbs;orb1++){
        Hopping_stream >> Hopping_NN(orb0,orb1);
        }
    }


    string Crystal_Field_string;
    Crystal_Field_string=matchstring2(inputfile_, "Crystal_Field");
    stringstream Crystal_Field_stream(Crystal_Field_string);

    Crystal_Field.resize(n_orbs);
    for(int n=0;n<n_orbs;n++){
        Crystal_Field_stream >> Crystal_Field[n];
    }

    pi=4.00*atan(double(1.0));
    Eav=0.0;

    Temperature=0.001;
    //beta=(11605.0/Temperature);
     beta=(1.0/Temperature);

    mus=0.25;
    mu_old=0.25;
    cout << "____________________________________" << endl;
}


double Parameters::matchstring(string file,string match) {
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass=false;
    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass==false) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }
    if (pass==false) {
        string errorout=match;
        errorout+="= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters::matchstring2(string file,string match) {

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if(readFile.is_open())
    {
        while(!readFile.eof())
        {
            getline(readFile,line);

            if ((offset = line.find(match, 0)) != string::npos) {
                amount = line.substr (offset+match.length()+1);				}

        }
        readFile.close();
    }
    else
    {cout<<"Unable to open input file while in the Parameters class."<<endl;}




    cout << match << " = " << amount << endl;
    return amount;
}

#endif



