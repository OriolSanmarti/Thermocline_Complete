#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <unordered_map>
using namespace std;

const double M_PI = 3.141592;

// Fluid and solid properties
int casemoltensalt = 1; // 1 -> solar salt 2->yaramost 3->hitecxl 4->air 5->Jarytherm oil

double calcRhoFluid(double T){
        double value = 0;
        if(casemoltensalt == 1){
            value = 2090 - 0.636*T;
        }else if(casemoltensalt == 2){
            value = 2085 - 0.74*T;
        }
        else if(casemoltensalt == 3){
            value = 2240 - 0.8266*T;
        }
        else if(casemoltensalt == 4){
            value = 0.5;
        }
        else if(casemoltensalt == 5){
            T += 273;
            value = 1261.569 - 0.7419*(T);
        }
        return value;
}
double calcCpFluid(double T){
        double value = 0;
        if(casemoltensalt == 1){
            value = 1443 + 0.172*T;
        }else if(casemoltensalt == 2){
            value = 1549 - 0.15*T;
        }
        else if(casemoltensalt == 3){
            value = 1536 - 0.2624*T;
        }
        else if(casemoltensalt == 4){
            value = 1075;
        }
        else if(casemoltensalt == 5){
            T += 273;
            value = 649.84+3.18*(T);
        }
        return value;
}
double calcLambdaFluid(double T){
        double value = 0;
        if(casemoltensalt == 1){
            value = 0.443 + 1.93E-4*T;
        }else if(casemoltensalt == 2){
            value = 0.697 - 4.61E-4*T;
        }
        else if(casemoltensalt == 3){
            value = 0.519;
        }
        else if(casemoltensalt == 4){
            value = 0.05;
        }
        else if(casemoltensalt == 5){
            T += 273;
            value = 0.15-8.2E-5*(T);
        }
        return value;
}
double calcMuFluid(double T){
        double value = 0;
        if(casemoltensalt == 1){
            value = 22.714E-3 - 0.120E-3*T + 2.281E-7*pow(T, 2.0) -1.474E-10*pow(T, 3.0);
        }else if(casemoltensalt == 2){
            value = (31.59 - 0.1948*T + 0.000425*pow(T,2.0) - 0.0000002122*pow(T,3.0) )/ 1000.0;
        }
        else if(casemoltensalt == 3){
            value = 1372000*pow(T,-3.364);
        }
        else if(casemoltensalt == 4){
            value = 3.4E-5;
        }
        else if(casemoltensalt == 5){
            T += 273;
            double A = 19.75102*pow(log(T), 4) - 492.2114*pow(log(T), 3) + 4602.036*pow(log(T), 2) - 19136.34*log(T) + 29858.54;
            value = exp(A);

        }
        return value;
}

// Objectes

class numerical
{
public:
	int N;
	int M;
	double tmax; // s
	double inct; // s
	double UA;
	double Qi;
	double Qbot;
	double Qtop;
	double alphafb;//26337.75; // W/m^2
	double t; // s
	int ite;
	double g;
	double errorGS; // Gauss-Seidel error
};
class physical
{
public:
	double L; // m
	double D; // m
	double incz; // m
	double voidf; //0.22;  // Vf/Vs
	double voidfPCM;
	double Vi;
	double Vbi;
	double VbiPCM;
	double Vfi;
	double VfiPCM;
	double Sfi;
	double SfiPCM;
	double dpi; // 0.01338; // m 
	double dpenc;  // Diametre encapsulat // m
	double Ltube;
	double LPCMbot;
	double NodesPCMbot;
	double LPCMtop;
	double NodesPCMtop;	
	double Text; // ºC
	double htop;
	double hbottom;
	double Thot;
	double Tcold;
	double Tcutoff;
	double Tmid;
	double mfc; // kg/s
	double mfd; // kg/s
	double mf;
};

class Fluid
{
public:
    // Constructor
    Fluid(){}
    Fluid(int N, double Tinicial) : Tf(N, Tinicial), Tantf(N, Tinicial), P(N, 0.0), alpha(N, 0.0), rhof(N, 0.0), rhofant(N, 0.0), rho_velf(N, 0.0), cpf(N, 0.0), lambdaf(N, 0.0), muf(N, 0.0)
    {
        // Puedes realizar otras operaciones de inicialización aquí si es necesario
    }

    std::vector<double> Tf;
    std::vector<double> Tantf;
    vector<double> P;
	vector<double> alpha;
	vector<double> rhof;
	vector<double> rhofant;
	vector<double> rho_velf;
	vector<double> cpf;
	vector<double> lambdaf;
	vector<double> muf;
};

class Solid
{
public:
	// Constructor
    Solid(){}
    Solid(int N, double Tinicial) : Ts(N, Tinicial), Tants(N, Tinicial), rhos(N, 0.0), cps(N, 0.0), lambdas(N, 0.0)
    {
        // Puedes realizar otras operaciones de inicialización aquí si es necesario
    }
    vector<double> Ts;
	vector<double> Tants;
	vector<double> rhos;
	vector<double> cps;
	vector<double> lambdas;
};

class PCM
{
public:
	PCM(){}
	PCM(int N, int M, double Tinicial) : Tsolid(N, vector<double> (M,Tinicial)), Tantsolid(N, vector<double> (M,Tinicial)), h(N, vector<double> (M,0)), hant(N, vector<double> (M,0)), mass(N, vector<double> (M,0)), massant(N, vector<double> (M,0)), h1D(N, 0), h1Dant(N, 0), Mass1D(N, 0),  Mass1Dant(N, 0)
	{
		
	}

	std::vector<std::vector<double>> Tsolid;
    std::vector<std::vector<double>> Tantsolid;
    std::vector<std::vector<double>> h;
    std::vector<std::vector<double>> hant;
    std::vector<std::vector<double>> mass;
    std::vector<std::vector<double>> massant;
	vector<double> h1D;
	vector<double> h1Dant;
	vector<double> Mass1D;
	vector<double> Mass1Dant;
	
	double cpsentalpia = 1340; // 1340;  // 1223; // J/kgK
	double rhosentalpia = 2040; // 2040; // 2250; // kg/m3
	double cplentalpia = 1340; // 1340;  // 1665; // J/kgK
	double rholentalpia = 2040; // 2040; // 1900; // kg/m3
	//double masstotalpcm = VbiPCM * rhosentalpia; // kg
	double Lentalpia = 1.34E5; // 1.34E5;  // 1.75E5; // J/kg
	double lambdaPCM = 0.5; // W/mK
	double Tsl = 300; // 306; // ºC
	double Tfusio = Tsl-0.1; // ºC
	double Tl = Tsl + 1.0;   // ºC
	double Tref = 25;        // ºC
};

class FluidCell{
    public:
        FluidCell(){}
        FluidCell(double T_, double incz_, double Qloss_,double D_)
        {
            setT(T_);
            incz = incz_;
            Qloss = Qloss_;
            vol = incz_ * M_PI*D_*D_/4.0;
            D = D_;
        }
        void setT(double T_)
        {
            T = T_;
            rho = calcRhoFluid(T_);
            cp = calcCpFluid(T_);
            mu = calcMuFluid(T_);
            lambda = calcLambdaFluid(T_);
        }
        void updateProps()
        {
            Tant = T;
            rhoant = rho;
            cpant = cp;
            muant = mu;
            lambdaant = lambda;
            inczant = incz;
            vol = incz * M_PI*D*D/4.0;
        }
        double T;
        double Tant;
        double incz;
        double inczant;
        double rho;
        double rhoant;
        double cp;
        double cpant;
        double mu;
        double muant;
        double lambda;
        double lambdaant;
        double Qloss = 0.0;
        double v = 0.0;
        double vol;
        double rho_vel_up;
        double rho_vel_down;
        double D;
        double Tmid = (560+290) / 2.0;
};

class simulation
{
public:
	simulation(){}
	Fluid fluid;
	Solid solid;
	PCM PCM;
	numerical num;
	physical phy;	
	FluidCell topCV;
	FluidCell bottomCV;
	
	simulation(int N, int M, double Tinicial, double htop, double hbottom, double D) : fluid(N, Tinicial), solid(N, Tinicial), PCM(N, M, Tinicial), topCV(Tinicial, htop, 0.0, D), bottomCV(Tinicial, hbottom, 0.0, D)
    {
        // Puedes realizar otras operaciones de inicialización aquí si es necesario
    }
    
    int estat = 2;
	int estatanterior = 1;
	string PCMtype = "tube"; // "tube" or "PB"
	string fillerType = "ST"; // "tube" or "PB" or "ST"
	string Modeltype = "2D"; // "1D" or "2D"
};

void setupNumerical(simulation &sim, std::unordered_map<std::string, double> atributos)
{
	sim.num.N = static_cast<int>(atributos["N"]);
    sim.num.M = static_cast<int>(atributos["M"]);
    sim.num.tmax = atributos["tmax"];
    sim.num.inct = atributos["inct"];
	sim.num.UA = atributos["UA"];
	sim.num.Qi = atributos["Qi"];
	sim.num.Qbot = atributos["Qbot"];
	sim.num.Qtop = atributos["Qtop"];
	sim.num.alphafb = atributos["alphafb"];
	sim.num.t = atributos["t"];
	sim.num.ite = atributos["ite"];
	sim.num.g = atributos["g"];
	sim.num.errorGS = atributos["errorGS"]; 
}
void setupPhysical(simulation &sim, std::unordered_map<std::string, double> atributos)
{
	sim.phy.L = atributos["L"];
	sim.phy.D = atributos["D"];
	sim.phy.voidf = atributos["voidf"];
	sim.phy.voidfPCM = atributos["voidfPCM"];
	sim.phy.dpi = atributos["dpi"]; 
	sim.phy.dpenc = atributos["dpenc"]; 
	sim.phy.LPCMbot = atributos["LPCMbot"];
	sim.phy.LPCMtop = atributos["LPCMtop"];
	sim.phy.Text = atributos["Text"];
	sim.phy.htop = atributos["htop"];
	sim.phy.hbottom = atributos["hbottom"];
	sim.phy.Thot = atributos["Thot"];
	sim.phy.Tcold = atributos["Tcold"];
	sim.phy.Tcutoff = atributos["Tcutoff"];	
	sim.phy.mfc = atributos["mfc"];
	sim.phy.mfd = atributos["mfd"];
	sim.phy.mf = atributos["mf"];
	
	double L = sim.phy.L;
	double D = sim.phy.D;
	sim.phy.incz = L / (double)sim.num.N;
	double incz = sim.phy.incz;
	double voidf = sim.phy.voidf;
	double voidfPCM = sim.phy.voidfPCM;
	sim.phy.Vi = M_PI*D*D/4.0*incz;
	sim.phy.Vbi = M_PI*D*D/4.0*incz*(1-voidf);
	sim.phy.VbiPCM = M_PI*D*D/4.0*incz*(1-voidfPCM);
	sim.phy.Vfi = M_PI*D*D/4.0*incz*(voidf);
	sim.phy.VfiPCM = M_PI*D*D/4.0*incz*(voidfPCM);
	sim.phy.Sfi = M_PI*D*D/4.0*(voidf);
	sim.phy.SfiPCM = M_PI*D*D/4.0*(voidfPCM);
	sim.phy.Ltube = D;
	double LPCMbot = sim.phy.LPCMbot;
	sim.phy.NodesPCMbot = LPCMbot/L * sim.num.N;
	double LPCMtop = sim.phy.LPCMtop;
	sim.phy.NodesPCMtop = LPCMtop/L * sim.num.N;
	cout << "Sfi: " << sim.phy.Sfi << " dpi: " << sim.phy.dpi << " voidf: " << sim.phy.voidf << " Vfi: " << sim.phy.Vfi << " incz: " << sim.phy.incz << " Vbi: " << sim.phy.Vbi << " N: " << (double)sim.num.N<< endl;
	sim.phy.Tmid = (sim.phy.Thot + sim.phy.Tcold) / 2.0;
}


double porosity = 0.345; 

double calcLambdaSolid(double T){
        double lambdaraw = 0.963;
        return (1-porosity)*lambdaraw + porosity*calcLambdaFluid(T);
        //return 5.69;
        //return 0.57;
}
double calcRhoSolid(double T){
        double rhoraw = 3053.4;
        return (1-porosity)*rhoraw + porosity*calcRhoFluid(T);
        //return 2500;
        //return 2785;
}
double calcCpSolid(double T){
        double rhoraw = 3053.4;
        double cpraw = 1280;
        return (1.0/(((1-porosity)*rhoraw)+(porosity*calcRhoFluid(T)))) * ((1-porosity)*(cpraw*rhoraw) + (porosity*calcCpFluid(T)*calcRhoFluid(T)));
        //return 830.0;
        //return 1148;
}




void calcContinuity(simulation &sim)
{
	double Sfi = sim.phy.Sfi;
	double incz = sim.phy.incz;
	double voidf = sim.phy.voidf;
	double inct = sim.num.inct;
	int N = sim.num.N;
	double mf = sim.phy.mf;
	
    double rho_in = sim.bottomCV.rho;
    double vfluid_in = mf/(rho_in*Sfi);
    // Bottom control volume
    // vfluid_in is is we had the section of the tube, so we correct the velocity because at the
    // bottom of the first c.v we have a diferent section by a factor of voidFraction
    sim.bottomCV.rho_vel_down = rho_in*vfluid_in * voidf;
    sim.bottomCV.rho_vel_up = (1.0/voidf) * (sim.bottomCV.rho_vel_down - (sim.bottomCV.rho - sim.bottomCV.rhoant)/inct * sim.bottomCV.incz);
    // Bed
    sim.fluid.rho_velf[0] = sim.bottomCV.rho_vel_up;
    for (uint32_t i = 1; i < N + 1; i++){
        sim.fluid.rho_velf[i] = sim.fluid.rho_velf[i-1] - (sim.fluid.rhof[i-1] - sim.fluid.rhofant[i-1])/inct * incz;
    }
    // Top control volume
    sim.topCV.rho_vel_down = sim.fluid.rho_velf.back();

    //We assume bottom in- top out mass is the same
    sim.topCV.rho_vel_up = rho_in*vfluid_in * voidf;
    // We calculate the height variation
    sim.topCV.incz = (1.0/sim.topCV.rho) * (sim.topCV.rhoant*sim.topCV.inczant + inct *((sim.topCV.rho_vel_down*voidf - sim.topCV.rho_vel_up)));
}
/*
void calculateProperties(simulation &sim)
{
	int N = sim.num.N;
	double mf = sim.phy.mf;
	double Tmid = sim.phy.Tmid;
    for(int i = 0; i<N+1; i++)
    {
        sim.fluid.rhof[i] = calcRhoFluid(Tmid);
        sim.fluid.cpf[i] = calcCpFluid(Tmid);
        sim.fluid.lambdaf[i] = calcLambdaFluid(Tmid);
        sim.fluid.muf[i] = calcMuFluid(Tmid);
        if(i<N)
        {
            sim.solid.rhos[i] = calcRhoSolid(Tmid);
            sim.solid.cps[i] = calcCpSolid(Tmid);
            sim.solid.lambdas[i] = calcLambdaSolid(Tmid);
        }
    }
}*/
void calculateProperties(simulation &sim)
{
	int N = sim.num.N;
	double mf = sim.phy.mf;
	double Tmid = sim.phy.Tmid;
    for(int i = 0; i<N+1; i++)
    {
        sim.fluid.rhof[i] = calcRhoFluid(sim.fluid.Tf[i]);
        sim.fluid.cpf[i] = calcCpFluid(sim.fluid.Tf[i]);
        sim.fluid.lambdaf[i] = calcLambdaFluid(sim.fluid.Tf[i]);
        sim.fluid.muf[i] = calcMuFluid(sim.fluid.Tf[i]);
        if(i<N)
        {
            sim.solid.rhos[i] = calcRhoSolid(sim.solid.Ts[i]);
            sim.solid.cps[i] = calcCpSolid(sim.solid.Ts[i]);
            sim.solid.lambdas[i] = calcLambdaSolid(sim.solid.Ts[i]);
        }
    }
}
void actualizeProperties(simulation &sim)
{
	/*
    for(int i = 0; i<N; i++)
    {
        rhofant[i] = rhof[i];
        Tantf[i] = Tf[i];
		Tants[i] = Ts[i];
    }*/
    sim.fluid.rhofant = sim.fluid.rhof;
    sim.fluid.Tantf = sim.fluid.Tf;
    sim.solid.Tants = sim.solid.Ts;
    sim.PCM.Tantsolid = sim.PCM.Tsolid;
    sim.PCM.hant = sim.PCM.h;
    sim.PCM.h1Dant = sim.PCM.h1D;
    sim.PCM.massant = sim.PCM.mass;
}


void CalcLosses(double Tf, double &UA, simulation &sim)
{
	double D = sim.phy.D;
	double incz = sim.phy.incz;
	
    double U = 0.0;
    double e = 0.4;
    double lambdae = 0.086;
    double alphaSurfo = 95.49;
    double alphaSurfi = 4.0;
    double Pi = M_PI * D;
    double Po = M_PI * (D+2.0*e);
    double Do = D + 2.0*e;
    double Di = D;
    double Rfi = 0.0;
    double Rfo = 0.0;
    
    double alphaRfi = 1.0 / (Rfi + 1.0/alphaSurfi);
    double alphaRfo = 1.0 / (Rfo + 1.0/alphaSurfo);

    U = 1.0 / ( 1.0/alphaRfi * (Po/Pi) + Po/(2*M_PI*lambdae)*log(Do/Di) + 1/alphaRfo );

    UA = U * Po * incz;
    UA = 0.0;
}
double calcEnergyEntalpia(int x, simulation &sim)
{
	double rhosentalpia = sim.PCM.rhosentalpia;
	double masstotalpcm = sim.phy.VbiPCM * rhosentalpia;
	double dpi = sim.phy.dpi;
	double dpenc = sim.phy.dpenc;
	double Ltube = sim.phy.Ltube;
	double VbiPCM = sim.phy.VbiPCM;
	int M = sim.num.M;
	
	double Qentalpia = 0.0;
	if(sim.Modeltype == "1D")
	{
		Qentalpia = (sim.PCM.h1D[x] - sim.PCM.h1Dant[x]) * masstotalpcm;
	}
	else if(sim.Modeltype == "2D")
	{
		double Vp = 0.0;
		if (sim.PCMtype == "PB")Vp = 4.0/3.0 * M_PI * pow((dpi/2.0), 3);
		else if(sim.PCMtype == "tube") Vp = M_PI * pow(dpenc, 2.0) / 4.0 * Ltube; 
		// La h són J/kg, per tant:
		double summas = 0.0;
		for(int i = 0; i<M-1; i++)
		{
			Qentalpia += (sim.PCM.h[x][i] - sim.PCM.hant[x][i])* (sim.PCM.mass[x][i]); 
			summas += sim.PCM.mass[x][i];
			//cout <<"m: " << i << " h: " << h[x][i] << " hant: " << hant[x][i] << " mass: " << mass[x][i] << endl;
		}
		Qentalpia = Qentalpia / (Vp/VbiPCM);
		summas = summas / (Vp/VbiPCM);
		
		//Qentalpia = (h[x][0] - hant[x][0]) * masstotalpcm;	
		
		//cout << " masstotalpcm: " << masstotalpcm << " sumams: " << summas << endl;
	}
	
	return Qentalpia;
}
void calcAlpha(double i, string fillerType_, double voidf_, simulation &sim)
{
	double D = sim.phy.D;
	double dpi = sim.phy.dpi;
	double Vi = sim.phy.Vi;
	double dpenc = sim.phy.dpenc;
	double Ltube = sim.phy.Ltube;
	double mf = sim.phy.mf;
	double Tmid = sim.phy.Tmid;
	
	
	double Sfi_ = M_PI*D*D/4.0*(voidf_);
	if(fillerType_ == "PB")
	{
		double velfi = mf / (sim.fluid.rhof[i]*Sfi_);
		double Re = fabs(sim.fluid.rhof[i]*velfi*dpi/sim.fluid.muf[i]);
	    double Pr = sim.fluid.muf[i]*sim.fluid.cpf[i]/sim.fluid.lambdaf[i];
	    double Nu = 2.0+1.1*pow(Re,0.6)*pow(Pr,(1.0/3.0));
	    sim.fluid.alpha[i] = Nu*6*(1-voidf_)*sim.fluid.lambdaf[i] / (pow(dpi,2.0)); // W/m3
        //cout << " alpha: " << alpha[i] << endl;
		sim.fluid.alpha[i] =  sim.fluid.alpha[i] * Vi;							   // W
	}
	else if(fillerType_ == "tube")
	{
		double velfi = mf / (sim.fluid.rhof[i]*Sfi_);
        double STtube = 0.06;
        double SDtube = 0.06;
        velfi = velfi * (STtube/(SDtube-dpenc)); 
		double Re = fabs(sim.fluid.rhof[i]*velfi*dpenc/sim.fluid.muf[i]);
	    double Prf = sim.fluid.muf[i]*sim.fluid.cpf[i]/sim.fluid.lambdaf[i];
	    double mufs = calcMuFluid(sim.solid.Ts[i]);
	    double cpfs = calcCpFluid(sim.solid.Ts[i]);
	    double lambdafs = calcLambdaFluid(sim.solid.Ts[i]);
	    double Prs = mufs*cpfs/lambdafs;					
	    double Crow = 0.95;
	    double Nu = 0.51*Crow*pow(Re, 0.5)*pow(Prf, 0.37) * pow((Prf/Prs), 0.25);
	    double apcm = (1-voidf_) * (M_PI*dpenc*Ltube + 2*M_PI*pow(dpenc, 2.0)/4.0 ) / (M_PI*pow(dpenc, 2.0)/4.0*Ltube);		// m2/m3 surface area per unit volume of PCM
        //cout << apcm << endl;
	    sim.fluid.alpha[i] = Nu * sim.fluid.lambdaf[i] * apcm/ dpenc;  // W/m3
        //cout << "apcm: " << apcm <<" alpha: " << alpha[i] << endl;
		sim.fluid.alpha[i] = sim.fluid.alpha[i] * Vi;							    // W
	}
	else if(fillerType_ == "ST")
	{
		double d = 0.01;
		double Ltp = 0.02;
		double C = 3.46;
		double D = sim.phy.D;
		double Ntt = M_PI*D*D/(C*Ltp*Ltp);
		double Ptubs = M_PI*d*Ntt;
		double incz = sim.phy.incz;
	    double Nu = 5.10;		// 
	    sim.fluid.alpha[i] = Nu*sim.fluid.lambdaf[i] / d; // W/m2
		sim.fluid.alpha[i] =  sim.fluid.alpha[i] * Ptubs*incz;							   // W
	}
	else 
	{
		cout << "INCORRECT PCM TYPE" << endl;
		return;	
	}
}
void calcAlphaPCM(double i, simulation &sim)
{
	double SfiPCM = sim.phy.SfiPCM;
	double dpi = sim.phy.dpi;
	double voidfPCM = sim.phy.voidfPCM;
	double Vi = sim.phy.Vi;
	double VbiPCM = sim.phy.VbiPCM;
	double dpenc = sim.phy.dpenc;
	double Ltube = sim.phy.Ltube;
	double mf = sim.phy.mf;
	double Tmid = sim.phy.Tmid;
	
	if(sim.PCMtype == "PB")
	{
		double velfi = mf / (sim.fluid.rhof[i]*SfiPCM);
		double Prf = sim.fluid.muf[i]*sim.fluid.cpf[i]/sim.fluid.lambdaf[i];
		double Re = fabs(sim.fluid.rhof[i]*velfi*dpi/sim.fluid.muf[i]);
		double Nu = 2.0+1.1*pow(Re,0.6)*pow(Prf,(1.0/3.0));
	    double Vp = 4.0/3.0 * M_PI * pow((dpi/2.0), 3);
	    sim.fluid.alpha[i] = Nu*6*(1-voidfPCM)*sim.fluid.lambdaf[i] / (pow(dpi,2.0)); // W/m3	
	    sim.fluid.alpha[i] =  sim.fluid.alpha[i] * Vi * (Vp/VbiPCM);  // Return W
	    //alpha[i] = 0;
        //cout << "alpha: " << alpha[i] << endl;
        //cout << "Vi: " << Vi << " VfiPCM: " << VfiPCM << " Vfi: " << Vfi << " Aw: " << Aw << " Aw*VbiPCM: " << Aw/VbiPCM<<endl;
	}
	else if(sim.PCMtype == "tube")
	{
		double velfi = mf / (sim.fluid.rhof[i]*SfiPCM);
        double STtube = 0.06;
        double SDtube = 0.06;
        velfi = velfi * (STtube/(SDtube-dpenc)); 
		double Re = fabs(sim.fluid.rhof[i]*velfi*dpenc/sim.fluid.muf[i]);
	    double Prf = sim.fluid.muf[i]*sim.fluid.cpf[i]/sim.fluid.lambdaf[i];
	    double mufs = calcMuFluid(sim.solid.Ts[i]);
	    double cpfs = calcCpFluid(sim.solid.Ts[i]);
	    double lambdafs = calcLambdaFluid(sim.solid.Ts[i]);
	    double Prs = mufs*cpfs/lambdafs;					
	    double Crow = 0.95;
	    double Nu = 0.51*Crow*pow(Re, 0.5)*pow(Prf, 0.37) * pow((Prf/Prs), 0.25);
	    double Vp = M_PI * pow(dpenc, 2.0) / 4.0 * Ltube;     	// Volum d'un tub m3
	    double apcm = (1-voidfPCM) * (M_PI*dpenc*Ltube + 2*M_PI*pow(dpenc, 2.0)/4.0 ) / (M_PI*pow(dpenc, 2.0)/4.0*Ltube);		// m2/m3 surface area per unit volume of PCM
	    sim.fluid.alpha[i] = Nu * sim.fluid.lambdaf[i] * apcm/ dpenc;  // W/m3
		sim.fluid.alpha[i] =  sim.fluid.alpha[i] * Vi * (Vp/VbiPCM);		// Return W		
	}
	else
	{
		cout << "INCORRECT PCM TYPE" << endl;
		return;
	}
}
// PostProcess variables for heat (Q) accumulated
double Qacc_f = 0.0;
double Qacc_s = 0.0;
double Qacc_pcm = 0.0;
double Qaccbuffer = 0.0;
double Qaccsump = 0.0;
// PostProcess variables for Energy accumulated
double Eacc_f = 0.0;
double Eacc_s = 0.0;
double Eacc_pcm = 0.0;
double Eaccbuffer = 0.0;
double Eaccsump = 0.0;
// PostProcess variables for balances
double Eaccfb = 0.0;
double Eaccsb = 0.0;
double Econvfb = 0.0;
double Econvsb = 0.0;
double Emfb = 0.0;
double Eaccpcm = 0.0;
double Econvpcm = 0.0;
double Qaccfb = 0.0;
double Qaccsb = 0.0;
double Qconvfb = 0.0;
double Qconvsb = 0.0;
double Qaccpcm = 0.0;
double Qconvpcm = 0.0;
double Qmfb = 0.0;
double Q_mass_f = 0.0;
double E_mass_f = 0.0;
double Qloss_t = 0.0;
double Eloss_t = 0.0;
// PostProcess variables for Temperature among height
bool printedAll = false; 
double tempsCiclePrint = 50.0;  // The minimum elapsed hours to start printing the fluid temperature
double factmwh = 2.777777777E-7; // to convert from [kJ] to [MWh]
int contCicless = 0; // Cycle counter
double tant = 0; // [s] Elapsed time on the previous cycle
void PostProcess(simulation &sim)
{
	double NodesPCMbot = sim.phy.NodesPCMbot;
	double voidfPCM = sim.phy.voidfPCM;
	double VfiPCM = sim.phy.VfiPCM;
	double NodesPCMtop = sim.phy.NodesPCMtop;
	double Vfi = sim.phy.Vfi;
	double Vbi = sim.phy.Vbi;
	double incz = sim.phy.incz;
	double inct = sim.num.inct;
	int N = sim.num.N;
	int ite = sim.num.ite;
	double t = sim.num.t;
	double mf = sim.phy.mf;
	double mfc = sim.phy.mfc;
	double mfd = sim.phy.mfd;
	double Tmid = sim.phy.Tmid;
	
	// File header initialization
    ofstream Archivo;
    ofstream Archivo34;
    if(ite == 1){
    	Archivo.open("Output.txt", ios::trunc);
	    Archivo << "#temps[s] cicle Eacc_f[MWh] Eacc_S[MWh] Eacc_PCM[MWh] Eacc_t[MWh]" << endl;
	    Archivo.close();
	}
    // File Tout header initialization
    ofstream Archivo3;
    if(ite == 1){
    	Archivo3.open("Tout.txt", ios::trunc);
	    Archivo3 << "#Time Tf[N-1] Tf[3*N/4] Tf[N/2] Ts[N-1]" << endl;
	    Archivo3.close();
	}
	// NOU checkear	
	double mdotg = 0;
	if(sim.estat == 1) mdotg = -mfc;
    if(sim.estat == 2) mdotg = mfd;
	for(int i = 0; i<N; i++){
        //double Tfii = (Tf[i]+Tf[i+1]) / 2.0;
        //double Tfiiant = (Tantf[i]+Tantf[i+1]) / 2.0;
        //double rhofii = (rhof[i]+rhof[i+1]) / 2.0;
        //double cpfii = (cpf[i]+cpf[i+1]) / 2.0;
        double Tfii = sim.fluid.Tf[i];
        double Tfiiant = sim.fluid.Tantf[i];
        double rhofii = sim.fluid.rhof[i];
        double cpfii = sim.fluid.cpf[i];
        if (i < NodesPCMbot) // PCM bottom
		{
		   calcAlpha(i, sim.PCMtype, voidfPCM, sim);
           Qaccfb += rhofii*cpfii*(Tfii-Tfiiant)*VfiPCM/inct;
           Qaccpcm += calcEnergyEntalpia(i, sim)/inct;
           Qconvpcm += sim.fluid.alpha[i] * (sim.solid.Ts[i] - Tfii);
		}
        else if (i > (N - NodesPCMtop)) // PCM top
        {
           calcAlpha(i, sim.PCMtype, voidfPCM, sim);
           Qaccfb += rhofii*cpfii*(Tfii-Tfiiant)*VfiPCM/inct;
           Qaccpcm += calcEnergyEntalpia(i,sim)/inct;
           Qconvpcm += sim.fluid.alpha[i] * (Tfii - sim.solid.Ts[i]);
        }
        else // Filler
        {
        	Qaccfb += rhofii*cpfii*(Tfii-Tfiiant)*Vfi/inct;
        	Qconvfb += sim.fluid.alpha[i] * (sim.solid.Ts[i]-Tfii); 
        	Qaccsb += sim.solid.rhos[i]*sim.solid.cps[i]*(sim.solid.Ts[i]-sim.solid.Tants[i])*Vbi/inct;
        	Qconvsb += sim.fluid.alpha[i] * (Tfii-sim.solid.Ts[i]);
        }
	        
        //Fluid
        //Qaccfb += rhofii*cpfii*(Tfii-Tfiiant)*Vfi/inct;
        //Qconvfb += alpha[i] * (Ts[i]-Tfii);
        if(sim.estat == 1 and i != (N-1)) Qmfb += mdotg*cpfii*(sim.fluid.Tf[i+1]-sim.fluid.Tf[i]);
        else if(sim.estat == 2 and i != 0) Qmfb += mdotg*cpfii*(sim.fluid.Tf[i]-sim.fluid.Tf[i-1]);
        
        //Solid
        //Qaccsb += rhos[i]*cps[i]*(Ts[i]-Tants[i])*Vbi/inct;
        //Qconvsb += alpha[i] * (Tfii-Ts[i]) ;
        //Losses
        //CalcLosses(Tfii, UA);
        //Qloss_t += UA*(Tfii - Text); 
    }

	Eaccfb = Qaccfb*inct*0.001*factmwh;
    Econvfb = Qconvfb*inct*0.001*factmwh;
    Emfb = Qmfb*inct*0.001*factmwh;
    Eaccsb = Qaccsb*inct*0.001*factmwh;
    Econvsb = Qconvsb*inct*0.001*factmwh;
    Eloss_t = Qloss_t * inct*0.001*factmwh;   
    Eaccpcm = Qaccpcm * inct*0.001*factmwh; 
    Econvpcm = Qconvpcm * inct*0.001*factmwh;  
	Q_mass_f = mdotg * ((sim.fluid.cpf[N] * sim.fluid.Tf[N]) - (sim.fluid.Tf[0]*sim.fluid.cpf[0]));
    E_mass_f += Q_mass_f * inct * 0.001*factmwh;
    if(ite%100 == 0){
        Archivo34.open("Bal.txt", ios::app);
        if(sim.estat == 1) Archivo34 << t/3600.0 << " " << Eaccfb << " " << Econvfb << " " << Emfb << " " << Eaccsb << " " << Econvsb <<" " << Eloss_t << " " << E_mass_f << " " << Q_mass_f <<" "<< Eaccpcm << " " << Econvpcm<<endl;
        if(sim.estat == 2) Archivo34 << t/3600.0 << " " << Eaccfb << " " << Econvfb << " " << Emfb << " " << Eaccsb << " " << Econvsb << " " << Eloss_t << " " << E_mass_f << " " << Q_mass_f<<" "<< Eaccpcm << " " << Econvpcm<<endl;
        Archivo3.close();
    }
    // NOU
    // If state changes, we save the accumulated energy and the duration of the cycle
    if(sim.estat != sim.estatanterior){
        Archivo.open("Output.txt", ios::app);
        Archivo << (t-tant)/3600.0 << " "  <<  sim.estatanterior << " " << Eacc_f + Eaccbuffer + Eaccsump << " " << Eacc_s << " " << Eacc_pcm << " "<< (Eacc_f+Eacc_s) + Eaccbuffer + Eaccsump + Eacc_pcm <<endl;
        tant = t;
    }
    Archivo.close();
    
    // Energy balances
    if(sim.estat != sim.estatanterior){ // We want to know how much energy is accumulated each cycle, then we set to 0 at the start of each cycle
    	Qacc_f = 0.0, Qacc_s = 0, Qacc_pcm = 0;
    	Qaccbuffer = 0.0, Qaccsump = 0.0;
	}

	// The accumulated energy for the fluid and solid is calculated as Qacc = rho*cp*(T - Tant)*Volume / time-step [W]
    for(int i = 0; i<N; i++){
    	 if (i < NodesPCMbot) // PCM bottom
			{
               Qacc_f += sim.fluid.rhof[i]*sim.fluid.cpf[i]*(sim.fluid.Tf[i]-sim.fluid.Tantf[i])*VfiPCM/inct;
               Qacc_pcm += calcEnergyEntalpia(i, sim)/inct;
			}
            else if (i > (N - NodesPCMtop)) // PCM top
            {
               Qacc_f += sim.fluid.rhof[i]*sim.fluid.cpf[i]*(sim.fluid.Tf[i]-sim.fluid.Tantf[i])*VfiPCM/inct;
               Qacc_pcm += calcEnergyEntalpia(i,sim)/inct;
            }
	        else // Filler
	        {
	        	Qacc_f += sim.fluid.rhof[i]*sim.fluid.cpf[i]*(sim.fluid.Tf[i]-sim.fluid.Tantf[i])*Vfi/inct;
        		Qacc_s += sim.solid.rhos[i]*sim.solid.cps[i]*(sim.solid.Ts[i]-sim.solid.Tants[i])*Vbi/inct; 
	        }
    }

    // Buffer and Sump accumulated heat calculated as Qacc = rho*cp*(T-Tant)*CVvol / time-step [W]
    double rhofi = sim.topCV.rho;
    double cpfi = sim.topCV.cp;
    Qaccbuffer += rhofi*cpfi*(sim.topCV.T - sim.topCV.Tant)*sim.topCV.vol/inct;
     rhofi = sim.bottomCV.rho;
     cpfi = sim.bottomCV.cp;
    Qaccsump += rhofi*cpfi*(sim.bottomCV.T - sim.bottomCV.Tant)*sim.bottomCV.vol/inct; 

    // Fluid and solid accumulated heat converted from [W] to [MWh]
    Eacc_f = Qacc_f*inct*0.001*factmwh; 
    Eacc_s = Qacc_s*inct*0.001*factmwh; 
    Eacc_pcm = Qacc_pcm*inct*0.001*factmwh;
    
	// Sump and buffer accumulated heat converted from [W] to [MWh]
    Eaccbuffer = Qaccbuffer*inct*0.001*factmwh;
    Eaccsump = Qaccsump*inct*0.001*factmwh;

    //------------------------------------------------------------------------------------------------------------------------------------------------------//
         
	// Here starts the PostProcess to keep the temeprature among height of a converged cycle
	// We find the start of the first cycle after de required time has elapsed
    if (sim.estat != sim.estatanterior && t/3600.0 >= tempsCiclePrint){
        printedAll = true;
        contCicless++; // To track how many cycles we have printed
        if(contCicless > 2) printedAll = false; // If we have printed 2 cycles we don't need printing anymore
    }
    
    if (printedAll){ // If the flag is true, we print the temperature of the cycle
        int everyIte = 100; // If we print every iteration it will make the code really slow and generate so much unnecessary information
        if(ite%everyIte == 0){ // We check that we are on the correct iteration to print
            std::ofstream Archivo2;
            char nombre[20000]; // This is a strategy to open a file with a generated number 
            int numero = ite/everyIte;
            if(sim.estat == 1) sprintf(nombre, "tfluid/Tc_%05d.csv", numero); // The filename will be Tc_nnnnn.csv if the cycle is charging
            if(sim.estat == 2) sprintf(nombre, "tfluid/Td_%05d.csv", numero); // The filename will be Tc_nnnnn.csv if the cycle is discharging
            Archivo2.open(nombre, std::ios::trunc); // Open the file that have been created
            Archivo2 << "#Time: " << t/3600 << " hours" <<  " estatglobal: " << sim.estat << std::endl;
            Archivo2 << "#z[m] Tf[C] Ts[C]" << std::endl;
            for(int i = 0; i<N; i++){
                Archivo2 << incz*i << " " << sim.fluid.Tf[i] << " " << sim.solid.Ts[i] << std::endl; // Print the fluid and solid temperature 
            } 
            Archivo.close();
        }
        ofstream Archivo3;
        Archivo3.open("Tout.txt", ios::app);
        Archivo3 << t/3600.0 << " " << sim.fluid.Tf[N-1] << " " << sim.fluid.Tf[N/2] <<" " << sim.fluid.Tf[30] << " " << sim.fluid.Tf[0]<< std::endl;
	    Archivo3.close();
    }

    //-------------------------------------------------------------------------------------------------------------------------------------------------------//
    
    

}

void InitializeFiles()
{
    ofstream Archivo;
    Archivo.open("Output.txt", ios::trunc);
    Archivo << "temps[s] cicle Eacc_f[MW] Eacc_S[MW] Eacc_t[MW]" << endl;
    Archivo.close();
    ofstream Archivo3;
    Archivo3.open("Bal.txt", ios::trunc);
    Archivo3 << "#Qaccf Qconvf Qmf Qaccs Qconvs Qlosst EaccPCM[MWh] EconvPCM[MWh]" << endl;
    Archivo3.close();

}

void CalculateAlphaInit(simulation &sim)
{
	double dpenc = sim.phy.dpenc;
	double dpi = sim.phy.dpi;
	double Ltube = sim.phy.Ltube;
	double VbiPCM = sim.phy.VbiPCM;
	double Sfi = sim.phy.Sfi;
	double voidf = sim.phy.voidf;
	double mf = sim.phy.mf;
	double mfc = sim.phy.mfc;
	double mfd = sim.phy.mfd;
	double Tmid = sim.phy.Tmid;
	double Thot = sim.phy.Thot;
	double Tcold = sim.phy.Tcold;
	
    double rhohot = calcRhoFluid(Thot);
    double rhocold = calcRhoFluid(Tcold);
    double lambdahot = calcLambdaFluid(Thot);
    double lambdacold = calcLambdaFluid(Tcold);
    double muhot = calcMuFluid(Thot);
    double mucold = calcMuFluid(Tcold);
    double cphot = calcCpFluid(Thot);
    double cpcold = calcCpFluid(Tcold);
    double velhot = mfc /rhohot/Sfi;
    double velcold = mfc /rhocold/Sfi;
    double Rehot = rhohot*velhot*dpi/muhot;
    double Recold = rhocold*velcold*dpi/mucold;
    double Prhot = muhot*cphot/lambdahot;
    double Prcold = mucold*cpcold/lambdacold;
    double Nuhot = 2 + 1.1*pow(Rehot, 0.6)*pow(Prhot, 1.0/3.0);
    double Nucold = 2 + 1.1*pow(Recold, 0.6)*pow(Prcold, 1.0/3.0);
    double alphahot = Nuhot * 6 * (1-voidf) / (pow(dpi, 2.0));
    double alphacold = Nucold * 6 * (1-voidf) / (pow(dpi, 2.0));
    double alphamid = (alphahot + alphacold) * 0.5;
    cout << "alphahot: " << alphahot << endl;
    cout << "alphacold: " << alphacold << endl;
    cout << "alphamid: " << alphamid << endl;
    cout << "Nuhot: " << Nuhot << " Prhot: " << Prhot << " Rehot: " << Rehot << endl;
}

void ChooseChargeDischarge(simulation &sim)
{
	double Thot = sim.phy.Thot;
	double Tcold = sim.phy.Tcold;
	double Tcutoff = sim.phy.Tcutoff;
    if(sim.estat == 1){
        //if(Tf[0] >= Tcold + Tcutoff) estat = 2; //Criteri primer node fluid
        if(sim.bottomCV.T >= Tcold + Tcutoff) sim.estat = 2; //Criteri buffer zone
        //if(t > 0.0*3600) estat = 2; //Criteri buffer zone
        else sim.estat = 1;
        //cout << "sim.estat " << sim.estat <<endl;
    }
    else if(sim.estat == 2){
        //if(Tf[N] <= Thot - Tcutoff) estat = 1; //Criteri ultim node fluid
        if(sim.topCV.T <= Thot - Tcutoff) sim.estat = 1;
        //if(topCV.T <= Thot - Tcutoffaux) estat = 1;
        else sim.estat = 2;
        //cout << "sim.estat "  << sim.estat << endl;
    }
}


void initentalpia(simulation &sim)
{
	double cpsentalpia = sim.PCM.cpsentalpia; // 1340;  // 1223; // J/kgK
	double rhosentalpia = sim.PCM.rhosentalpia; // 2040; // 2250; // kg/m3
	double cplentalpia = sim.PCM.cplentalpia; // 1340;  // 1665; // J/kgK
	double rholentalpia = sim.PCM.rholentalpia; // 2040; // 1900; // kg/m3
	//double masstotalpcm = VbiPCM * rhosentalpia; // kg
	double Lentalpia = sim.PCM.Lentalpia; // 1.34E5;  // 1.75E5; // J/kg
	double lambdaPCM = sim.PCM.lambdaPCM; // W/mK
	double Tsl = sim.PCM.Tsl; // 306; // ºC
	double Tfusio = sim.PCM.Tfusio; // ºC
	double Tl = sim.PCM.Tl;   // ºC
	double Tref = sim.PCM.Tref;        // ºC
	int N = sim.num.N;
	int M = sim.num.M;
	
    double fracliquid = 0;
    double hmaxsolid = cpsentalpia*(Tfusio - Tref);
    double hmaxsl = hmaxsolid + cpsentalpia*(Tsl - Tfusio) + ((Tsl-Tfusio) / (Tl-Tfusio))*Lentalpia;
    double hminliquid = cpsentalpia*(Tsl-Tref) + cplentalpia*(Tl-Tsl) + Lentalpia;
    cout << "hmaxsolid:  " << hmaxsolid << " hmaxsl: " << hmaxsl << " hminliquid: " << hminliquid <<endl;

	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			double Tp = sim.PCM.Tsolid[i][j];
			fracliquid = (Tp-Tfusio) / (Tl-Tfusio);
			// Calculem on estem
		    if(Tp <= Tfusio)
		    {
		    	sim.PCM.h[i][j] = cpsentalpia * (Tp-Tref);
		    	sim.PCM.hant[i][j] = sim.PCM.h[i][j];
		    	sim.PCM.h1D[i] = cpsentalpia * (Tp-Tref);
		    	sim.PCM.h1Dant[i] = sim.PCM.h1D[i];
			}
			else if(Tp > Tfusio and Tp < Tsl)
			{
				sim.PCM.h[i][j] = cpsentalpia * (Tp-Tref) + fracliquid*Lentalpia;
				sim.PCM.hant[i][j] = sim.PCM.h[i][j];
				sim.PCM.h1D[i] = cpsentalpia * (Tp-Tref) + fracliquid*Lentalpia;
				sim.PCM.h1Dant[i] = sim.PCM.h1D[i];
			}
			else if(Tp > Tsl and Tp < Tl)
			{
				sim.PCM.h[i][j] = cpsentalpia * (Tsl-Tref) + cplentalpia*(Tp-Tsl) + fracliquid*Lentalpia;
				sim.PCM.hant[i][j] = sim.PCM.h[i][j];
				sim.PCM.h1D[i] = cpsentalpia * (Tsl-Tref) + cplentalpia*(Tp-Tsl) + fracliquid*Lentalpia;
				sim.PCM.h1Dant[i] = sim.PCM.h1D[i];
			}
			else
			{
				sim.PCM.h[i][j] = cpsentalpia*(Tsl-Tref) + cplentalpia*(Tp-Tsl) + Lentalpia;
				sim.PCM.hant[i][j] = sim.PCM.h[i][j];
				sim.PCM.h1D[i] = cpsentalpia*(Tsl-Tref) + cplentalpia*(Tp-Tsl) + Lentalpia;
				sim.PCM.h1Dant[i] = sim.PCM.h1D[i];
			}
			//cout << "i: " << i <<" j: " << j << " h: " << h[i][j] << endl;
		}
	}

}

void initmass(simulation &sim)
{
	double cpsentalpia = sim.PCM.cpsentalpia; // 1340;  // 1223; // J/kgK
	double rhosentalpia = sim.PCM.rhosentalpia; // 2040; // 2250; // kg/m3
	double cplentalpia = sim.PCM.cplentalpia; // 1340;  // 1665; // J/kgK
	double rholentalpia = sim.PCM.rholentalpia; // 2040; // 1900; // kg/m3
	//double masstotalpcm = VbiPCM * rhosentalpia; // kg
	double Lentalpia = sim.PCM.Lentalpia; // 1.34E5;  // 1.75E5; // J/kg
	double lambdaPCM = sim.PCM.lambdaPCM; // W/mK
	double Tsl = sim.PCM.Tsl; // 306; // ºC
	double Tfusio = sim.PCM.Tfusio; // ºC
	double Tl = sim.PCM.Tl;   // ºC
	double Tref = sim.PCM.Tref;        // ºC
	double dpenc = sim.phy.dpenc;
	double dpi = sim.phy.dpi;
	double Ltube = sim.phy.Ltube;
	double VbiPCM = sim.phy.VbiPCM;
	int N = sim.num.N;
	int M = sim.num.M;
	
    double incr = 0.0;
    if(sim.PCMtype == "tube") incr = (dpenc/(M-1)) / 2.0;
	if(sim.PCMtype == "PB") incr = (dpi/(M-1)) / 2.0;

	for(int j = 0; j<N; j++)
	{
		double masstotaltube = 0.0;
		double Vttube = 0.0;
		for(int i = 0; i<M-1; i++)
		{
			// Valors fisics
			double dest = i * incr + incr;
			double dwest = i * incr;
            double Vp = 0;
            if(sim.PCMtype == "PB") Vp = 4.0/3.0 * M_PI * (pow(dest, 3.0) - pow(dwest, 3.0));
            if(sim.PCMtype == "tube") Vp = M_PI * (pow(dest, 2.0) - pow(dwest, 2.0)) * Ltube;
            Vttube += Vp;
            //cout <<"i: "<< i <<" Vp: " << Vp << " Vt: " << Vttube<< endl;
			sim.PCM.mass[j][i] = Vp * rhosentalpia;
			sim.PCM.massant[j][i] = sim.PCM.mass[j][i];
			masstotaltube += sim.PCM.mass[j][i];
		}
		sim.PCM.Mass1D[j] = rhosentalpia * VbiPCM;
		sim.PCM.Mass1Dant[j] = sim.PCM.Mass1D[j];
		//cout << "total volume tube: " << Vttube <<"mass total tube: " << masstotaltube << endl;
	}
}

double calcentalpia(double h, simulation &sim)
{
	double cpsentalpia = sim.PCM.cpsentalpia; // 1340;  // 1223; // J/kgK
	double rhosentalpia = sim.PCM.rhosentalpia; // 2040; // 2250; // kg/m3
	double cplentalpia = sim.PCM.cplentalpia; // 1340;  // 1665; // J/kgK
	double rholentalpia = sim.PCM.rholentalpia; // 2040; // 1900; // kg/m3
	//double masstotalpcm = VbiPCM * rhosentalpia; // kg
	double Lentalpia = sim.PCM.Lentalpia; // 1.34E5;  // 1.75E5; // J/kg
	double lambdaPCM = sim.PCM.lambdaPCM; // W/mK
	double Tsl = sim.PCM.Tsl; // 306; // ºC
	double Tfusio = sim.PCM.Tfusio; // ºC
	double Tl = sim.PCM.Tl;   // ºC
	double Tref = sim.PCM.Tref;        // ºC
	
	double Tp = 0;
    double fracliquid = 0;
    double hmaxsolid = cpsentalpia*(Tfusio - Tref);
    double hmaxsl = hmaxsolid + cpsentalpia*(Tsl - Tfusio) + ((Tsl-Tfusio) / (Tl-Tfusio))*Lentalpia;
    double hminliquid = cpsentalpia*(Tsl-Tref) + cplentalpia*(Tl-Tsl) + Lentalpia;
    //cout << "hmaxsolid:  " << hmaxsolid << " hmaxsl: " << hmaxsl << " hminliquid: " << hminliquid <<endl;

    // Calculem on estem
    if(h <= hmaxsolid)
    {
        Tp = Tref + h/cpsentalpia;
        //cout << " Em sento solid a temperatura: " << Tp <<endl;
    }
    else if(h > hmaxsolid and h < hmaxsl)
    {
        Tp = (h + cpsentalpia*Tref + (Tfusio/(Tl-Tfusio))*Lentalpia) / (cpsentalpia + (1/(Tl-Tfusio))*Lentalpia);
        fracliquid = (Tp-Tfusio) / (Tl-Tfusio);
        //cout << " Em sento entre solid i sl a temperatura: " << Tp << " amb fracliquid: " << fracliquid << endl;
    }
    else if(h >= hmaxsl and h < hminliquid)
    {
        Tp = (h - cpsentalpia*(Tsl-Tref) + Tsl*cplentalpia + (Tfusio/(Tl-Tfusio))*Lentalpia) / (cplentalpia + (1/(Tl-Tfusio))*Lentalpia);
        fracliquid = (Tp-Tfusio) / (Tl-Tfusio);
        //cout << " Em sento entre sl i liquid a temperatura: " << Tp << " amb fracliquid: " << fracliquid << endl;
    }
    else
    {
        Tp = (h - Lentalpia + cplentalpia*Tsl - cpsentalpia*(Tsl-Tref)) / cplentalpia;
        //cout << " Em sento liquid a temperatura: " << Tp << endl;
    }
    
    return Tp;
}
void ChargeSolidEntalpiaCopia(int x, simulation &sim)
{
	double cpsentalpia = sim.PCM.cpsentalpia; // 1340;  // 1223; // J/kgK
	double rhosentalpia = sim.PCM.rhosentalpia; // 2040; // 2250; // kg/m3
	double cplentalpia = sim.PCM.cplentalpia; // 1340;  // 1665; // J/kgK
	double rholentalpia = sim.PCM.rholentalpia; // 2040; // 1900; // kg/m3
	//double masstotalpcm = VbiPCM * rhosentalpia; // kg
	double Lentalpia = sim.PCM.Lentalpia; // 1.34E5;  // 1.75E5; // J/kg
	double lambdaPCM = sim.PCM.lambdaPCM; // W/mK
	double Tsl = sim.PCM.Tsl; // 306; // ºC
	double Tfusio = sim.PCM.Tfusio; // ºC
	double Tl = sim.PCM.Tl;   // ºC
	double Tref = sim.PCM.Tref;        // ºC
	
	double dpenc = sim.phy.dpenc;
	double dpi = sim.phy.dpi;
	double Ltube = sim.phy.Ltube;
	double VbiPCM = sim.phy.VbiPCM;
	double inct = sim.num.inct;
	int N = sim.num.N;
	int M = sim.num.M;
	double errorGS = sim.num.errorGS;
	
    // Info entalpia
    double fracliquid = 0;
    double incr = 0.0;
    if(sim.PCMtype == "tube") incr = (dpenc/(M-1)) / 2.0;
    if(sim.PCMtype == "PB") incr = (dpi/(M-1)) / 2.0;
	bool itefinGS = false;	
	double beta = 0.5;
    int contiteracions = 0;
	while(!itefinGS){
		itefinGS = true;
        contiteracions += 1;
	    //for(int i = 0; i<M; i++)
	    for(int i = M-1; i>=0; i--)
	    {
	    	double Tsolidactual = sim.PCM.Tsolid[x][i];
	        // Valors fisics
	        double dpws = incr;
	        double dpes = incr;
	        double dest = i * incr + incr;
	        double dwest = i * incr ;
	        if(i == 0) dpes =  incr + incr/2.0;
	        if(i == M-1) dpws = incr/2.0;
	        if(i == M-1) dwest = i*incr - incr/2.0;
	        double Sws = 0;
	        double Ses = 0;
	        double Vp = 0;
	        double Aw = 0;
	        if(sim.PCMtype == "PB")
	        {
		        Sws = pow(dwest, 2.0) * 4 * M_PI; 
		        Ses = pow(dest, 2.0) * 4 * M_PI;
				Vp = 4.0/3.0 * M_PI * (pow(dest, 3.0) - pow(dwest, 3.0));
				Aw = 4.0/3.0 * M_PI * pow((dpi/2.0), 3);
			}
			else if(sim.PCMtype == "tube")
			{
				Sws = 2*dwest * M_PI * Ltube; // m2
		        Ses = 2*dest * M_PI * Ltube; // m2
				Vp = M_PI * (pow(dest, 2.0) - pow(dwest, 2.0)) / 1.0 * Ltube;
				Aw = M_PI*dpenc * Ltube;
				//cout << "i: " << i << " dpw: " << dpws << " dpes: " << dpes << " dest: " << dest << " dwest: " << dwest << " Sws: " << Sws << " Ses: " << Ses << endl;
			}
			else
			{
				cout << "INCORRECT PCM TYPE" << endl;
				return;
			}

	        if(i == (M-1))
	        {
	            double aws = lambdaPCM * Sws / (dpws);
	            //double bps = alpha[x]*Aw;
	            double bps = sim.fluid.alpha[x];
	            //
	            if(sim.PCMtype == "tube") {
	            	bps = sim.fluid.alpha[x]/Aw;	
	            	aws = lambdaPCM / (dpws);
				}
				//
	            sim.PCM.Tsolid[x][i] = (beta * (aws*sim.PCM.Tsolid[x][i-1] + bps*sim.fluid.Tf[x]) + (1-beta)*(-aws*sim.PCM.Tantsolid[x][i] + aws*sim.PCM.Tantsolid[x][i-1] - bps*sim.PCM.Tantsolid[x][i] + bps*sim.fluid.Tantf[x])) / (beta*(aws + bps));
	            //cout <<"x: " << x << " i: " << i << " T: " << Tsolid[x][i] << " Tant: " << Tantsolid[x][i] << endl;
	        }
	        else
	        {
	        	// t = n+1
	        	double aws = 0;
	        	if(i != 0) aws = lambdaPCM*Sws/dpws * (sim.PCM.Tsolid[x][i] - sim.PCM.Tsolid[x][i-1]);
	        	double aes = lambdaPCM*Ses/dpes * (sim.PCM.Tsolid[x][i+1] - sim.PCM.Tsolid[x][i]);
	        	// t = n
	        	double awsant = 0;
	        	if(i != 0) awsant = lambdaPCM*Sws/dpws * (sim.PCM.Tantsolid[x][i] - sim.PCM.Tantsolid[x][i-1]);
	        	double aesant = lambdaPCM*Ses/dpes * (sim.PCM.Tantsolid[x][i+1] - sim.PCM.Tantsolid[x][i]);
	        	
				sim.PCM.h[x][i] = (sim.PCM.massant[x][i]*sim.PCM.hant[x][i] + inct * (beta*(-aws + aes) + (1-beta)*(-awsant + aesant))) / sim.PCM.mass[x][i];
	            sim.PCM.Tsolid[x][i] = calcentalpia(sim.PCM.h[x][i], sim);
	            //cout <<"x: " << x << " i: " << i << " h: " << h[x][i] << " T: " << Tsolid[x][i] << " hant: " << hant[x][i] << " aws: " << aws << " aes: " << aes << " t: " << t / 3600.0 << endl;
	        }
	        
	        if(errorGS < fabs(sim.PCM.Tsolid[x][i]-Tsolidactual)) itefinGS = false;
	    }
        if(contiteracions > 10) itefinGS = true;
    
	}
	sim.solid.Ts[x] = sim.PCM.Tsolid[x][M-1];
}
void ChargeSolidEntalpia(int x, simulation &sim)
{
	double cpsentalpia = sim.PCM.cpsentalpia; // 1340;  // 1223; // J/kgK
	double rhosentalpia = sim.PCM.rhosentalpia; // 2040; // 2250; // kg/m3
	double cplentalpia = sim.PCM.cplentalpia; // 1340;  // 1665; // J/kgK
	double rholentalpia = sim.PCM.rholentalpia; // 2040; // 1900; // kg/m3
	//double masstotalpcm = VbiPCM * rhosentalpia; // kg
	double Lentalpia = sim.PCM.Lentalpia; // 1.34E5;  // 1.75E5; // J/kg
	double lambdaPCM = sim.PCM.lambdaPCM; // W/mK
	double Tsl = sim.PCM.Tsl; // 306; // ºC
	double Tfusio = sim.PCM.Tfusio; // ºC
	double Tl = sim.PCM.Tl;   // ºC
	double Tref = sim.PCM.Tref;        // ºC
	
	double dpenc = sim.phy.dpenc;
	double dpi = sim.phy.dpi;
	double Ltube = sim.phy.Ltube;
	double VbiPCM = sim.phy.VbiPCM;
	double inct = sim.num.inct;
	int N = sim.num.N;
	int M = sim.num.M;
	double errorGS = sim.num.errorGS;
	
    // Info entalpia
    double fracliquid = 0;
    double incr = 0.0;
    if(sim.PCMtype == "tube") incr = (dpenc/(M-1)) / 2.0;
    if(sim.PCMtype == "PB") incr = (dpi/(M-1)) / 2.0;
	bool itefinGS = false;	
    int contiteracions = 0;
	while(!itefinGS){
		itefinGS = true;
        contiteracions += 1;
        double volumtotalpcm = 0.0, masstotalpcm = 0.0;
	    for(int i = M-1; i>=0; i--)
	    {
	    	double Tsolidactual = sim.PCM.Tsolid[x][i];
	    	double hactual = sim.PCM.h[x][i];
	        // Valors fisics
	        double dpws = incr;
	        double dpes = incr;
	        double rest = i * incr + incr;
	        double rwest = i * incr ;
	        if(i == 0) dpes =  incr + incr/2.0;
	        if(i == 1) dpws = incr + incr/2.0; //
	        if(i == M-1) dpws = incr/2.0;
	        if(i == M-2) dpes = incr/2.0; //
	        double Sws = 0;
	        double Ses = 0;
	        double Vp = 0;
	        double Aw = 0;
	        if(sim.PCMtype == "PB")
	        {
		        Sws = pow(rwest, 2.0) * 4 * M_PI; 
		        Ses = pow(rest, 2.0) * 4 * M_PI;
				Vp = 4.0/3.0 * M_PI * (pow(rest, 3.0) - pow(rwest, 3.0));
				if (i != M-1)volumtotalpcm += Vp;
				if (i != M-1)masstotalpcm += sim.PCM.mass[x][i];
				Aw = 4 * M_PI * pow((dpi/2.0), 2.0);
			}
			else if(sim.PCMtype == "tube")
			{
				Sws = 2*rwest * M_PI * Ltube; // m2
		        Ses = 2*rest * M_PI * Ltube; // m2
				Vp = M_PI * (pow(rest, 2.0) - pow(rwest, 2.0)) / 1.0 * Ltube;
				Aw = M_PI*dpenc * Ltube;
				//cout << "i: " << i << " dpw: " << dpws << " dpes: " << dpes << " dest: " << rest << " dwest: " << rwest << " Sws: " << Sws << " Ses: " << Ses << endl;
			}
			else
			{
				cout << "INCORRECT PCM TYPE" << endl;
				std::exit(0);
			}

	        if(i == (M-1))
	        {
	            //double aws = lambdaPCM * Sws / (dpws);
	            //double bps = alpha[x];
	            double aws = lambdaPCM / (dpws);
	            double bps = sim.fluid.alpha[x]/Aw;
	            
	            sim.PCM.Tsolid[x][i] =  (aws*sim.PCM.Tsolid[x][i-1] + bps*sim.fluid.Tf[x])  / (aws + bps);
	            //Tsolid[x][i] = Tcold;
	            //cout <<"x: " << x << " i: " << i << " h: " << h[x][i] << " T: " << Tsolid[x][i] << " bps: " << bps <<" aws: " << aws << " t: " << t / 3600.0 << endl;
	        }
	        else if(i == 0)
	        {
	        	double aes = lambdaPCM*Ses/dpes * (sim.PCM.Tsolid[x][i+1] - sim.PCM.Tsolid[x][i]);
	        	
				sim.PCM.h[x][i] = (sim.PCM.massant[x][i]*sim.PCM.hant[x][i] + inct * (aes)) / sim.PCM.mass[x][i];
	            sim.PCM.Tsolid[x][i] = calcentalpia(sim.PCM.h[x][i], sim);
	            
	            //cout <<"x: " << x << " i: " << i << " h: " << h[x][i] << " T: " << Tsolid[x][i] << " hant: " << hant[x][i] << " Massant: " << massant[x][i]<<" aes: " << aes << " t: " << t / 3600.0 << endl;
			}
	        else
	        {
	        	double aws = lambdaPCM*Sws/dpws * (sim.PCM.Tsolid[x][i] - sim.PCM.Tsolid[x][i-1]);
	        	double aes = lambdaPCM*Ses/dpes * (sim.PCM.Tsolid[x][i+1] - sim.PCM.Tsolid[x][i]);
	        	
				sim.PCM.h[x][i] = (sim.PCM.massant[x][i]*sim.PCM.hant[x][i] + inct * (-aws + aes)) / sim.PCM.mass[x][i];
	            sim.PCM.Tsolid[x][i] = calcentalpia(sim.PCM.h[x][i], sim);
	            //cout <<"x: " << x << " i: " << i << " h: " << h[x][i] << " T: " << Tsolid[x][i] << " Tsolide: " << Tsolid[x][i+1]<< " Tsolidw: " << Tsolid[x][i-1] << " hant: " << hant[x][i] << " Massant: " << massant[x][i]<< " aws: " << aws << " aes: " << aes << " t: " << t / 3600.0 << endl;
	            //cout << "Sws: " << Sws << " dpws: " << dpws << " Ses: " << Ses << " dpes: " << dpes << endl;
	        }
	        
	        if(errorGS < fabs(sim.PCM.Tsolid[x][i]-Tsolidactual) or errorGS < fabs(sim.PCM.h[x][i])-hactual) itefinGS = false;
	        //cout << itefinGS << endl;
	    }
	    //cout << "Volum PCM: " << volumtotalpcm << "  volum analitic: "<< 4.0/3.0 * M_PI * pow((dpi/2.0), 3) << " mass total: " << masstotalpcm << " massanalitica: " << 4.0/3.0 * M_PI * pow((dpi/2.0), 3) * rhosentalpia<< " ite: " << contiteracions <<endl;
        //if(contiteracions > 10) itefinGS = true;
    
	}
	sim.solid.Ts[x] = sim.PCM.Tsolid[x][M-1];
	/*
	cout << " Tf: " << Tf[x];
	cout << " ||Ts-> ";
	for(int i = M-1; i>=0; i--)
	{
		cout << " - "<< Tsolid[x][i] << " ";
	}
	cout << endl;*/
}


void ChargeTopCV(simulation &sim)
{
	double inct = sim.num.inct;
	double UA = sim.num.UA;
	double Qi = sim.num.Qi;
	double mf = sim.phy.mf;
	double mfc = sim.phy.mfc;
	double mfd = sim.phy.mfd;
	double Tmid = sim.phy.Tmid;
	double Text = sim.phy.Text;
	double Thot = sim.phy.Thot;
	
    double rhofi = sim.topCV.rho;
    double cpfi = sim.topCV.cp;
    double lambdafi =sim.topCV.lambda;
    double mufi = sim.topCV.mu;

    double afitop = + mf*cpfi + (UA)/2.0;
    double bfitop = rhofi*cpfi*sim.topCV.vol/inct - mf*cpfi + (UA)/2.0;
    double cfitop = rhofi*cpfi*sim.topCV.vol/(inct)*sim.topCV.Tant + (UA)*Text - Qi;
    sim.topCV.T = (cfitop - afitop*Thot) / bfitop;
    
}

void ChargeBottomCV(simulation &sim)
{
	double inct = sim.num.inct;
	double UA = sim.num.UA;
	double Qi = sim.num.Qi;
    double rhofi = sim.bottomCV.rho;
    double cpfi = sim.bottomCV.cp;
    double lambdafi = sim.bottomCV.lambda;
    double mufi = sim.bottomCV.mu;
	double mf = sim.phy.mf;
	double mfc = sim.phy.mfc;
	double mfd = sim.phy.mfd;
	double Tmid = sim.phy.Tmid;
	double Text = sim.phy.Text;
	double Thot = sim.phy.Thot;

    double afibottom = mf*cpfi + (UA)/2.0;
    double bfibottom = rhofi*cpfi*sim.bottomCV.vol/inct - mf*cpfi + (UA)/2.0;
    double cfibottom = rhofi*cpfi*sim.bottomCV.vol/(inct)*sim.bottomCV.Tant + (UA)*Text - Qi;
    sim.bottomCV.T = (cfibottom - afibottom*sim.fluid.Tf[0]) / bfibottom;
}

void ChargeMid(simulation &sim)
{
	double Sfi = sim.phy.Sfi;
	double dpi = sim.phy.dpi;
	double voidf = sim.phy.voidf;
	double Vfi = sim.phy.Vfi;
	double incz = sim.phy.incz;
	double NodesPCMbot = sim.phy.NodesPCMbot;
	double voidfPCM = sim.phy.voidfPCM;
	double VfiPCM = sim.phy.VfiPCM;
	double NodesPCMtop = sim.phy.NodesPCMtop;
	double Vbi = sim.phy.Vbi;
	double inct = sim.num.inct;
	int N = sim.num.N;
	double g = sim.num.g;
	double UA = sim.num.UA;
	double Qi = sim.num.Qi;
	double Qbot = sim.num.Qbot;
	double Qtop = sim.num.Qtop;
	double errorGS = sim.num.errorGS;
	double mf = sim.phy.mf;
	double mfc = sim.phy.mfc;
	double mfd = sim.phy.mfd;
	double Tmid = sim.phy.Tmid;
	double Text = sim.phy.Text;
	double Thot = sim.phy.Thot;
	
    for(int i = N-1; i>=0; i--){
    	bool itefinGS = false;
        int contiteracions = 0;
	    while(!itefinGS){
            contiteracions += 1;
	    	double Tsact = sim.solid.Ts[i];
	        // Actualize properties
	        double rhofi = sim.fluid.rhof[i]; 
	        double cpfi = sim.fluid.cpf[i];
	        double lambdafi = sim.fluid.lambdaf[i];
	        double mufi = sim.fluid.muf[i];
	
	        // Fluid state
	        double velfi = mf / (sim.fluid.rhof[i]*Sfi);
	        
 			// Fluid momentum
	        double Rei = rhofi*fabs(velfi)*dpi / (6.0*(1.0-voidf)*mufi);
	        double K = 5.0/Rei + 0.4/pow(Rei, 0.1);
	        sim.fluid.P[i] = ( K * (6.0*(1-voidf)*rhofi*fabs(velfi)*Vfi / (dpi*pow(voidf,3.0))) + rhofi*g ) * incz + sim.fluid.P[i+1];
	        
	        // Fluid energy
	        if(i == 0) Qi = Qbot;
	        else if(i == N-1) Qi = Qtop;
	        else Qi = 0.0;
	        //CalcLosses(Tf[i], UA);
	         
	        double Tfn = 0.0;
	        if(i == N-1) Tfn = sim.topCV.T;
	        else Tfn = sim.fluid.Tf[i+1];
	         
	        double afi = 0.0;
	        double bfi = 0.0;
	        double cfi = 0.0;
	        double abi = 0.0;
	        double bbi = 0.0;
	        
	        if (i < NodesPCMbot) // PCM bottom
			{
                // Fluid
	        	calcAlpha(i, sim.PCMtype, voidfPCM, sim);
                afi = mf*cpfi + (UA)/2.0;
                bfi = rhofi*cpfi*VfiPCM/inct - mf*cpfi + sim.fluid.alpha[i] + (UA)/2.0;
                cfi = rhofi*cpfi*VfiPCM/(inct)*sim.fluid.Tantf[i] + sim.fluid.alpha[i]*sim.solid.Ts[i] + (UA)*Text - Qi;
                sim.fluid.Tf[i] = (cfi - afi*Tfn) / bfi;
	        	// Solid
	        	if(sim.Modeltype == "1D")
	        	{
	        		abi = sim.PCM.Mass1D[i]/inct;
                	bbi = sim.fluid.alpha[i]*(sim.fluid.Tf[i]-sim.solid.Ts[i]) + sim.PCM.Mass1Dant[i]*sim.PCM.h1Dant[i]/inct;
                	sim.PCM.h1D[i] = bbi / abi;
                	sim.solid.Ts[i] = calcentalpia(sim.PCM.h1D[i], sim);
				}
				else if(sim.Modeltype == "2D")
				{
					calcAlphaPCM(i, sim);
	        		ChargeSolidEntalpia(i, sim);
				}
			}
            else if (i > (N - NodesPCMtop)) // PCM top
            {
                // Fluid
	        	calcAlpha(i, sim.PCMtype, voidfPCM, sim);
                afi = mf*cpfi + (UA)/2.0;
                bfi = rhofi*cpfi*VfiPCM/inct - mf*cpfi + sim.fluid.alpha[i] + (UA)/2.0;
                cfi = rhofi*cpfi*VfiPCM/(inct)*sim.fluid.Tantf[i] + sim.fluid.alpha[i]*sim.solid.Ts[i] + (UA)*Text - Qi;
                sim.fluid.Tf[i] = (cfi - afi*Tfn) / bfi;
                // Solid
                if(sim.Modeltype == "1D")
	        	{
	        		abi = sim.PCM.Mass1D[i]/inct;
                	bbi = sim.fluid.alpha[i]*(sim.fluid.Tf[i]-sim.solid.Ts[i]) + sim.PCM.Mass1Dant[i]*sim.PCM.h1Dant[i]/inct;
                	sim.PCM.h1D[i] = bbi / abi;
                	sim.solid.Ts[i] = calcentalpia(sim.PCM.h1D[i], sim);
				}
				else if(sim.Modeltype == "2D")
				{
					calcAlphaPCM(i, sim);
	        		ChargeSolidEntalpia(i, sim);
				}
            }
	        else // Filler
	        {
                // Fluid
				calcAlpha(i, sim.fillerType, voidf, sim);	
                afi = mf*cpfi + (UA)/2.0;
                bfi = rhofi*cpfi*Vfi/inct - mf*cpfi + sim.fluid.alpha[i] + (UA)/2.0;
                cfi = rhofi*cpfi*Vfi/(inct)*sim.fluid.Tantf[i] + sim.fluid.alpha[i]*sim.solid.Ts[i] + (UA)*Text - Qi;
                sim.fluid.Tf[i] = (cfi - afi*Tfn) / bfi;
                
                // Solid
                abi = sim.solid.rhos[i]*Vbi*sim.solid.cps[i]/inct + sim.fluid.alpha[i];
                bbi = sim.solid.rhos[i]*Vbi*sim.solid.cps[i]/inct*sim.solid.Tants[i] + sim.fluid.alpha[i]*sim.fluid.Tf[i];
				sim.solid.Ts[i] = bbi / abi;	
			}
	        
	        if(errorGS > fabs(sim.solid.Ts[i]-Tsact)) itefinGS = true;
            if(contiteracions > 10) itefinGS = true;
            if (isnan(sim.solid.Ts[i])) cout <<"i: " << i << endl;
            //break;
	    }
    }
}

void DischargeBottomCV(simulation &sim)
{
	double inct = sim.num.inct;
	double g = sim.num.g;
	double UA = sim.num.UA;
	double Qi = sim.num.Qi;
	double Qbot = sim.num.Qbot;
	double Qtop = sim.num.Qtop;
	double mf = sim.phy.mf;
	double mfc = sim.phy.mfc;
	double mfd = sim.phy.mfd;
	double Tmid = sim.phy.Tmid;
	double Text = sim.phy.Text;
	double Thot = sim.phy.Thot;
	double Tcold = sim.phy.Tcold;
	//double errorGS = sim.num.errorGS;
    double rhofi = sim.bottomCV.rho;
    double cpfi = sim.bottomCV.cp;
    double lambdafi = sim.bottomCV.lambda;
    double mufi = sim.bottomCV.mu;

    double afibottom = rhofi*cpfi*sim.bottomCV.vol/inct + mf*cpfi + (UA)/2.0;
    double bfibottom = - mf*cpfi + (UA)/2.0;
    double cfibottom = rhofi*cpfi*sim.bottomCV.vol/(inct)*sim.bottomCV.Tant + (UA)*Text - Qi;
    sim.bottomCV.T = (cfibottom - bfibottom*Tcold) / afibottom;
}

void DischargeTopCV(simulation &sim)
{
	int N = sim.num.N;
	double inct = sim.num.inct;
	double g = sim.num.g;
	double UA = sim.num.UA;
	double Qi = sim.num.Qi;
	double Qbot = sim.num.Qbot;
	double Qtop = sim.num.Qtop;
	double mf = sim.phy.mf;
	double mfc = sim.phy.mfc;
	double mfd = sim.phy.mfd;
	double Tmid = sim.phy.Tmid;
	double Text = sim.phy.Text;
	double Thot = sim.phy.Thot;
	double Tcold = sim.phy.Tcold;
	//double errorGS = sim.num.errorGS;
    double rhofi = sim.topCV.rho; 
    double cpfi = sim.topCV.cp;
    double lambdafi = sim.topCV.lambda;
    double mufi = sim.topCV.mu;

    double afitop = rhofi*cpfi*sim.topCV.vol/inct + mf*cpfi + (UA)/2.0;
    double bfitop = - mf*cpfi + (UA)/2.0;
    double cfitop = rhofi*cpfi*sim.topCV.vol/(inct)*sim.topCV.Tant + (UA)*Text - Qi;
    sim.topCV.T = (cfitop - bfitop*sim.fluid.Tf[N-1]) / afitop;
}

void DischargeMid(simulation &sim)
{
	double Sfi = sim.phy.Sfi;
	double dpi = sim.phy.dpi;
	double voidf = sim.phy.voidf;
	double Vfi = sim.phy.Vfi;
	double incz = sim.phy.incz;
	double NodesPCMbot = sim.phy.NodesPCMbot;
	double voidfPCM = sim.phy.voidfPCM;
	double VfiPCM = sim.phy.VfiPCM;
	double NodesPCMtop = sim.phy.NodesPCMtop;
	double Vbi = sim.phy.Vbi;	
	double inct = sim.num.inct;
	int N = sim.num.N;
	double g = sim.num.g;
	double UA = sim.num.UA;
	double Qi = sim.num.Qi;
	double Qbot = sim.num.Qbot;
	double Qtop = sim.num.Qtop;
	double errorGS = sim.num.errorGS;
	double mf = sim.phy.mf;
	double mfc = sim.phy.mfc;
	double mfd = sim.phy.mfd;
	double Tmid = sim.phy.Tmid;
	double Text = sim.phy.Text;
	double Thot = sim.phy.Thot;
	double Tcold = sim.phy.Tcold;
	//cout << "Sfi: " << Sfi << " dpi: " << dpi << " voidf: " << voidf << " Vfi: " << Vfi << " incz: " << incz << " Vbi: " << Vbi << endl;
	
    for(int i = 0; i<N; i++){
    	bool itefinGS = false;
        int contiteracions = 0;
    	while(!itefinGS){
            contiteracions += 1;
    		double Tsact = sim.solid.Ts[i];
    		double Tfact = sim.fluid.Tf[i];
	        // Actualize properties
	        double rhofi = sim.fluid.rhof[i]; 
	        double cpfi = sim.fluid.cpf[i];
	        double lambdafi = sim.fluid.lambdaf[i];
	        double mufi = sim.fluid.muf[i];
	        // Fluid state
	        double velfi = mf / (sim.fluid.rhof[i]*Sfi);
	        
		    // Fluid momentum
	        double Rei = rhofi*velfi*dpi / (6.0*(1.0-voidf)*mufi);
	        double K = 5.0/Rei + 0.4/pow(Rei, 0.1);
	        sim.fluid.P[i] = ( K * (6.0*(1-voidf)*rhofi*fabs(velfi)*Vfi / (dpi*pow(voidf,3.0))) + rhofi*g ) * incz + sim.fluid.P[i+1];

	        // Fluid energy
	        if(i == 0) Qi = Qbot;
	        else if(i == N) Qi = Qtop;
	        else Qi = 0.0;
	        //CalcLosses(Tf[i], UA);
	        
	        double Tsud = 0.0;
	        if(i == 0) Tsud = sim.bottomCV.T;
	        else Tsud = sim.fluid.Tf[i-1];
	
	        double afi = 0.0;
	        double bfi = 0.0;
	        double cfi = 0.0;
	        double abi = 0.0;
	        double bbi = 0.0;
	
	        
	        if (i < NodesPCMbot)
	        {	
                // Fluid
	        	calcAlpha(i, sim.PCMtype, voidfPCM, sim);
                afi = rhofi*cpfi*VfiPCM/inct + mf*cpfi + sim.fluid.alpha[i] + (UA)/2.0;
                bfi = - mf*cpfi + (UA)/2.0;
                cfi = rhofi*cpfi*VfiPCM/(inct)*sim.fluid.Tantf[i] + sim.fluid.alpha[i]*sim.solid.Ts[i] + (UA)*Text - Qi;
                sim.fluid.Tf[i] = (cfi - bfi*Tsud) / afi;
                //Tf[i] = Tcold;
                //cout <<  "i: " << i << " Tf: " << Tf[i];
        		// Solid
                if(sim.Modeltype == "1D")
	        	{
	        		abi = sim.PCM.Mass1D[i]/inct;
                	bbi = sim.fluid.alpha[i]*(sim.fluid.Tf[i]-sim.solid.Ts[i]) + sim.PCM.Mass1Dant[i]*sim.PCM.h1Dant[i]/inct;
                	sim.PCM.h1D[i] = bbi / abi;
                    //cout <<"alphai: " << alpha[i]<<" Mass1D: " << Mass1D[i] << " abi: " << abi << " bbi: " << bbi << " h1d: "<<h1D[i] << endl;
                	sim.solid.Ts[i] = calcentalpia(sim.PCM.h1D[i], sim);
                	//cout << " Ts: " << Ts[i] << endl;
				}
				else if(sim.Modeltype == "2D")
				{
					calcAlphaPCM(i, sim);
	        		ChargeSolidEntalpia(i, sim);
				}
				//cout <<  " i: " << i << " Tf: " << Tf[i] << " Ts: " << Ts[i] << endl;
			}
            else if (i > (N - NodesPCMtop))
            {
                // Fluid
	        	calcAlpha(i, sim.PCMtype, voidfPCM, sim);
                afi = rhofi*cpfi*VfiPCM/inct + mf*cpfi + sim.fluid.alpha[i] + (UA)/2.0;
                bfi = - mf*cpfi + (UA)/2.0;
                cfi = rhofi*cpfi*VfiPCM/(inct)*sim.fluid.Tantf[i] + sim.fluid.alpha[i]*sim.solid.Ts[i] + (UA)*Text - Qi;
                sim.fluid.Tf[i] = (cfi - bfi*Tsud) / afi;
        		// Solid
               	if(sim.Modeltype == "1D")
	        	{
	        		abi = sim.PCM.Mass1D[i]/inct;
                	bbi = sim.fluid.alpha[i]*(sim.fluid.Tf[i]-sim.solid.Ts[i]) + sim.PCM.Mass1Dant[i]*sim.PCM.h1Dant[i]/inct;
                	sim.PCM.h1D[i] = bbi / abi;
                	sim.solid.Ts[i] = calcentalpia(sim.PCM.h1D[i], sim);
				}
				else if(sim.Modeltype == "2D")
				{
					calcAlphaPCM(i, sim);
	        		ChargeSolidEntalpia(i, sim);
				}
            }
	        else
			{
                // Fluid
				calcAlpha(i, sim.fillerType, voidf, sim);
                afi = rhofi*cpfi*Vfi/inct + mf*cpfi + sim.fluid.alpha[i] + (UA)/2.0;
                bfi = - mf*cpfi + (UA)/2.0;
                cfi = rhofi*cpfi*Vfi/(inct)*sim.fluid.Tantf[i] + sim.fluid.alpha[i]*sim.solid.Ts[i] + (UA)*Text - Qi;
                sim.fluid.Tf[i] = (cfi - bfi*Tsud) / afi;
        
                // Solid
                abi = sim.solid.rhos[i]*Vbi*sim.solid.cps[i]/inct + sim.fluid.alpha[i];
                bbi = sim.solid.rhos[i]*Vbi*sim.solid.cps[i]/inct*sim.solid.Tants[i] + sim.fluid.alpha[i]*sim.fluid.Tf[i];
				sim.solid.Ts[i] = bbi / abi;	
			}
			       
	        if(errorGS > fabs(sim.solid.Ts[i] - Tsact) and errorGS > fabs(sim.fluid.Tf[i] - Tfact)) itefinGS = true;
	        //cout << "Contite: " << contiteracions << " Ts[i]: "<< Ts[i] << " Tsact: " << Tsact << " Tf[i]: " << Tf[i] << " Tfact: " << Tfact << endl;
            //if(contiteracions > 10) itefinGS = true;
            //if (isnan(Ts[i])) cout <<"i: " << i << endl;
            //break;
	    }
    }
}

void Print(int ite, simulation &sim)
{
	double t = sim.num.t;
	//if(ite%1000==0)cout << endl << "Temps: " << t/3660.0 << " [h] " << "estat: " << estat << " htop: " << topCV.incz << endl;	
	if(ite%10==0)cout << endl << "Temps: " << t/3660.0 << " [h] " << "sim.estat: " << sim.estat << " htop: " << sim.topCV.incz << endl;	
}

void validation(simulation &sim)
{
	int N = sim.num.N;
	int M = sim.num.M;
    std::ifstream myfile ("Tboundary.txt");
    std::string myline;
    int i = 0;
    if ( myfile.is_open() ) {
        while ( myfile ) {
            std::getline (myfile, myline);
            std::cout << myline << ": " << i << '\n';
            sim.fluid.Tf[i] = stod(myline);
            sim.fluid.Tantf[i] = stod(myline);
            sim.solid.Ts[i] = stod(myline);
            sim.solid.Tants[i] = stod(myline);
            // Initialize bolitas
            for(int j = 0; j<M; j++)
            {
            	sim.PCM.Tsolid[i][j] = stod(myline);
            	sim.PCM.Tantsolid[i][j] = stod(myline);
			}
            //
            cout << "Tf: " << sim.fluid.Tf[i] << " Ts: " << sim.solid.Ts[i] << endl;
            i++;
            if(i>N-1) break;
        }
    }
}


void FillParams(std::unordered_map<std::string, double> &atributos, string name)
{
	std::ifstream archivo("Input/" + name + ".txt");
	if (archivo.is_open()) {
        std::string nombre;
        double valor;

        // Leer pares nombre-valor desde el archivo y asignarlos a los atributos
        while (archivo >> nombre >> valor) {
            atributos[nombre] = valor;
        }

        // Cerrar el archivo después de leer
        archivo.close();
    }else {
        std::cerr << "Error opening the file: " << "input.txt"<< std::endl;
    }
}

int main(){
	// Entrem parametres per inicialitzar l'objecte simulation
	std::unordered_map<std::string, double> atributos;
	FillParams(atributos, "input");

	int N = static_cast<int>(atributos["N"]);
    int M = static_cast<int>(atributos["M"]);
    double Tinicial = atributos["Tinicial"];
    double htop = atributos["htop"];
    double hbottom = atributos["hbottom"];
    double D = atributos["D"];
    
	simulation sim(N, M, Tinicial, htop, hbottom, D);
	
    setupNumerical(sim, atributos);
	setupPhysical(sim, atributos);
	
	cout << "L: " << sim.phy.L << endl;
	cout << "D: " << sim.phy.D  << endl;
    //cout << "Vi: " << Vi << " Vbi: " << Vbi << " Vfi: " << Vfi << " Sfi: " << Sfi << " dpi= " << dpi <<endl;

	double inct = sim.num.inct;
	double tmax = sim.num.tmax;
	double  t = sim.num.t;
	int ite = sim.num.ite;
	double mf = sim.phy.mf;
	double mfc = sim.phy.mfc;
	double mfd = sim.phy.mfd;
	
    //Start temperature ma
    //validation();
    for(auto x: sim.fluid.Tf) cout << x << endl;
    cout << "********SOLID********" << endl;
    for(auto x: sim.solid.Ts) cout << x << endl;

    //Initialize and actualize Top and Bottom CV
    sim.topCV.updateProps();
    sim.bottomCV.updateProps();

    // Initialize Files and alpha mid
    InitializeFiles();
    CalculateAlphaInit(sim);

    //Initialize and actualize fluid and solid bed properties
    calculateProperties(sim);
    actualizeProperties(sim);
    
	initentalpia(sim);
	initmass(sim);
	//Borrar
	/*
	for(int i = 0; i<N; i++)
	{
		for(int j = 0; j<M; j++)
		{
			cout << "Tsolid: " << Tsolid[i][j] << " hsolid: " << h[i][j] << " hant: " << hant[i][j] << " ";
		}
		cout << endl;
	}
	*/
    //
	// Temporal discretization
	while(t < tmax){
        ite++;
        sim.num.ite = ite;
		
		//cout << "inct: " << inct << " tmax: " << tmax << " t: " << t << " ite: " << ite << endl;
        // Print actual state
        Print(ite, sim);

        // Choose if the cycle is charge or discharge
        ChooseChargeDischarge(sim); 
        

        // Calculate new properties
        calculateProperties(sim);
        sim.topCV.setT(sim.topCV.T);
        sim.bottomCV.setT(sim.bottomCV.T);

        // Charging
        if(sim.estat == 1)
        {
            mf = -mfc;
            sim.phy.mf = mf;
            ChargeTopCV(sim);
            ChargeMid(sim);
            ChargeBottomCV(sim);
        }
        // Discharging
        else if(sim.estat == 2)
        {
            mf = mfd;
            sim.phy.mf = mf;
            DischargeBottomCV(sim);
            DischargeMid(sim);
            DischargeTopCV(sim);
        } 

        // Apply continuity equation
        //calcContinuity();
        
		// Post Process
        PostProcess(sim);
        
        // Actualize solid and fluid properties
        actualizeProperties(sim);
        sim.topCV.updateProps();
        sim.bottomCV.updateProps();

        // Time-step and cycle actualization
        t += inct;
        sim.num.t = t;
		if(sim.estat != sim.estatanterior) sim.estatanterior = sim.estat;
		
		
		// *********Proves solid PCM*********
		//for(int i = 0; i<N; i++) ChargeSolid(i);
		//actualizeProperties();
		// **********************************
	}	
	
	
}
