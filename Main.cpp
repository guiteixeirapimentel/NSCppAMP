#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cassert>
#include <string>

#include <amp.h>
#include <amp_math.h>

#include <Windows.h>

struct StateData
{
	float u;
	float v;
	float P;
	float T;

	
	StateData(float u, float v, float P, float T) :u(u), v(v), P(P), T(T) {};
	StateData() :u(0.0f), v(0.0f), P(0.0f), T(0.0f) {};
};

//State CalcState(const State& U);
//State CalcU(const State& estado);

void CalcExtrapolationsBC(concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn,
	concurrency::array_view<float, 2> Pn, concurrency::array_view<float, 2> Tn);

float CalcDensity(float P, float R, float T) restrict(amp);
float CalcViscSutherland(float mu0, float T0, float T) restrict(amp);
float CalcEt(float rho, float Cv, float T, float Vquad) restrict(amp);
float Calckk(float mu, float Cp, float prandtlNumber) restrict(amp);

float CalcDensity(float P, float R, float T);
float CalcViscSutherland(float mu0, float T0, float T);
float CalcEt(float rho, float Cv, float T, float Vquad);
float Calck(float mu, float Cp, float prandtlNumber);

float GetE1Predictor(int i, int j, concurrency::array_view<float, 2> rhon,
	concurrency::array_view<float, 2> un) restrict(amp);

float GetE2Predictor(int i, int j, concurrency::array_view<float, 2> rhon,
	concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> Pn, concurrency::array_view<float, 2> mun) restrict(amp);

float GetE3Predictor(int i, int j, concurrency::array_view<float, 2> rhon,
	concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun)restrict(amp);

float GetE4Predictor(int i, int j, concurrency::array_view<float, 2> Etn,
	concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun,
	concurrency::array_view<float, 2> Tn, concurrency::array_view<float, 2> Kn, concurrency::array_view<float, 2> Pn)restrict(amp);

float GetF1Predictor(int i, int j, concurrency::array_view<float, 2> rhon, 
	concurrency::array_view<float, 2> vn)restrict(amp);

float GetF2Predictor(int i, int j, concurrency::array_view<float, 2> rhon,
	concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun) restrict(amp);

float GetF3Predictor(int i, int j, concurrency::array_view<float, 2> rhon, concurrency::array_view<float, 2> un,
	concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun,
	concurrency::array_view<float, 2> Pn) restrict(amp);

float GetF4Predictor(int i, int j, concurrency::array_view<float, 2> Etn,
	concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun,
	concurrency::array_view<float, 2> Tn, concurrency::array_view<float, 2> Kn, concurrency::array_view<float, 2> Pn) restrict(amp);

float CalcTauXXPredictor(int i, int j, concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn,
	concurrency::array_view<float, 2> mun)restrict(amp);

float CalcTauYYPredictor(int i, int j, concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn,
	concurrency::array_view<float, 2> mun) restrict(amp);

float CalcTauXYPredictorE(int i, int j, concurrency::array_view<float, 2> un,
	concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun) restrict(amp);

float CalcTauXYPredictorF(int i, int j, concurrency::array_view<float, 2> un,
	concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun) restrict(amp);

float CalcQxPredictor(int i, int j, concurrency::array_view<float, 2> Tn,
	concurrency::array_view<float, 2> Kn) restrict(amp);

float CalcQyPredictor(int i, int j, concurrency::array_view<float, 2> Tn,
	concurrency::array_view<float, 2> Kn) restrict(amp);

float GetE1Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn,
	concurrency::array_view<float, 2> uPredn) restrict(amp);

float GetE2Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn,
	concurrency::array_view<float, 2> uPredn, concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> PPredn, concurrency::array_view<float, 2> muPredn) restrict(amp);

float GetE3Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn,
	concurrency::array_view<float, 2> uPredn, concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn)restrict(amp);

float GetE4Corrector(int i, int j, concurrency::array_view<float, 2> EtPredn,
	concurrency::array_view<float, 2> uPredn, concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn,
	concurrency::array_view<float, 2> TPredn, concurrency::array_view<float, 2> KPredn, concurrency::array_view<float, 2> PPredn)restrict(amp);


float GetF1Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn,
	concurrency::array_view<float, 2> vPredn)restrict(amp);

float GetF2Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn,
	concurrency::array_view<float, 2> uPredn, concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn) restrict(amp);

float GetF3Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn, concurrency::array_view<float, 2> uPredn,
	concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn,
	concurrency::array_view<float, 2> PPredn) restrict(amp);

float GetF4Corrector(int i, int j, concurrency::array_view<float, 2> EtPredn,
	concurrency::array_view<float, 2> uPredn, concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn,
	concurrency::array_view<float, 2> TPredn, concurrency::array_view<float, 2> KPredn, concurrency::array_view<float, 2> PPredn) restrict(amp);

float CalcTauXXCorrector(int i, int j, concurrency::array_view<float, 2> uPredn, 
	concurrency::array_view<float, 2> vPredn,
	concurrency::array_view<float, 2> muPredn) restrict(amp);

float CalcTauYYCorrector(int i, int j, concurrency::array_view<float, 2> uPredn,
	concurrency::array_view<float, 2> vPredn,
	concurrency::array_view<float, 2> muPredn) restrict(amp);

float CalcTauXYCorrectorE(int i, int j, concurrency::array_view<float, 2> uPredn,
	concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn) restrict(amp);

float CalcTauXYCorrectorF(int i, int j, concurrency::array_view<float, 2> uPredn,
	concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn) restrict(amp);

float CalcQxCorrector(int i, int j, concurrency::array_view<float, 2> TPredn,
	concurrency::array_view<float, 2> KPredn) restrict(amp);

float CalcQyCorrector(int i, int j, concurrency::array_view<float, 2> TPredn,
	concurrency::array_view<float, 2> KPredn) restrict(amp);

float CalcSrho(int i, int j, concurrency::array_view<float, 2> rhon, 
	float xcoef, float ycoef) restrict(amp);

float CalcSrhou(int i, int j, concurrency::array_view<float, 2> U2n,
	float xcoef, float ycoef) restrict(amp);

float CalcSrhov(int i, int j, concurrency::array_view<float, 2> U3n, 
	float xcoef, float ycoef) restrict(amp);

float CalcSEt(int i, int j, concurrency::array_view<float, 2> Etn,
	float xcoef, float ycoef) restrict(amp);

#define machRefD float(0.15f/340.0f)
#define plateLengthD float(0.05f)
#define fieldLengthD float(plateLengthD * 3.0f)
#define aSeaLevelD float(340.28f)
#define velRefD float(machRefD * aSeaLevelD)
#define alfaD float(0.0f* M_PI / 180.0f)
#define pressureSeaLevelD float(101325.0f)
#define temperatureSeaLevelD float(288.16f)
#define gamaD float(1.4f) // Cp/Cv
#define prandtlNumberD float(0.71f)
#define muRefD float(1.7894e-5F * 2.0)// * 0.2f)
#define tempRefD float(temperatureSeaLevelD)
#define RD float(287.0f)
#define CvD float(RD / (gamaD - 1.0f))
#define CpD float(gamaD * CvD)
#define RHOrefD float(pressureSeaLevelD / (RD*tempRefD))
#define NPOINTSXD int(1000)
#define NPOINTSYD int(1000)//1000;
#define NITERATIONSD int(10000000)
#define LHORID float(fieldLengthD)
#define ReLD float(RHOrefD * velRefD*plateLengthD / muRefD)
#define LVERTD float(3.0f*plateLengthD)
#define dxD float(fieldLengthD / (NPOINTSXD - 1))
#define dyD float(LVERTD / (NPOINTSYD - 1))

#define NPOINTSPXD int(plateLengthD / dxD)


const float velURef = velRefD * cos(alfaD);
const float velVRef = velRefD * sin(alfaD);

#define CxD (0.1f)
#define CyD (0.1f)

#define saveEvery int(500) // iterations

#define SALVARASCII (true)

int main()
{
	std::cout << "Re: " << ReLD << std::endl;
	//std::vector<float> U1n; // rho
	//std::vector<float> U2n; // rho*u
	//std::vector<float> U3n; // rho*v
	//std::vector<float> U4n; // Et

	std::vector<float> un;
	std::vector<float> vn;
	std::vector<float> Pn;
	std::vector<float> Tn;

	//std::vector<float> upn;
	//std::vector<float> vpn;
	//std::vector<float> Ppn;
	//std::vector<float> Tpn;

	//std::vector<float> U1pn; // rho
	//std::vector<float> U2pn; // rho*u
	//std::vector<float> U3pn; // rho*v
	//std::vector<float> U4pn; // Et

	//std::vector<float> dU1pdt;
	//std::vector<float> dU2pdt;
	//std::vector<float> dU3pdt;
	//std::vector<float> dU4pdt;
	//
	//std::vector<float> dU1cdt;
	//std::vector<float> dU2cdt;
	//std::vector<float> dU3cdt;
	//std::vector<float> dU4cdt;

	std::vector<float> mun;
	std::vector<float> kn;

	mun.resize(NPOINTSYD* NPOINTSXD, CalcViscSutherland(muRefD, tempRefD, tempRefD));
	kn.resize(NPOINTSYD* NPOINTSXD, Calck(muRefD, CpD, prandtlNumberD));

	un.resize(NPOINTSYD * NPOINTSXD);
	vn.resize(NPOINTSYD * NPOINTSXD);
	Pn.resize(NPOINTSYD * NPOINTSXD);
	Tn.resize(NPOINTSYD * NPOINTSXD);

	// SET INITIAL AND BOUNDARY CONDITIONS

	char b = 0;
	std::cout << "Reiniciar do zero? (y/n): ";
	std::cin >> b;
	std::cout << "\n";

	for (size_t i = 0; i < NPOINTSYD; i++)
		for (size_t j = 0; j < NPOINTSXD; j++)
		{
			un[j + (i * NPOINTSXD)] = velURef;
			vn[j + (i * NPOINTSXD)] = velVRef;
			Pn[j + (i * NPOINTSXD)] = pressureSeaLevelD;
			Tn[j + (i * NPOINTSXD)] = tempRefD;
		}

	unsigned int itInicial = 0;

	if (b == 'n')
	{
		std::cout << "Digite a ultima iteracao calculada: ";
		std::cin >> itInicial;
		std::cout << "\n";


		std::string nomeArq = "out/state";
		nomeArq += std::to_string(itInicial);
		nomeArq += ".dat";

		itInicial += 1;

		std::vector<StateData> estado;
		estado.resize(NPOINTSXD * NPOINTSYD);

		std::ifstream arqConv(nomeArq, std::ios::binary);
		char buf[500];
		int width, height;
		arqConv.getline(buf, 500);
		arqConv >> b; //%
		arqConv >> b; //D
		arqConv >> b; //A
		arqConv >> b; //T
		arqConv >> b; //#

		arqConv.read(reinterpret_cast<char*>(&width), sizeof(int));
		arqConv >> b;
		arqConv.read(reinterpret_cast<char*>(&height), sizeof(int));
		arqConv.get();
		
		arqConv.read(reinterpret_cast<char*>(estado.data()), sizeof(StateData)*height*width);
		
		arqConv.close();

		std::cout << "Lido arquivo com campo de largura: " << width << "\nAltura: " << height << std::endl;

		for (size_t idx = 0; idx < NPOINTSXD * NPOINTSYD; idx++)
		{
			mun[idx] = CalcViscSutherland(muRefD, tempRefD, tempRefD);
			kn[idx] = Calck(mun[idx], CpD, prandtlNumberD);
		}
	}

	const float PINF = pressureSeaLevelD;
	const float TINF = temperatureSeaLevelD;
	const float TWALL = tempRefD;
	const float UINF = velURef;
	const float VINF = velVRef;

	// CASE 1 (LEADING EDGES AND TRAILING EDGE)
	Pn[(NPOINTSXD/2) - (NPOINTSPXD/2) + (((NPOINTSYD/2) - 1) * NPOINTSXD)] = PINF;
	Tn[(NPOINTSXD/2) - (NPOINTSPXD/2) + (((NPOINTSYD/2) - 1) * NPOINTSXD)] = TINF;
	un[(NPOINTSXD/2) - (NPOINTSPXD/2) + (((NPOINTSYD/2) - 1) * NPOINTSXD)] = 0.0f;
	vn[(NPOINTSXD/2) - (NPOINTSPXD/2) + (((NPOINTSYD/2) - 1) * NPOINTSXD)] = 0.0f;
				
	Pn[(NPOINTSXD / 2) - (NPOINTSPXD / 2) + (((NPOINTSYD / 2)) * NPOINTSXD)] = PINF;
	Tn[(NPOINTSXD / 2) - (NPOINTSPXD / 2) + (((NPOINTSYD / 2)) * NPOINTSXD)] = TINF;
	un[(NPOINTSXD / 2) - (NPOINTSPXD / 2) + (((NPOINTSYD / 2)) * NPOINTSXD)] = 0.0f;
	vn[(NPOINTSXD / 2) - (NPOINTSPXD / 2) + (((NPOINTSYD / 2)) * NPOINTSXD)] = 0.0f;
				
	Pn[(NPOINTSXD / 2) + (NPOINTSPXD / 2) + (((NPOINTSYD / 2) - 1) * NPOINTSXD)] = PINF;
	Tn[(NPOINTSXD / 2) + (NPOINTSPXD / 2) + (((NPOINTSYD / 2) - 1) * NPOINTSXD)] = TINF;
	un[(NPOINTSXD / 2) + (NPOINTSPXD / 2) + (((NPOINTSYD / 2) - 1) * NPOINTSXD)] = 0.0f;
	vn[(NPOINTSXD / 2) + (NPOINTSPXD / 2) + (((NPOINTSYD / 2) - 1) * NPOINTSXD)] = 0.0f;
				
	Pn[(NPOINTSXD / 2) + (NPOINTSPXD / 2) + (((NPOINTSYD / 2)) * NPOINTSXD)] = PINF;
	Tn[(NPOINTSXD / 2) + (NPOINTSPXD / 2) + (((NPOINTSYD / 2)) * NPOINTSXD)] = TINF;
	un[(NPOINTSXD / 2) + (NPOINTSPXD / 2) + (((NPOINTSYD / 2)) * NPOINTSXD)] = 0.0f;
	vn[(NPOINTSXD / 2) + (NPOINTSPXD / 2) + (((NPOINTSYD / 2)) * NPOINTSXD)] = 0.0f;

	// CASE 2 (INFLOW/LEFT.TOP.BOTTOM)

	for (unsigned int i = 0; i < NPOINTSYD; i++)
	{
		Pn[0 + (i*NPOINTSXD)] = PINF;
		Tn[0 + (i*NPOINTSXD)] = TINF;
		un[0 + (i*NPOINTSXD)] = UINF;
		vn[0 + (i*NPOINTSXD)] = VINF;
	}

	for (unsigned int j = 0; j < NPOINTSXD; j++)
	{
		Pn[j + (0 * NPOINTSXD)] = PINF;
		Tn[j + (0 * NPOINTSXD)] = TINF;
		un[j + (0 * NPOINTSXD)] = UINF;
		vn[j + (0 * NPOINTSXD)] = VINF;

		Pn[j + ((NPOINTSYD - 1) * NPOINTSXD)] = PINF;
		Tn[j + ((NPOINTSYD - 1) * NPOINTSXD)] = TINF;
		un[j + ((NPOINTSYD - 1) * NPOINTSXD)] = UINF;
		vn[j + ((NPOINTSYD - 1) * NPOINTSXD)] = VINF;
	}

	// CASE 3 (SURFACE EXCEPT LEADING EDGE)

	for (unsigned int j = (NPOINTSXD/2)- (NPOINTSPXD / 2); j < (NPOINTSXD/2)+(NPOINTSPXD/2); j++)
	{
		Pn[j + (((NPOINTSYD / 2) - 1) * NPOINTSXD)] = PINF;
		Tn[j + (((NPOINTSYD / 2) - 1) * NPOINTSXD)] = TWALL;
		un[j + (((NPOINTSYD / 2) - 1) * NPOINTSXD)] = 0.0f;
		vn[j + (((NPOINTSYD/2) - 1)*NPOINTSXD)] = 0.0f;

		Pn[j + (((NPOINTSYD / 2)) * NPOINTSXD)] = PINF;
		Tn[j + (((NPOINTSYD / 2)) * NPOINTSXD)] = TWALL;
		un[j + (((NPOINTSYD / 2)) * NPOINTSXD)] = 0.0f;
		vn[j + (((NPOINTSYD / 2)) * NPOINTSXD)] = 0.0f;
	}

	concurrency::array_view<float, 2> U1nGPU(NPOINTSYD, NPOINTSXD); // rho
	concurrency::array_view<float, 2> U2nGPU(NPOINTSYD, NPOINTSXD); // rho*u
	concurrency::array_view<float, 2> U3nGPU(NPOINTSYD, NPOINTSXD); // rho*v
	concurrency::array_view<float, 2> U4nGPU(NPOINTSYD, NPOINTSXD); // Et

	concurrency::array_view<float, 2> unGPU(NPOINTSYD, NPOINTSXD, un.data());
	concurrency::array_view<float, 2> vnGPU(NPOINTSYD, NPOINTSXD, vn.data());
	concurrency::array_view<float, 2> PnGPU(NPOINTSYD, NPOINTSXD, Pn.data());
	concurrency::array_view<float, 2> TnGPU(NPOINTSYD, NPOINTSXD, Tn.data());

	concurrency::array_view<float, 2> upnGPU(NPOINTSYD, NPOINTSXD);
	concurrency::array_view<float, 2> vpnGPU(NPOINTSYD, NPOINTSXD);
	concurrency::array_view<float, 2> PpnGPU(NPOINTSYD, NPOINTSXD);
	concurrency::array_view<float, 2> TpnGPU(NPOINTSYD, NPOINTSXD);

	concurrency::array_view<float, 2> U1pnGPU(NPOINTSYD, NPOINTSXD); // rho
	concurrency::array_view<float, 2> U2pnGPU(NPOINTSYD, NPOINTSXD); // rho*u
	concurrency::array_view<float, 2> U3pnGPU(NPOINTSYD, NPOINTSXD); // rho*v
	concurrency::array_view<float, 2> U4pnGPU(NPOINTSYD, NPOINTSXD); // Et

	concurrency::array_view<float, 2> dU1pdtGPU(NPOINTSYD, NPOINTSXD);
	concurrency::array_view<float, 2> dU2pdtGPU(NPOINTSYD, NPOINTSXD);
	concurrency::array_view<float, 2> dU3pdtGPU(NPOINTSYD, NPOINTSXD);
	concurrency::array_view<float, 2> dU4pdtGPU(NPOINTSYD, NPOINTSXD);
											   
	concurrency::array_view<float, 2> dU1cdtGPU(NPOINTSYD, NPOINTSXD);
	concurrency::array_view<float, 2> dU2cdtGPU(NPOINTSYD, NPOINTSXD);
	concurrency::array_view<float, 2> dU3cdtGPU(NPOINTSYD, NPOINTSXD);
	concurrency::array_view<float, 2> dU4cdtGPU(NPOINTSYD, NPOINTSXD);

	concurrency::array_view<float, 2> munGPU(NPOINTSYD, NPOINTSXD, mun.data());
	concurrency::array_view<float, 2> knGPU(NPOINTSYD, NPOINTSXD, kn.data());


	CalcExtrapolationsBC(unGPU, vnGPU, PnGPU, TnGPU);
	
	//estadoPredicted = estado;

	//U.resize(NPOINTSXD*NPOINTSYD);

	concurrency::parallel_for_each(U1nGPU.extent,
		[=](concurrency::index<2> idx)restrict(amp)
		{
			U1nGPU[idx] = CalcDensity(PnGPU[idx], RD, TnGPU[idx]);
			U2nGPU[idx] = U1nGPU[idx] * unGPU[idx];
			U3nGPU[idx] = U1nGPU[idx] * vnGPU[idx];

			const float VQuad = (unGPU[idx] * unGPU[idx]) + (vnGPU[idx] * vnGPU[idx]);

			U4nGPU[idx] = CalcEt(U1nGPU[idx], CvD, TnGPU[idx], VQuad);
		}
	);
	
	//for (unsigned int i = 0; i < NPOINTSY; i++)
	//{
	//	for (unsigned int j = 0; j < NPOINTSX; j++)
	//	{
	//		U[j + (i*NPOINTSX)] = CalcU(estado[j + (i*NPOINTSX)]);
	//	}
	//}


	float totalTime = 0.0f;


	concurrency::extent<2> topComputeDomain = concurrency::extent<2>((NPOINTSYD / 2) - 2, (NPOINTSXD - 2));

	concurrency::extent<2> bottomComputeDomain = concurrency::extent<2>((NPOINTSYD / 2) - 2, (NPOINTSXD - 2));

	concurrency::extent<2> leftComputeDomain = concurrency::extent<2>(2, (NPOINTSXD / 2) - (NPOINTSPXD / 2) - 1);

	concurrency::extent<2> rightComputeDomain = concurrency::extent<2>(2, (NPOINTSXD / 2) - (NPOINTSPXD / 2) - 1);

	try 
	{

		for (unsigned int it = itInicial; it < NITERATIONSD; it++)
		{
			// Calculate time step;

			float dt = 1e-8f;

			//float maxvijl = 0.0f;
			//		
			//for (unsigned int i = 1; i < NPOINTSYD - 1; i++)
			//{
			//	for (unsigned int j = 1; j < NPOINTSXD - 1; j++)
			//	{
			//		const float vijl = (4.0 / 3.0)*(gamaD)*(mun[j + (i * NPOINTSX)] / prandtlNumberD) / U[j + (i*NPOINTSX)].u;
			//
			//		if (vijl > maxvijl)
			//			maxvijl = vijl;
			//	}
			//}
			//
			//for (unsigned int i = 1; i < NPOINTSY - 1; i++)
			//{
			//	for (unsigned int j = 1; j < NPOINTSX - 1; j++)
			//	{
			//		const State& estadoij = estado[j + (i*NPOINTSX)];
			//
			//		const float aij = sqrt(gama * R * estadoij.T);
			//
			//		constexpr float suminverse = (1.0 / (dx*dx)) + (1.0 / (dy*dy));
			//
			//		//const float deltatcfl = (0.9/(1.0+(2.0/minMeshReynolds))) / ((fabs(estadoij.u) / dx) + (fabs(estadoij.v) / dy) + (aij*sqrt(suminverse)));
			//		const float deltatcfl = 0.5 /( (fabs(estadoij.u) / dx) + (fabs(estadoij.v) / dy) + (aij*sqrt(suminverse)) + (2.0 * maxvijl * suminverse) );
			//		if (deltatcfl < dt)
			//		{
			//			dt = deltatcfl;
			//		}
			//	}
			//}
			//


			// CALCULATE ALL PREDICTED DERIVATIVES

			// Plus in book -> minus in the array.

			// TOP PART

			concurrency::parallel_for_each(topComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + 1;
					const int j = idx[1] + 1;

					const float E1ipj = GetE1Predictor(i, j + 1, U1nGPU, unGPU);
					const float E2ipj = GetE2Predictor(i, j + 1, U1nGPU, unGPU, vnGPU, PnGPU, munGPU);
					const float E3ipj = GetE3Predictor(i, j + 1, U1nGPU, unGPU, vnGPU, munGPU);
					const float E4ipj = GetE4Predictor(i, j + 1, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float E1ij = GetE1Predictor(i, j, U1nGPU, unGPU);
					const float E2ij = GetE2Predictor(i, j, U1nGPU, unGPU, vnGPU, PnGPU, munGPU);
					const float E3ij = GetE3Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float E4ij = GetE4Predictor(i, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float F1ijp = GetF1Predictor(i - 1, j, U1nGPU, vnGPU);
					const float F2ijp = GetF2Predictor(i - 1, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float F3ijp = GetF3Predictor(i - 1, j, U1nGPU, unGPU, vnGPU, munGPU, PnGPU);
					const float F4ijp = GetF4Predictor(i - 1, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float F1ij = GetF1Predictor(i, j, U1nGPU, vnGPU);
					const float F2ij = GetF2Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float F3ij = GetF3Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU, PnGPU);
					const float F4ij = GetF4Predictor(i, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					dU1pdtGPU[i][j] = -((E1ipj - E1ij) / dxD) - ((F1ijp - F1ij) / dyD);
					dU2pdtGPU[i][j] = -((E2ipj - E2ij) / dxD) - ((F2ijp - F2ij) / dyD);
					dU3pdtGPU[i][j] = -((E3ipj - E3ij) / dxD) - ((F3ijp - F3ij) / dyD);
					dU4pdtGPU[i][j] = -((E4ipj - E4ij) / dxD) - ((F4ijp - F4ij) / dyD);

				});

			//for (int i = 1; i < (NPOINTSY/2) - 1; i++)
			//{
			//	for (int j = 1; j < NPOINTSX - 1; j++)
			//	{
			//		const float E1ipj = GetE1Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E2ipj = GetE2Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E3ipj = GetE3Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E4ipj = GetE4Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//
			//		const float E1ij = GetE1Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E2ij = GetE2Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E3ij = GetE3Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E4ij = GetE4Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ijp = GetF1Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F2ijp = GetF2Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F3ijp = GetF3Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F4ijp = GetF4Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ij = GetF1Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F2ij = GetF2Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F3ij = GetF3Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F4ij = GetF4Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		const State Eipj = { E1ipj, E2ipj, E3ipj, E4ipj };
			//		const State Eij = { E1ij, E2ij, E3ij, E4ij };
			//
			//		const State Fijp = { F1ijp, F2ijp, F3ijp, F4ijp };
			//		const State Fij = { F1ij, F2ij, F3ij, F4ij };
			//
			//		dUdtpredicted[j + (i * NPOINTSX)] = -((Eipj - Eij) / dx) - ((Fijp - Fij) / dy);
			//	}
			//}

			// BOTTOM PART

			concurrency::parallel_for_each(bottomComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) + 1;
					const int j = idx[1] + 1;

					const float E1ipj = GetE1Predictor(i, j + 1, U1nGPU, unGPU);
					const float E2ipj = GetE2Predictor(i, j + 1, U1nGPU, unGPU, vnGPU, PnGPU, munGPU);
					const float E3ipj = GetE3Predictor(i, j + 1, U1nGPU, unGPU, vnGPU, munGPU);
					const float E4ipj = GetE4Predictor(i, j + 1, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float E1ij = GetE1Predictor(i, j, U1nGPU, unGPU);
					const float E2ij = GetE2Predictor(i, j, U1nGPU, unGPU, vnGPU, PnGPU, munGPU);
					const float E3ij = GetE3Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float E4ij = GetE4Predictor(i, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float F1ijp = GetF1Predictor(i - 1, j, U1nGPU, vnGPU);
					const float F2ijp = GetF2Predictor(i - 1, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float F3ijp = GetF3Predictor(i - 1, j, U1nGPU, unGPU, vnGPU, munGPU, PnGPU);
					const float F4ijp = GetF4Predictor(i - 1, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float F1ij = GetF1Predictor(i, j, U1nGPU, vnGPU);
					const float F2ij = GetF2Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float F3ij = GetF3Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU, PnGPU);
					const float F4ij = GetF4Predictor(i, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					dU1pdtGPU[i][j] = -((E1ipj - E1ij) / dxD) - ((F1ijp - F1ij) / dyD);
					dU2pdtGPU[i][j] = -((E2ipj - E2ij) / dxD) - ((F2ijp - F2ij) / dyD);
					dU3pdtGPU[i][j] = -((E3ipj - E3ij) / dxD) - ((F3ijp - F3ij) / dyD);
					dU4pdtGPU[i][j] = -((E4ipj - E4ij) / dxD) - ((F4ijp - F4ij) / dyD);
				}
			);

			//for (int i = (NPOINTSY/2)+1; i < NPOINTSY - 1; i++)
			//{
			//	for (int j = 1; j < NPOINTSX - 1; j++)
			//	{
			//		const float E1ipj = GetE1Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E2ipj = GetE2Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E3ipj = GetE3Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E4ipj = GetE4Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//
			//		const float E1ij = GetE1Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E2ij = GetE2Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E3ij = GetE3Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E4ij = GetE4Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ijp = GetF1Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F2ijp = GetF2Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F3ijp = GetF3Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F4ijp = GetF4Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ij = GetF1Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F2ij = GetF2Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F3ij = GetF3Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F4ij = GetF4Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		const State Eipj = { E1ipj, E2ipj, E3ipj, E4ipj };
			//		const State Eij = { E1ij, E2ij, E3ij, E4ij };
			//
			//		const State Fijp = { F1ijp, F2ijp, F3ijp, F4ijp };
			//		const State Fij = { F1ij, F2ij, F3ij, F4ij };
			//
			//		dUdtpredicted[j + (i * NPOINTSX)] = -((Eipj - Eij) / dx) - ((Fijp - Fij) / dy);
			//	}
			//}

			// LEFT PART

			concurrency::parallel_for_each(leftComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) - 1;
					const int j = idx[1] + 1;

					const float E1ipj = GetE1Predictor(i, j + 1, U1nGPU, unGPU);
					const float E2ipj = GetE2Predictor(i, j + 1, U1nGPU, unGPU, vnGPU, PnGPU, munGPU);
					const float E3ipj = GetE3Predictor(i, j + 1, U1nGPU, unGPU, vnGPU, munGPU);
					const float E4ipj = GetE4Predictor(i, j + 1, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float E1ij = GetE1Predictor(i, j, U1nGPU, unGPU);
					const float E2ij = GetE2Predictor(i, j, U1nGPU, unGPU, vnGPU, PnGPU, munGPU);
					const float E3ij = GetE3Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float E4ij = GetE4Predictor(i, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float F1ijp = GetF1Predictor(i - 1, j, U1nGPU, vnGPU);
					const float F2ijp = GetF2Predictor(i - 1, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float F3ijp = GetF3Predictor(i - 1, j, U1nGPU, unGPU, vnGPU, munGPU, PnGPU);
					const float F4ijp = GetF4Predictor(i - 1, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float F1ij = GetF1Predictor(i, j, U1nGPU, vnGPU);
					const float F2ij = GetF2Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float F3ij = GetF3Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU, PnGPU);
					const float F4ij = GetF4Predictor(i, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					dU1pdtGPU[i][j] = -((E1ipj - E1ij) / dxD) - ((F1ijp - F1ij) / dyD);
					dU2pdtGPU[i][j] = -((E2ipj - E2ij) / dxD) - ((F2ijp - F2ij) / dyD);
					dU3pdtGPU[i][j] = -((E3ipj - E3ij) / dxD) - ((F3ijp - F3ij) / dyD);
					dU4pdtGPU[i][j] = -((E4ipj - E4ij) / dxD) - ((F4ijp - F4ij) / dyD);
				}
			);

			//for (unsigned int i = (NPOINTSY / 2) -1; i < (NPOINTSY / 2) +1; i++)
			//{
			//	for (unsigned int j = 1; j < (NPOINTSX / 2) - (NPOINTSPX / 2); j++)
			//	{
			//		const float E1ipj = GetE1Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E2ipj = GetE2Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E3ipj = GetE3Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E4ipj = GetE4Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//
			//		const float E1ij = GetE1Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E2ij = GetE2Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E3ij = GetE3Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E4ij = GetE4Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ijp = GetF1Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F2ijp = GetF2Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F3ijp = GetF3Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F4ijp = GetF4Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ij = GetF1Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F2ij = GetF2Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F3ij = GetF3Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F4ij = GetF4Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		const State Eipj = { E1ipj, E2ipj, E3ipj, E4ipj };
			//		const State Eij = { E1ij, E2ij, E3ij, E4ij };
			//
			//		const State Fijp = { F1ijp, F2ijp, F3ijp, F4ijp };
			//		const State Fij = { F1ij, F2ij, F3ij, F4ij };
			//
			//		dUdtpredicted[j + (i * NPOINTSX)] = -((Eipj - Eij) / dx) - ((Fijp - Fij) / dy);
			//	}
			//}

			// RIGHT PART
			concurrency::parallel_for_each(rightComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) - 1;
					const int j = idx[1] + (NPOINTSXD / 2) + (NPOINTSPXD / 2);

					const float E1ipj = GetE1Predictor(i, j + 1, U1nGPU, unGPU);
					const float E2ipj = GetE2Predictor(i, j + 1, U1nGPU, unGPU, vnGPU, PnGPU, munGPU);
					const float E3ipj = GetE3Predictor(i, j + 1, U1nGPU, unGPU, vnGPU, munGPU);
					const float E4ipj = GetE4Predictor(i, j + 1, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float E1ij = GetE1Predictor(i, j, U1nGPU, unGPU);
					const float E2ij = GetE2Predictor(i, j, U1nGPU, unGPU, vnGPU, PnGPU, munGPU);
					const float E3ij = GetE3Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float E4ij = GetE4Predictor(i, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float F1ijp = GetF1Predictor(i - 1, j, U1nGPU, vnGPU);
					const float F2ijp = GetF2Predictor(i - 1, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float F3ijp = GetF3Predictor(i - 1, j, U1nGPU, unGPU, vnGPU, munGPU, PnGPU);
					const float F4ijp = GetF4Predictor(i - 1, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);

					const float F1ij = GetF1Predictor(i, j, U1nGPU, vnGPU);
					const float F2ij = GetF2Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU);
					const float F3ij = GetF3Predictor(i, j, U1nGPU, unGPU, vnGPU, munGPU, PnGPU);
					const float F4ij = GetF4Predictor(i, j, U4nGPU, unGPU, vnGPU, munGPU, TnGPU, knGPU, PnGPU);


					dU1pdtGPU[i][j] = -((E1ipj - E1ij) / dxD) - ((F1ijp - F1ij) / dyD);
					dU2pdtGPU[i][j] = -((E2ipj - E2ij) / dxD) - ((F2ijp - F2ij) / dyD);
					dU3pdtGPU[i][j] = -((E3ipj - E3ij) / dxD) - ((F3ijp - F3ij) / dyD);
					dU4pdtGPU[i][j] = -((E4ipj - E4ij) / dxD) - ((F4ijp - F4ij) / dyD);
				}
			);

			dU1pdtGPU.synchronize();
			dU2pdtGPU.synchronize();
			dU3pdtGPU.synchronize();
			dU4pdtGPU.synchronize();
			
			//for (unsigned int i = (NPOINTSY / 2) - 1; i < (NPOINTSY / 2) + 1; i++)
			//{
			//	for (unsigned int j = (NPOINTSX / 2) + (NPOINTSPX / 2); j < NPOINTSX - 1; j++)
			//	{
			//		const float E1ipj = GetE1Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E2ipj = GetE2Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E3ipj = GetE3Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//		const float E4ipj = GetE4Predictor(estado, i, j + 1, NPOINTSX, NPOINTSY);
			//
			//		const float E1ij = GetE1Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E2ij = GetE2Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E3ij = GetE3Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float E4ij = GetE4Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ijp = GetF1Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F2ijp = GetF2Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F3ijp = GetF3Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//		const float F4ijp = GetF4Predictor(estado, i - 1, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ij = GetF1Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F2ij = GetF2Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F3ij = GetF3Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//		const float F4ij = GetF4Predictor(estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		const State Eipj = { E1ipj, E2ipj, E3ipj, E4ipj };
			//		const State Eij = { E1ij, E2ij, E3ij, E4ij };
			//
			//		const State Fijp = { F1ijp, F2ijp, F3ijp, F4ijp };
			//		const State Fij = { F1ij, F2ij, F3ij, F4ij };
			//
			//		dUdtpredicted[j + (i * NPOINTSX)] = -((Eipj - Eij) / dx) - ((Fijp - Fij) / dy);
			//	}
			//}

			// CALCULATE ALL THE PREDICTED VALUES

			//TOP PART
			//FALTA IMPLEMENTAR VISCOSIDADE NUMÉRICA;

			concurrency::parallel_for_each(topComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + 1;
					const int j = idx[1] + 1;

					U1pnGPU[i][j] = U1nGPU[i][j] + (dU1pdtGPU[i][j] * dt);
					U2pnGPU[i][j] = U2nGPU[i][j] + (dU2pdtGPU[i][j] * dt);
					U3pnGPU[i][j] = U3nGPU[i][j] + (dU3pdtGPU[i][j] * dt);
					U4pnGPU[i][j] = U4nGPU[i][j] + (dU4pdtGPU[i][j] * dt);

					const float PXCoef = concurrency::fast_math::fabsf(PnGPU[i][j + 1] - (2.0f*PnGPU[i][j]) + PnGPU[i][j - 1]) / 
						(PnGPU[i][j + 1] + (2.0f*PnGPU[i][j]) + PnGPU[i][j - 1]);
					const float PYCoef = concurrency::fast_math::fabsf(PnGPU[i - 1][j] - (2.0f * PnGPU[i][j]) + PnGPU[i+1][j]) / 
						(PnGPU[i - 1][j] + (2.0f * PnGPU[i][j]) + PnGPU[i + 1][j]);

					const float xcoef = CxD * PXCoef;
					const float ycoef = CyD * PYCoef;
					
					const float Srho = CalcSrho(i, j, U1pnGPU, xcoef, ycoef);
					const float Srhou = CalcSrhou(i, j, U2pnGPU, xcoef, ycoef);
					const float Srhov = CalcSrhov(i, j, U3pnGPU, xcoef, ycoef);
					const float SEt = CalcSEt(i, j, U4pnGPU, xcoef, ycoef);

					U1pnGPU[i][j] += Srho;
					U2pnGPU[i][j] += Srhou;
					U3pnGPU[i][j] += Srhov;
					U4pnGPU[i][j] += SEt;

					upnGPU[i][j] = U2pnGPU[i][j] / U1pnGPU[i][j];
					vpnGPU[i][j] = U3pnGPU[i][j] / U1pnGPU[i][j];
					TpnGPU[i][j] = ((U4pnGPU[i][j] / U1pnGPU[i][j]) - (((upnGPU[i][j] * upnGPU[i][j]) + (vpnGPU[i][j] * vpnGPU[i][j])) / 2.0f)) / CvD;
					PpnGPU[i][j] = U1pnGPU[i][j] * RD * TpnGPU[i][j];

					munGPU[i][j] = CalcViscSutherland(muRefD, tempRefD, TpnGPU[i][j]);
					knGPU[i][j] = Calckk(munGPU[i][j], CpD, prandtlNumberD);
				}
			);

			//for (int i = 1; i < (NPOINTSY / 2) - 1; i++)
			//{
			//	for (int j = 1; j < NPOINTSX - 1; j++)
			//	{
			//		Upredicted[j + (i*NPOINTSX)] = U[j + (i*NPOINTSX)] +
			//			(dUdtpredicted[j + (i * NPOINTSX)] * dt);
			//
			//		const float Srho = CalcSrho(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhou = CalcSrhou(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhov = CalcSrhov(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float SEt = CalcSEt(U, estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		Upredicted[j + (i*NPOINTSX)].u += Srho;
			//		Upredicted[j + (i*NPOINTSX)].v += Srhou;
			//		Upredicted[j + (i*NPOINTSX)].P += Srhov;
			//		Upredicted[j + (i*NPOINTSX)].T += SEt;
			//
			//		estadoPredicted[j + (i*NPOINTSX)] = CalcState(Upredicted[j + (i*NPOINTSX)]);
			//
			//		mu[j + (i*NPOINTSX)] = CalcViscSutherland(muRef, tempRef, estadoPredicted[j + (i*NPOINTSX)].T);
			//		k[j + (i*NPOINTSX)] = Calck(mu[j + (i*NPOINTSX)]);
			//		
			//		if (estadoPredicted[j + (i*NPOINTSX)].T < 0.0)
			//		{
			//			std::cout << "Erro i: " << i << " j: " << j << " on predicted values" << std::endl;
			//		}
			//	}
			//}

			// BOTTOM PART
			//FALTA IMPLEMENTAR VISCOSIDADE NUMÉRICA;
			concurrency::parallel_for_each(bottomComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) + 1;
					const int j = idx[1] + 1;

					U1pnGPU[i][j] = U1nGPU[i][j] + (dU1pdtGPU[i][j] * dt);
					U2pnGPU[i][j] = U2nGPU[i][j] + (dU2pdtGPU[i][j] * dt);
					U3pnGPU[i][j] = U3nGPU[i][j] + (dU3pdtGPU[i][j] * dt);
					U4pnGPU[i][j] = U4nGPU[i][j] + (dU4pdtGPU[i][j] * dt);

					const float PXCoef = concurrency::fast_math::fabsf(PnGPU[i][j + 1] - (2.0f * PnGPU[i][j]) + PnGPU[i][j - 1]) / 
						(PnGPU[i][j + 1] + (2.0f * PnGPU[i][j]) + PnGPU[i][j - 1]);
					const float PYCoef = concurrency::fast_math::fabsf(PnGPU[i - 1][j] - (2.0f * PnGPU[i][j]) + PnGPU[i + 1][j]) / 
						(PnGPU[i - 1][j] + (2.0f * PnGPU[i][j]) + PnGPU[i + 1][j]);

					const float xcoef = CxD * PXCoef;
					const float ycoef = CyD * PYCoef;

					const float Srho = CalcSrho(i, j, U1pnGPU, xcoef, ycoef);
					const float Srhou = CalcSrhou(i, j, U2pnGPU, xcoef, ycoef);
					const float Srhov = CalcSrhov(i, j, U3pnGPU, xcoef, ycoef);
					const float SEt = CalcSEt(i, j, U4pnGPU, xcoef, ycoef);

					U1pnGPU[i][j] += Srho;
					U2pnGPU[i][j] += Srhou;
					U3pnGPU[i][j] += Srhov;
					U4pnGPU[i][j] += SEt;

					upnGPU[i][j] = U2pnGPU[i][j] / U1pnGPU[i][j];
					vpnGPU[i][j] = U3pnGPU[i][j] / U1pnGPU[i][j];
					TpnGPU[i][j] = ((U4pnGPU[i][j] / U1pnGPU[i][j]) - ((upnGPU[i][j] * upnGPU[i][j]) + (vpnGPU[i][j] * vpnGPU[i][j]) / 2.0f)) / CvD;
					PpnGPU[i][j] = U1pnGPU[i][j] * RD * TpnGPU[i][j];

					munGPU[i][j] = CalcViscSutherland(muRefD, tempRefD, TpnGPU[i][j]);
					knGPU[i][j] = Calckk(munGPU[i][j], CpD, prandtlNumberD);
				}
			);

			//for (int i = (NPOINTSY / 2) + 1; i < NPOINTSY - 1; i++)
			//{
			//	for (int j = 1; j < NPOINTSX - 1; j++)
			//	{
			//		Upredicted[j + (i*NPOINTSX)] = U[j + (i*NPOINTSX)] +
			//			(dUdtpredicted[j + (i * NPOINTSX)] * dt);
			//
			//		const float Srho = CalcSrho(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhou = CalcSrhou(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhov = CalcSrhov(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float SEt = CalcSEt(U, estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		Upredicted[j + (i*NPOINTSX)].u += Srho;
			//		Upredicted[j + (i*NPOINTSX)].v += Srhou;
			//		Upredicted[j + (i*NPOINTSX)].P += Srhov;
			//		Upredicted[j + (i*NPOINTSX)].T += SEt;
			//
			//		estadoPredicted[j + (i*NPOINTSX)] = CalcState(Upredicted[j + (i*NPOINTSX)]);
			//
			//		mu[j + (i*NPOINTSX)] = CalcViscSutherland(muRef, tempRef, estadoPredicted[j + (i*NPOINTSX)].T);
			//		k[j + (i*NPOINTSX)] = Calck(mu[j + (i*NPOINTSX)]);
			//		
			//		if (estadoPredicted[j + (i*NPOINTSX)].T < 0.0)
			//		{
			//			std::cout << "Erro i: " << i << " j: " << j << " on predicted values" << std::endl;
			//		}
			//	}
			//}

			// LEFT PART

			//FALTA IMPLEMENTAR VISCOSIDADE NUMÉRICA;
			concurrency::parallel_for_each(bottomComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) - 1;
					const int j = idx[1] + 1;

					U1pnGPU[i][j] = U1nGPU[i][j] + (dU1pdtGPU[i][j] * dt);
					U2pnGPU[i][j] = U2nGPU[i][j] + (dU2pdtGPU[i][j] * dt);
					U3pnGPU[i][j] = U3nGPU[i][j] + (dU3pdtGPU[i][j] * dt);
					U4pnGPU[i][j] = U4nGPU[i][j] + (dU4pdtGPU[i][j] * dt);

					const float PXCoef = concurrency::fast_math::fabsf(PnGPU[i][j + 1] - (2.0f * PnGPU[i][j]) + PnGPU[i][j - 1]) / 
						(PnGPU[i][j + 1] + (2.0f * PnGPU[i][j]) + PnGPU[i][j - 1]);
					const float PYCoef = concurrency::fast_math::fabsf(PnGPU[i - 1][j] - (2.0f * PnGPU[i][j]) + PnGPU[i + 1][j]) / 
						(PnGPU[i - 1][j] + (2.0f * PnGPU[i][j]) + PnGPU[i + 1][j]);

					const float xcoef = CxD * PXCoef;
					const float ycoef = CyD * PYCoef;

					const float Srho = CalcSrho(i, j, U1pnGPU, xcoef, ycoef);
					const float Srhou = CalcSrhou(i, j, U2pnGPU, xcoef, ycoef);
					const float Srhov = CalcSrhov(i, j, U3pnGPU, xcoef, ycoef);
					const float SEt = CalcSEt(i, j, U4pnGPU, xcoef, ycoef);

					U1pnGPU[i][j] += Srho;
					U2pnGPU[i][j] += Srhou;
					U3pnGPU[i][j] += Srhov;
					U4pnGPU[i][j] += SEt;

					upnGPU[i][j] = U2pnGPU[i][j] / U1pnGPU[i][j];
					vpnGPU[i][j] = U3pnGPU[i][j] / U1pnGPU[i][j];
					TpnGPU[i][j] = ((U4pnGPU[i][j] / U1pnGPU[i][j]) - ((upnGPU[i][j] * upnGPU[i][j]) + (vpnGPU[i][j] * vpnGPU[i][j]) / 2.0f)) / CvD;
					PpnGPU[i][j] = U1pnGPU[i][j] * RD * TpnGPU[i][j];

					munGPU[i][j] = CalcViscSutherland(muRefD, tempRefD, TpnGPU[i][j]);
					knGPU[i][j] = Calckk(munGPU[i][j], CpD, prandtlNumberD);
				}
			);
			//for (unsigned int i = (NPOINTSY / 2) - 1; i < (NPOINTSY / 2) + 1; i++)
			//{
			//	for (unsigned int j = 1; j < (NPOINTSX / 2) - (NPOINTSPX / 2); j++)
			//	{
			//		Upredicted[j + (i*NPOINTSX)] = U[j + (i*NPOINTSX)] +
			//			(dUdtpredicted[j + (i * NPOINTSX)] * dt);
			//
			//		const float Srho = CalcSrho(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhou = CalcSrhou(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhov = CalcSrhov(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float SEt = CalcSEt(U, estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		Upredicted[j + (i*NPOINTSX)].u += Srho;
			//		Upredicted[j + (i*NPOINTSX)].v += Srhou;
			//		Upredicted[j + (i*NPOINTSX)].P += Srhov;
			//		Upredicted[j + (i*NPOINTSX)].T += SEt;
			//
			//		estadoPredicted[j + (i*NPOINTSX)] = CalcState(Upredicted[j + (i*NPOINTSX)]);
			//
			//		mu[j + (i*NPOINTSX)] = CalcViscSutherland(muRef, tempRef, estadoPredicted[j + (i*NPOINTSX)].T);
			//		k[j + (i*NPOINTSX)] = Calck(mu[j + (i*NPOINTSX)]);
			//
			//
			//		if (estadoPredicted[j + (i*NPOINTSX)].T < 0.0)
			//		{
			//			std::cout << "Erro i: " << i << " j: " << j << " on predicted values" << std::endl;
			//		}
			//	}
			//}

			// RIGHT PART
			//FALTA IMPLEMENTAR VISCOSIDADE NUMÉRICA;
			concurrency::parallel_for_each(bottomComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) - 1;
					const int j = idx[1] + (NPOINTSXD / 2) + (NPOINTSPXD / 2);

					U1pnGPU[i][j] = U1nGPU[i][j] + (dU1pdtGPU[i][j] * dt);
					U2pnGPU[i][j] = U2nGPU[i][j] + (dU2pdtGPU[i][j] * dt);
					U3pnGPU[i][j] = U3nGPU[i][j] + (dU3pdtGPU[i][j] * dt);
					U4pnGPU[i][j] = U4nGPU[i][j] + (dU4pdtGPU[i][j] * dt);

					const float PXCoef = concurrency::fast_math::fabsf(PnGPU[i][j + 1] - (2.0f * PnGPU[i][j]) + PnGPU[i][j - 1]) / 
						(PnGPU[i][j + 1] + (2.0f * PnGPU[i][j]) + PnGPU[i][j - 1]);
					const float PYCoef = concurrency::fast_math::fabsf(PnGPU[i - 1][j] - (2.0f * PnGPU[i][j]) + PnGPU[i + 1][j]) / 
						(PnGPU[i - 1][j] + (2.0f * PnGPU[i][j]) + PnGPU[i + 1][j]);

					const float xcoef = CxD * PXCoef;
					const float ycoef = CyD * PYCoef;

					const float Srho = CalcSrho(i, j, U1pnGPU, xcoef, ycoef);
					const float Srhou = CalcSrhou(i, j, U2pnGPU, xcoef, ycoef);
					const float Srhov = CalcSrhov(i, j, U3pnGPU, xcoef, ycoef);
					const float SEt = CalcSEt(i, j, U4pnGPU, xcoef, ycoef);

					U1pnGPU[i][j] += Srho;
					U2pnGPU[i][j] += Srhou;
					U3pnGPU[i][j] += Srhov;
					U4pnGPU[i][j] += SEt;

					upnGPU[i][j] = U2pnGPU[i][j] / U1pnGPU[i][j];
					vpnGPU[i][j] = U3pnGPU[i][j] / U1pnGPU[i][j];
					TpnGPU[i][j] = ((U4pnGPU[i][j] / U1pnGPU[i][j]) - (((upnGPU[i][j] * upnGPU[i][j]) + (vpnGPU[i][j] * vpnGPU[i][j])) / 2.0f)) / CvD;
					PpnGPU[i][j] = U1pnGPU[i][j] * RD * TpnGPU[i][j];

					munGPU[i][j] = CalcViscSutherland(muRefD, tempRefD, TpnGPU[i][j]);
					knGPU[i][j] = Calckk(munGPU[i][j], CpD, prandtlNumberD);
				}
			);

			//for (int i = 1; i < NPOINTSYD-1; i++)
			//	for (int j = 1; j < NPOINTSXD-1; j++)
			//	{
			//		if (TpnGPU[i][j] <= 0.0f)
			//			std::cout << "pi: " << i << " j: " << j << " val1: " << TpnGPU[i][j] << std::endl;
			//
			//		if (TpnGPU[i][j] <= 0.0f)
			//			std::cout << "pi: " << i << " j: " << j << " val2: " << TpnGPU[i][j] << std::endl;
			//
			//		if (TpnGPU[i][j] <= 0.0f)
			//			std::cout << "pi: " << i << " j: " << j << " val3: " << TpnGPU[i][j] << std::endl;
			//
			//		if (TpnGPU[i][j] <= 0.0f)
			//			std::cout << "pi: " << i << " j: " << j << " val4: " << TpnGPU[i][j] << std::endl;
			//	}


			//for (unsigned int i = (NPOINTSY / 2) - 1; i < (NPOINTSY / 2) + 1; i++)
			//{
			//	for (unsigned int j = (NPOINTSX / 2) + (NPOINTSPX / 2); j < NPOINTSX - 1; j++)
			//	{
			//		Upredicted[j + (i*NPOINTSX)] = U[j + (i*NPOINTSX)] +
			//			(dUdtpredicted[j + (i * NPOINTSX)] * dt);
			//
			//		const float Srho = CalcSrho(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhou = CalcSrhou(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhov = CalcSrhov(U, estado, i, j, NPOINTSX, NPOINTSY);
			//		const float SEt = CalcSEt(U, estado, i, j, NPOINTSX, NPOINTSY);
			//
			//		Upredicted[j + (i*NPOINTSX)].u += Srho;
			//		Upredicted[j + (i*NPOINTSX)].v += Srhou;
			//		Upredicted[j + (i*NPOINTSX)].P += Srhov;
			//		Upredicted[j + (i*NPOINTSX)].T += SEt;
			//
			//		estadoPredicted[j + (i*NPOINTSX)] = CalcState(Upredicted[j + (i*NPOINTSX)]);
			//
			//		mu[j + (i*NPOINTSX)] = CalcViscSutherland(muRef, tempRef, estadoPredicted[j + (i*NPOINTSX)].T);
			//		k[j + (i*NPOINTSX)] = Calck(mu[j + (i*NPOINTSX)]);
			//		
			//		if (estadoPredicted[j + (i*NPOINTSX)].T < 0.0)
			//		{
			//			std::cout << "Erro i: " << i << " j: " << j << " on predicted values" << std::endl;
			//		}
			//	}
			//}

			// Calcula extrapolações para os valores preditos

			U1pnGPU.synchronize();
			U2pnGPU.synchronize();
			U3pnGPU.synchronize();
			U4pnGPU.synchronize();

			upnGPU.synchronize();
			vpnGPU.synchronize();
			PpnGPU.synchronize();
			TpnGPU.synchronize();

			CalcExtrapolationsBC(upnGPU, vpnGPU, PpnGPU, TpnGPU);

			upnGPU.synchronize();
			vpnGPU.synchronize();
			PpnGPU.synchronize();
			TpnGPU.synchronize();

			// CALCULATE ALL THE CORRECTED DERIVATIVES

			// TOP PART

			concurrency::parallel_for_each(topComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + 1;
					const int j = idx[1] + 1;

					const float E1imj = GetE1Corrector(i, j - 1, U1pnGPU, upnGPU);
					const float E2imj = GetE2Corrector(i, j - 1, U1pnGPU, upnGPU, vpnGPU, PpnGPU, munGPU);
					const float E3imj = GetE3Corrector(i, j - 1, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float E4imj = GetE4Corrector(i, j - 1, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float E1ij = GetE1Corrector(i, j, U1pnGPU, upnGPU);
					const float E2ij = GetE2Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, PpnGPU, munGPU);
					const float E3ij = GetE3Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float E4ij = GetE4Corrector(i, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float F1ijm = GetF1Corrector(i + 1, j, U1pnGPU, vpnGPU);
					const float F2ijm = GetF2Corrector(i + 1, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float F3ijm = GetF3Corrector(i + 1, j, U1pnGPU, upnGPU, vpnGPU, munGPU, PpnGPU);
					const float F4ijm = GetF4Corrector(i + 1, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float F1ij = GetF1Corrector(i, j, U1pnGPU, vpnGPU);
					const float F2ij = GetF2Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float F3ij = GetF3Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU, PpnGPU);
					const float F4ij = GetF4Corrector(i, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					dU1cdtGPU[i][j] = -((E1ij - E1imj) / dxD) - ((F1ij - F1ijm) / dyD);
					dU2cdtGPU[i][j] = -((E2ij - E2imj) / dxD) - ((F2ij - F2ijm) / dyD);
					dU3cdtGPU[i][j] = -((E3ij - E3imj) / dxD) - ((F3ij - F3ijm) / dyD);
					dU4cdtGPU[i][j] = -((E4ij - E4imj) / dxD) - ((F4ij - F4ijm) / dyD);
				}
			);

			//for (int i = 1; i < (NPOINTSY / 2) - 1; i++)
			//{
			//	for (int j = 1; j < NPOINTSX - 1; j++)
			//	{
			//		const float E1imj = GetE1Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E2imj = GetE2Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E3imj = GetE3Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E4imj = GetE4Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//
			//		const float E1ij = GetE1Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E2ij = GetE2Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E3ij = GetE3Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E4ij = GetE4Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ijm = GetF1Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F2ijm = GetF2Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F3ijm = GetF3Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F4ijm = GetF4Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ij = GetF1Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F2ij = GetF2Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F3ij = GetF3Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F4ij = GetF4Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		const State Eimj = { E1imj, E2imj, E3imj, E4imj };
			//		const State Eij = { E1ij, E2ij, E3ij, E4ij };
			//
			//		const State Fijm = { F1ijm, F2ijm, F3ijm, F4ijm };
			//		const State Fij = { F1ij, F2ij, F3ij, F4ij };
			//
			//		dUdtcorrected[j + (i * NPOINTSX)] = -((Eij - Eimj) / dx) - ((Fij - Fijm) / dy);
			//	}
			//}

			// BOTTOM PART
			concurrency::parallel_for_each(bottomComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) + 1;
					const int j = idx[1] + 1;

					const float E1imj = GetE1Corrector(i, j - 1, U1pnGPU, upnGPU);
					const float E2imj = GetE2Corrector(i, j - 1, U1pnGPU, upnGPU, vpnGPU, PpnGPU, munGPU);
					const float E3imj = GetE3Corrector(i, j - 1, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float E4imj = GetE4Corrector(i, j - 1, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float E1ij = GetE1Corrector(i, j, U1pnGPU, upnGPU);
					const float E2ij = GetE2Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, PpnGPU, munGPU);
					const float E3ij = GetE3Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float E4ij = GetE4Corrector(i, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float F1ijm = GetF1Corrector(i + 1, j, U1pnGPU, vpnGPU);
					const float F2ijm = GetF2Corrector(i + 1, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float F3ijm = GetF3Corrector(i + 1, j, U1pnGPU, upnGPU, vpnGPU, munGPU, PpnGPU);
					const float F4ijm = GetF4Corrector(i + 1, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float F1ij = GetF1Corrector(i, j, U1pnGPU, vpnGPU);
					const float F2ij = GetF2Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float F3ij = GetF3Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU, PpnGPU);
					const float F4ij = GetF4Corrector(i, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					dU1cdtGPU[i][j] = -((E1ij - E1imj) / dxD) - ((F1ij - F1ijm) / dyD);
					dU2cdtGPU[i][j] = -((E2ij - E2imj) / dxD) - ((F2ij - F2ijm) / dyD);
					dU3cdtGPU[i][j] = -((E3ij - E3imj) / dxD) - ((F3ij - F3ijm) / dyD);
					dU4cdtGPU[i][j] = -((E4ij - E4imj) / dxD) - ((F4ij - F4ijm) / dyD);
				}
			);

			//for (int i = (NPOINTSY / 2) + 1; i < NPOINTSY - 1; i++)
			//{
			//	for (int j = 1; j < NPOINTSX - 1; j++)
			//	{
			//		const float E1imj = GetE1Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E2imj = GetE2Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E3imj = GetE3Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E4imj = GetE4Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//
			//		const float E1ij = GetE1Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E2ij = GetE2Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E3ij = GetE3Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E4ij = GetE4Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ijm = GetF1Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F2ijm = GetF2Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F3ijm = GetF3Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F4ijm = GetF4Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ij = GetF1Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F2ij = GetF2Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F3ij = GetF3Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F4ij = GetF4Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		const State Eimj = { E1imj, E2imj, E3imj, E4imj };
			//		const State Eij = { E1ij, E2ij, E3ij, E4ij };
			//
			//		const State Fijm = { F1ijm, F2ijm, F3ijm, F4ijm };
			//		const State Fij = { F1ij, F2ij, F3ij, F4ij };
			//
			//		dUdtcorrected[j + (i * NPOINTSX)] = -((Eij - Eimj) / dx) - ((Fij - Fijm) / dy);
			//	}
			//}

			// LEFT PART
			concurrency::parallel_for_each(leftComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) - 1;
					const int j = idx[1] + 1;

					const float E1imj = GetE1Corrector(i, j - 1, U1pnGPU, upnGPU);
					const float E2imj = GetE2Corrector(i, j - 1, U1pnGPU, upnGPU, vpnGPU, PpnGPU, munGPU);
					const float E3imj = GetE3Corrector(i, j - 1, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float E4imj = GetE4Corrector(i, j - 1, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float E1ij = GetE1Corrector(i, j, U1pnGPU, upnGPU);
					const float E2ij = GetE2Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, PpnGPU, munGPU);
					const float E3ij = GetE3Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float E4ij = GetE4Corrector(i, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float F1ijm = GetF1Corrector(i + 1, j, U1pnGPU, vpnGPU);
					const float F2ijm = GetF2Corrector(i + 1, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float F3ijm = GetF3Corrector(i + 1, j, U1pnGPU, upnGPU, vpnGPU, munGPU, PpnGPU);
					const float F4ijm = GetF4Corrector(i + 1, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float F1ij = GetF1Corrector(i, j, U1pnGPU, vpnGPU);
					const float F2ij = GetF2Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float F3ij = GetF3Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU, PpnGPU);
					const float F4ij = GetF4Corrector(i, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					dU1cdtGPU[i][j] = -((E1ij - E1imj) / dxD) - ((F1ij - F1ijm) / dyD);
					dU2cdtGPU[i][j] = -((E2ij - E2imj) / dxD) - ((F2ij - F2ijm) / dyD);
					dU3cdtGPU[i][j] = -((E3ij - E3imj) / dxD) - ((F3ij - F3ijm) / dyD);
					dU4cdtGPU[i][j] = -((E4ij - E4imj) / dxD) - ((F4ij - F4ijm) / dyD);
				}
			);
			//for (unsigned int i = (NPOINTSY / 2)-1; i < (NPOINTSY / 2) + 1; i++)
			//{
			//	for (unsigned int j = 1; j < (NPOINTSX / 2) - (NPOINTSPX / 2); j++)
			//	{
			//		const float E1imj = GetE1Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E2imj = GetE2Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E3imj = GetE3Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E4imj = GetE4Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//
			//		const float E1ij = GetE1Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E2ij = GetE2Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E3ij = GetE3Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E4ij = GetE4Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ijm = GetF1Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F2ijm = GetF2Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F3ijm = GetF3Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F4ijm = GetF4Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ij = GetF1Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F2ij = GetF2Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F3ij = GetF3Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F4ij = GetF4Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		const State Eimj = { E1imj, E2imj, E3imj, E4imj };
			//		const State Eij = { E1ij, E2ij, E3ij, E4ij };
			//
			//		const State Fijm = { F1ijm, F2ijm, F3ijm, F4ijm };
			//		const State Fij = { F1ij, F2ij, F3ij, F4ij };
			//
			//		dUdtcorrected[j + (i * NPOINTSX)] = -((Eij - Eimj) / dx) - ((Fij - Fijm) / dy);
			//	}
			//}

			// RIGHT PART

			concurrency::parallel_for_each(rightComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) - 1;
					const int j = idx[1] + (NPOINTSXD / 2) + (NPOINTSPXD / 2);

					const float E1imj = GetE1Corrector(i, j - 1, U1pnGPU, upnGPU);
					const float E2imj = GetE2Corrector(i, j - 1, U1pnGPU, upnGPU, vpnGPU, PpnGPU, munGPU);
					const float E3imj = GetE3Corrector(i, j - 1, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float E4imj = GetE4Corrector(i, j - 1, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float E1ij = GetE1Corrector(i, j, U1pnGPU, upnGPU);
					const float E2ij = GetE2Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, PpnGPU, munGPU);
					const float E3ij = GetE3Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float E4ij = GetE4Corrector(i, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float F1ijm = GetF1Corrector(i + 1, j, U1pnGPU, vpnGPU);
					const float F2ijm = GetF2Corrector(i + 1, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float F3ijm = GetF3Corrector(i + 1, j, U1pnGPU, upnGPU, vpnGPU, munGPU, PpnGPU);
					const float F4ijm = GetF4Corrector(i + 1, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					const float F1ij = GetF1Corrector(i, j, U1pnGPU, vpnGPU);
					const float F2ij = GetF2Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU);
					const float F3ij = GetF3Corrector(i, j, U1pnGPU, upnGPU, vpnGPU, munGPU, PpnGPU);
					const float F4ij = GetF4Corrector(i, j, U4pnGPU, upnGPU, vpnGPU, munGPU, TpnGPU, knGPU, PpnGPU);

					dU1cdtGPU[i][j] = -((E1ij - E1imj) / dxD) - ((F1ij - F1ijm) / dyD);
					dU2cdtGPU[i][j] = -((E2ij - E2imj) / dxD) - ((F2ij - F2ijm) / dyD);
					dU3cdtGPU[i][j] = -((E3ij - E3imj) / dxD) - ((F3ij - F3ijm) / dyD);
					dU4cdtGPU[i][j] = -((E4ij - E4imj) / dxD) - ((F4ij - F4ijm) / dyD);
				}
			);

			dU1cdtGPU.synchronize();
			dU2cdtGPU.synchronize();
			dU3cdtGPU.synchronize();
			dU4cdtGPU.synchronize();

			//for(int i = 0; i < NPOINTSYD; i++)
			//	for (int j = 0; j < NPOINTSXD; j++)
			//	{
			//		if (dU1cdtGPU[i][j] != 0.0f)
			//			std::cout << "ci: " << i << " j: " << j << " val1: " << dU1cdtGPU[i][j] << std::endl;
			//
			//		if (dU2cdtGPU[i][j] != 0.0f)
			//			std::cout << "ci: " << i << " j: " << j << " val2: " << dU2cdtGPU[i][j] << std::endl;
			//
			//		if (dU3cdtGPU[i][j] != 0.0f)
			//			std::cout << "ci: " << i << " j: " << j << " val3: " << dU3cdtGPU[i][j] << std::endl;
			//
			//		if (dU4cdtGPU[i][j] != 0.0f)
			//			std::cout << "ci: " << i << " j: " << j << " val4: " << dU4cdtGPU[i][j] << std::endl;
			//	}
			//for (unsigned int i = (NPOINTSY / 2) - 1; i < (NPOINTSY / 2) + 1; i++)
			//{
			//	for (unsigned int j = (NPOINTSX / 2) + (NPOINTSPX / 2); j < NPOINTSX - 1; j++)
			//	{
			//		const float E1imj = GetE1Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E2imj = GetE2Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E3imj = GetE3Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//		const float E4imj = GetE4Corrector(estadoPredicted, i, j - 1, NPOINTSX, NPOINTSY);
			//
			//		const float E1ij = GetE1Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E2ij = GetE2Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E3ij = GetE3Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float E4ij = GetE4Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ijm = GetF1Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F2ijm = GetF2Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F3ijm = GetF3Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//		const float F4ijm = GetF4Corrector(estadoPredicted, i + 1, j, NPOINTSX, NPOINTSY);
			//
			//		const float F1ij = GetF1Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F2ij = GetF2Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F3ij = GetF3Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float F4ij = GetF4Corrector(estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		const State Eimj = { E1imj, E2imj, E3imj, E4imj };
			//		const State Eij = { E1ij, E2ij, E3ij, E4ij };
			//
			//		const State Fijm = { F1ijm, F2ijm, F3ijm, F4ijm };
			//		const State Fij = { F1ij, F2ij, F3ij, F4ij };
			//
			//		dUdtcorrected[j + (i * NPOINTSX)] = -((Eij - Eimj) / dx) - ((Fij - Fijm) / dy);
			//	}
			//}

			// CALCULATE ALL THE CORRECTED VALUES

			// TOP PART

			concurrency::parallel_for_each(topComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + 1;
					const int j = idx[1] + 1;

					U1nGPU[i][j] = U1nGPU[i][j] + ((dU1pdtGPU[i][j] + dU1cdtGPU[i][j]) * (dt / 2.0f));
					U2nGPU[i][j] = U2nGPU[i][j] + ((dU2pdtGPU[i][j] + dU2cdtGPU[i][j]) * (dt / 2.0f));
					U3nGPU[i][j] = U3nGPU[i][j] + ((dU3pdtGPU[i][j] + dU3cdtGPU[i][j]) * (dt / 2.0f));
					U4nGPU[i][j] = U4nGPU[i][j] + ((dU4pdtGPU[i][j] + dU4cdtGPU[i][j]) * (dt / 2.0f));

					const float PXCoef = concurrency::fast_math::fabsf(PpnGPU[i][j + 1] - (2.0f * PpnGPU[i][j]) + PpnGPU[i][j - 1]) /
						(PpnGPU[i][j + 1] + (2.0f * PpnGPU[i][j]) + PpnGPU[i][j - 1]);
					const float PYCoef = concurrency::fast_math::fabsf(PpnGPU[i - 1][j] - (2.0f * PpnGPU[i][j]) + PpnGPU[i + 1][j]) /
						(PpnGPU[i - 1][j] + (2.0f * PpnGPU[i][j]) + PpnGPU[i + 1][j]);

					const float xcoef = CxD * PXCoef;
					const float ycoef = CyD * PYCoef;

					const float Srho = CalcSrho(i, j, U1nGPU, xcoef, ycoef);
					const float Srhou = CalcSrhou(i, j, U2nGPU, xcoef, ycoef);
					const float Srhov = CalcSrhov(i, j, U3nGPU, xcoef, ycoef);
					const float SEt = CalcSEt(i, j, U4nGPU, xcoef, ycoef);

					U1pnGPU[i][j] += Srho;
					U2pnGPU[i][j] += Srhou;
					U3pnGPU[i][j] += Srhov;
					U4pnGPU[i][j] += SEt;
					
					unGPU[i][j] = U2nGPU[i][j] / U1nGPU[i][j];
					vnGPU[i][j] = U3nGPU[i][j] / U1nGPU[i][j];
					TnGPU[i][j] = ((U4nGPU[i][j] / U1nGPU[i][j]) - ((unGPU[i][j] * unGPU[i][j]) + (vnGPU[i][j] * vnGPU[i][j]) / 2.0f)) / CvD;
					PnGPU[i][j] = U1nGPU[i][j] * RD * TnGPU[i][j];

					munGPU[i][j] = CalcViscSutherland(muRefD, tempRefD, TnGPU[i][j]);
					knGPU[i][j] = Calckk(munGPU[i][j], CpD, prandtlNumberD);
				}
			);
			//for (int i = 1; i < (NPOINTSY / 2)-1; i++)
			//{
			//	for (int j = 1; j < NPOINTSX - 1; j++)
			//	{
			//		Ucorrected[j + (i*NPOINTSX)] = U[j + (i*NPOINTSX)] +
			//			((dUdtpredicted[j + (i * NPOINTSX)] + dUdtcorrected[j + (i * NPOINTSX)])*(dt / 2.0));
			//
			//		const float Srho = CalcSrho(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhou = CalcSrhou(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhov = CalcSrhov(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float SEt = CalcSEt(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		Ucorrected[j + (i*NPOINTSX)].u += Srho;
			//		Ucorrected[j + (i*NPOINTSX)].v += Srhou;
			//		Ucorrected[j + (i*NPOINTSX)].P += Srhov;
			//		Ucorrected[j + (i*NPOINTSX)].T += SEt;
			//
			//		estado[j + (i*NPOINTSX)] = CalcState(Ucorrected[j + (i*NPOINTSX)]);
			//
			//		mu[j + (i*NPOINTSX)] = CalcViscSutherland(muRef, tempRef, estado[j + (i*NPOINTSX)].T);
			//		k[j + (i*NPOINTSX)] = Calck(mu[j + (i*NPOINTSX)]);
			//
			//
			//		if (estado[j + (i*NPOINTSX)].T < 0.0)
			//		{
			//			std::cout << "Erro i: " << i << " j: " << j << " on corrected values" << std::endl;
			//		}
			//
			//	}
			//}

			// BOTTOM PART
			concurrency::parallel_for_each(bottomComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) + 1;
					const int j = idx[1] + 1;

					U1nGPU[i][j] = U1nGPU[i][j] + ((dU1pdtGPU[i][j] + dU1cdtGPU[i][j]) * (dt / 2.0f));
					U2nGPU[i][j] = U2nGPU[i][j] + ((dU2pdtGPU[i][j] + dU2cdtGPU[i][j]) * (dt / 2.0f));
					U3nGPU[i][j] = U3nGPU[i][j] + ((dU3pdtGPU[i][j] + dU3cdtGPU[i][j]) * (dt / 2.0f));
					U4nGPU[i][j] = U4nGPU[i][j] + ((dU4pdtGPU[i][j] + dU4cdtGPU[i][j]) * (dt / 2.0f));

					const float PXCoef = concurrency::fast_math::fabsf(PpnGPU[i][j + 1] - (2.0f * PpnGPU[i][j]) + PpnGPU[i][j - 1]) /
						(PpnGPU[i][j + 1] + (2.0f * PpnGPU[i][j]) + PpnGPU[i][j - 1]);
					const float PYCoef = concurrency::fast_math::fabsf(PpnGPU[i - 1][j] - (2.0f * PpnGPU[i][j]) + PpnGPU[i + 1][j]) /
						(PpnGPU[i - 1][j] + (2.0f * PpnGPU[i][j]) + PpnGPU[i + 1][j]);

					const float xcoef = CxD * PXCoef;
					const float ycoef = CyD * PYCoef;

					const float Srho = CalcSrho(i, j, U1nGPU, xcoef, ycoef);
					const float Srhou = CalcSrhou(i, j, U2nGPU, xcoef, ycoef);
					const float Srhov = CalcSrhov(i, j, U3nGPU, xcoef, ycoef);
					const float SEt = CalcSEt(i, j, U4nGPU, xcoef, ycoef);

					U1pnGPU[i][j] += Srho;
					U2pnGPU[i][j] += Srhou;
					U3pnGPU[i][j] += Srhov;
					U4pnGPU[i][j] += SEt;

					unGPU[i][j] = U2nGPU[i][j] / U1nGPU[i][j];
					vnGPU[i][j] = U3nGPU[i][j] / U1nGPU[i][j];
					TnGPU[i][j] = ((U4nGPU[i][j] / U1nGPU[i][j]) - ((unGPU[i][j] * unGPU[i][j]) + (vnGPU[i][j] * vnGPU[i][j]) / 2.0f)) / CvD;
					PnGPU[i][j] = U1nGPU[i][j] * RD * TnGPU[i][j];

					munGPU[i][j] = CalcViscSutherland(muRefD, tempRefD, TnGPU[i][j]);
					knGPU[i][j] = Calckk(munGPU[i][j], CpD, prandtlNumberD);
				}
			);

			//for (int i = (NPOINTSY / 2) + 1; i < NPOINTSY - 1; i++)
			//{
			//	for (int j = 1; j < NPOINTSX - 1; j++)
			//	{
			//		Ucorrected[j + (i*NPOINTSX)] = U[j + (i*NPOINTSX)] +
			//			((dUdtpredicted[j + (i * NPOINTSX)] + dUdtcorrected[j + (i * NPOINTSX)])*(dt / 2.0));
			//
			//		const float Srho = CalcSrho(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhou = CalcSrhou(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhov = CalcSrhov(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float SEt = CalcSEt(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		Ucorrected[j + (i*NPOINTSX)].u += Srho;
			//		Ucorrected[j + (i*NPOINTSX)].v += Srhou;
			//		Ucorrected[j + (i*NPOINTSX)].P += Srhov;
			//		Ucorrected[j + (i*NPOINTSX)].T += SEt;
			//
			//		estado[j + (i*NPOINTSX)] = CalcState(Ucorrected[j + (i*NPOINTSX)]);
			//
			//		mu[j + (i*NPOINTSX)] = CalcViscSutherland(muRef, tempRef, estado[j + (i*NPOINTSX)].T);
			//		k[j + (i*NPOINTSX)] = Calck(mu[j + (i*NPOINTSX)]);
			//
			//
			//		if (estado[j + (i*NPOINTSX)].T < 0.0)
			//		{
			//			std::cout << "Erro i: " << i << " j: " << j << " on corrected values" << std::endl;
			//		}
			//
			//	}
			//}

			// LEFT PART
			concurrency::parallel_for_each(leftComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) - 1;
					const int j = idx[1] + 1;

					U1nGPU[i][j] = U1nGPU[i][j] + ((dU1pdtGPU[i][j] + dU1cdtGPU[i][j]) * (dt / 2.0f));
					U2nGPU[i][j] = U2nGPU[i][j] + ((dU2pdtGPU[i][j] + dU2cdtGPU[i][j]) * (dt / 2.0f));
					U3nGPU[i][j] = U3nGPU[i][j] + ((dU3pdtGPU[i][j] + dU3cdtGPU[i][j]) * (dt / 2.0f));
					U4nGPU[i][j] = U4nGPU[i][j] + ((dU4pdtGPU[i][j] + dU4cdtGPU[i][j]) * (dt / 2.0f));

					const float PXCoef = concurrency::fast_math::fabsf(PpnGPU[i][j + 1] - (2.0f * PpnGPU[i][j]) + PpnGPU[i][j - 1]) /
						(PpnGPU[i][j + 1] + (2.0f * PpnGPU[i][j]) + PpnGPU[i][j - 1]);
					const float PYCoef = concurrency::fast_math::fabsf(PpnGPU[i - 1][j] - (2.0f * PpnGPU[i][j]) + PpnGPU[i + 1][j]) /
						(PpnGPU[i - 1][j] + (2.0f * PpnGPU[i][j]) + PpnGPU[i + 1][j]);

					const float xcoef = CxD * PXCoef;
					const float ycoef = CyD * PYCoef;

					const float Srho = CalcSrho(i, j, U1nGPU, xcoef, ycoef);
					const float Srhou = CalcSrhou(i, j, U2nGPU, xcoef, ycoef);
					const float Srhov = CalcSrhov(i, j, U3nGPU, xcoef, ycoef);
					const float SEt = CalcSEt(i, j, U4nGPU, xcoef, ycoef);

					U1pnGPU[i][j] += Srho;
					U2pnGPU[i][j] += Srhou;
					U3pnGPU[i][j] += Srhov;
					U4pnGPU[i][j] += SEt;

					unGPU[i][j] = U2nGPU[i][j] / U1nGPU[i][j];
					vnGPU[i][j] = U3nGPU[i][j] / U1nGPU[i][j];
					TnGPU[i][j] = ((U4nGPU[i][j] / U1nGPU[i][j]) - ((unGPU[i][j] * unGPU[i][j]) + (vnGPU[i][j] * vnGPU[i][j]) / 2.0f)) / CvD;
					PnGPU[i][j] = U1nGPU[i][j] * RD * TnGPU[i][j];

					munGPU[i][j] = CalcViscSutherland(muRefD, tempRefD, TnGPU[i][j]);
					knGPU[i][j] = Calckk(munGPU[i][j], CpD, prandtlNumberD);
				}
			);
			//for (unsigned int i = (NPOINTSY / 2) - 1; i < (NPOINTSY / 2) + 1; i++)
			//{
			//	for (unsigned int j = 1; j < (NPOINTSX / 2) - (NPOINTSPX / 2); j++)
			//	{
			//		Ucorrected[j + (i*NPOINTSX)] = U[j + (i*NPOINTSX)] +
			//			((dUdtpredicted[j + (i * NPOINTSX)] + dUdtcorrected[j + (i * NPOINTSX)])*(dt / 2.0));
			//
			//		const float Srho = CalcSrho(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhou = CalcSrhou(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhov = CalcSrhov(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float SEt = CalcSEt(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		Ucorrected[j + (i*NPOINTSX)].u += Srho;
			//		Ucorrected[j + (i*NPOINTSX)].v += Srhou;
			//		Ucorrected[j + (i*NPOINTSX)].P += Srhov;
			//		Ucorrected[j + (i*NPOINTSX)].T += SEt;
			//
			//		estado[j + (i*NPOINTSX)] = CalcState(Ucorrected[j + (i*NPOINTSX)]);
			//
			//		mu[j + (i*NPOINTSX)] = CalcViscSutherland(muRef, tempRef, estado[j + (i*NPOINTSX)].T);
			//		k[j + (i*NPOINTSX)] = Calck(mu[j + (i*NPOINTSX)]);
			//		
			//		if (estado[j + (i*NPOINTSX)].T < 0.0)
			//		{
			//			std::cout << "Erro i: " << i << " j: " << j << " on corrected values" << std::endl;
			//		}
			//
			//	}
			//}

			// RIGHT PART
			concurrency::parallel_for_each(rightComputeDomain,
				[=](concurrency::index<2> idx) restrict(amp)
				{
					const int i = idx[0] + (NPOINTSYD / 2) - 1;
					const int j = idx[1] + (NPOINTSXD / 2) + (NPOINTSPXD / 2);

					U1nGPU[i][j] = U1nGPU[i][j] + ((dU1pdtGPU[i][j] + dU1cdtGPU[i][j]) * (dt / 2.0f));
					U2nGPU[i][j] = U2nGPU[i][j] + ((dU2pdtGPU[i][j] + dU2cdtGPU[i][j]) * (dt / 2.0f));
					U3nGPU[i][j] = U3nGPU[i][j] + ((dU3pdtGPU[i][j] + dU3cdtGPU[i][j]) * (dt / 2.0f));
					U4nGPU[i][j] = U4nGPU[i][j] + ((dU4pdtGPU[i][j] + dU4cdtGPU[i][j]) * (dt / 2.0f));

					const float PXCoef = concurrency::fast_math::fabsf(PpnGPU[i][j + 1] - (2.0f * PpnGPU[i][j]) + PpnGPU[i][j - 1]) /
						(PpnGPU[i][j + 1] + (2.0f * PpnGPU[i][j]) + PpnGPU[i][j - 1]);
					const float PYCoef = concurrency::fast_math::fabsf(PpnGPU[i - 1][j] - (2.0f * PpnGPU[i][j]) + PpnGPU[i + 1][j]) /
						(PpnGPU[i - 1][j] + (2.0f * PpnGPU[i][j]) + PpnGPU[i + 1][j]);

					const float xcoef = CxD * PXCoef;
					const float ycoef = CyD * PYCoef;

					const float Srho = CalcSrho(i, j, U1nGPU, xcoef, ycoef);
					const float Srhou = CalcSrhou(i, j, U2nGPU, xcoef, ycoef);
					const float Srhov = CalcSrhov(i, j, U3nGPU, xcoef, ycoef);
					const float SEt = CalcSEt(i, j, U4nGPU, xcoef, ycoef);

					U1pnGPU[i][j] += Srho;
					U2pnGPU[i][j] += Srhou;
					U3pnGPU[i][j] += Srhov;
					U4pnGPU[i][j] += SEt;

					unGPU[i][j] = U2nGPU[i][j] / U1nGPU[i][j];
					vnGPU[i][j] = U3nGPU[i][j] / U1nGPU[i][j];
					TnGPU[i][j] = ((U4nGPU[i][j] / U1nGPU[i][j]) - ((unGPU[i][j] * unGPU[i][j]) + (vnGPU[i][j] * vnGPU[i][j]) / 2.0f)) / CvD;
					PnGPU[i][j] = U1nGPU[i][j] * RD * TnGPU[i][j];

					munGPU[i][j] = CalcViscSutherland(muRefD, tempRefD, TnGPU[i][j]);
					knGPU[i][j] = Calckk(munGPU[i][j], CpD, prandtlNumberD);
				}
			);

			U1pnGPU.synchronize();
			U2pnGPU.synchronize();
			U3pnGPU.synchronize();
			U4pnGPU.synchronize();

			//for (unsigned int i = (NPOINTSY / 2) - 1; i < (NPOINTSY / 2) + 1; i++)
			//{
			//	for (unsigned int j = (NPOINTSX / 2) + (NPOINTSPX / 2); j < NPOINTSX - 1; j++)
			//	{
			//		Ucorrected[j + (i*NPOINTSX)] = U[j + (i*NPOINTSX)] +
			//			((dUdtpredicted[j + (i * NPOINTSX)] + dUdtcorrected[j + (i * NPOINTSX)])*(dt / 2.0));
			//
			//		const float Srho = CalcSrho(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhou = CalcSrhou(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float Srhov = CalcSrhov(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//		const float SEt = CalcSEt(Upredicted, estadoPredicted, i, j, NPOINTSX, NPOINTSY);
			//
			//		Ucorrected[j + (i*NPOINTSX)].u += Srho;
			//		Ucorrected[j + (i*NPOINTSX)].v += Srhou;
			//		Ucorrected[j + (i*NPOINTSX)].P += Srhov;
			//		Ucorrected[j + (i*NPOINTSX)].T += SEt;
			//
			//		estado[j + (i*NPOINTSX)] = CalcState(Ucorrected[j + (i*NPOINTSX)]);
			//
			//		mu[j + (i*NPOINTSX)] = CalcViscSutherland(muRef, tempRef, estado[j + (i*NPOINTSX)].T);
			//		k[j + (i*NPOINTSX)] = Calck(mu[j + (i*NPOINTSX)]);
			//		
			//		if (estado[j + (i*NPOINTSX)].T < 0.0)
			//		{
			//			std::cout << "Erro i: " << i << " j: " << j << " on corrected values" << std::endl;
			//		}
			//
			//	}
			//}

			unGPU.synchronize();
			vnGPU.synchronize();
			PnGPU.synchronize();
			TnGPU.synchronize();

			CalcExtrapolationsBC(unGPU, vnGPU, PnGPU, TnGPU);




			totalTime += dt;

			if (it % saveEvery == 0)
			{
				std::vector<StateData> estado;

				std::cout << std::endl;
				std::cout << "Iteration number: " << it << std::endl;
				std::cout << "dt: " << dt << std::endl;

				std::cout << "Total time elapsed: " << totalTime << std::endl;

				//std::cout << "Min mesh Reynolds: " << minMeshReynolds << std::endl;

				const int NPOINTSXINT = static_cast<int>(NPOINTSXD);
				const int NPOINTSYINT = static_cast<int>(NPOINTSYD);

				std::string nend = std::to_string(it) + ".datf";

				std::ofstream arqstate("out/state" + nend, std::ios::binary);

				arqstate << "Total time: ;" << std::to_string(totalTime * 1e6) << " us\n";
				arqstate << "%DAT#";
				arqstate.write(reinterpret_cast<char const*>(&NPOINTSXINT), sizeof(int));
				arqstate << ";";
				arqstate.write(reinterpret_cast<char const*>(&NPOINTSYINT), sizeof(int));
				arqstate << "\n";

				for (int i = 0; i < NPOINTSYD; i++)
				{
					for (int j = 0; j < NPOINTSXD; j++)
					{
						estado.push_back({ unGPU[i][j], vnGPU[i][j], PnGPU[i][j], TnGPU[i][j] });
						//std::cout << "i: " << i << " j: " << j << " unGPU: " << unGPU[i][j] << std::endl;
					}
				}


				arqstate.write(reinterpret_cast<char const*>(estado.data()), sizeof(StateData)* size_t(NPOINTSXD * NPOINTSYD));

				arqstate.close();

				std::cout << "Salvo arquivo " << "out/state" + nend << std::endl;

				if (SALVARASCII)
				{
					std::ofstream arqu("out/u/u" + nend);

					arqu << "Total time: ;" << std::to_string(totalTime * 1e6) << " us\n";

					for (unsigned int i = 0; i < NPOINTSYD; i++)
					{
						for (unsigned int j = 0; j < NPOINTSXD; j++)
						{
							arqu << estado[j + (i * NPOINTSXD)].u << ";";

						}
						arqu << "\n";
					}


					arqu.close();

					std::ofstream arqv("out/v/v" + nend);

					arqv << "Total time: ;" << std::to_string(totalTime * 1e6) << " us\n";

					for (unsigned int i = 0; i < NPOINTSYD; i++)
					{
						for (unsigned int j = 0; j < NPOINTSXD; j++)
						{
							arqv << estado[j + (i * NPOINTSXD)].v << ";";

						}
						arqv << "\n";
					}

					arqv.close();

					std::ofstream arqp("out/p/P" + nend);

					arqp << "Total time: ;" << std::to_string(totalTime * 1e6) << " us\n";

					for (unsigned int i = 0; i < NPOINTSYD; i++)
					{
						for (unsigned int j = 0; j < NPOINTSXD; j++)
						{
							arqp << estado[j + (i * NPOINTSXD)].P << ";";

						}
						arqp << "\n";
					}

					arqp.close();

					std::ofstream arqT("out/t/T" + nend);

					arqT << "Total time: ;" << std::to_string(totalTime * 1e6) << " us\n";

					for (unsigned int i = 0; i < NPOINTSYD; i++)
					{
						for (unsigned int j = 0; j < NPOINTSXD; j++)
						{
							arqT << estado[j + (i * NPOINTSXD)].T << ";";

						}
						arqT << "\n";
					}

					arqT.close();
				}
			}

		}
	}
	catch (...)
	{
		std::cout << "thrown" << std::endl;
	}

	

	

	return 0;
}

// TODO: Conferir Equações abaixo.
// TODO²: Terminar de paralelizar tudo.

float CalcDensity(float P, float R, float T) restrict(amp)
{
	return P / (R*T);
}
float CalcViscSutherland(float mu0, float T0, float T) restrict(amp)
{
	return mu0 * concurrency::fast_math::pow(T / T0, 1.5f)*((T0 + 110.0f) / (T + 110.0f));
}

float Calckk(float mu, float Cp, float prandtlNumber) restrict(amp)
{
	return mu * Cp / prandtlNumber;
}

float CalcEt(float rho, float Cv, float T, float Vquad) restrict(amp)
{
	return rho * ((Cv * T) + (Vquad / 2.0f));
}

float CalcDensity(float P, float R, float T)
{
	return P / (R * T);
}
float CalcViscSutherland(float mu0, float T0, float T)
{
	return mu0 * concurrency::fast_math::pow(T / T0, 1.5f) * ((T0 + 110.0f) / (T + 110.0f));
}

float Calck(float mu, float Cp, float prandtlNumber)
{
	return mu * CpD / prandtlNumberD;
}

float CalcEt(float rho, float Cv, float T, float Vquad)
{
	return rho * ((CvD * T) + (Vquad / 2.0f));
}

void CalcExtrapolationsBC(concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn,
	concurrency::array_view<float, 2> Pn, concurrency::array_view<float, 2> Tn)
{	

	concurrency::extent<1> computeDomainWallsJ = concurrency::extent<1>(NPOINTSPXD);

	concurrency::parallel_for_each(computeDomainWallsJ,
		[=](concurrency::index<1> idx) restrict(amp)
		{
			const int j = idx[0] + (NPOINTSXD / 2) - (NPOINTSPXD / 2);

			Pn[(NPOINTSYD / 2) - 1][j] = Pn[(NPOINTSYD / 2) - 2][j];

			Pn[(NPOINTSYD / 2)][j] = Pn[(NPOINTSYD / 2) + 1][j];
		}
	);
	
	//for (unsigned int j = (NPOINTSXD / 2) - (NPOINTSPXD/2); j < (NPOINTSXD / 2) + (NPOINTSPXD/2); j++)
	//{
	//	// CHANGED EXTRAPOLATION TO ZERO GRADIENT!
	//
	//	Pn[(NPOINTSYD / 2) - 1][j] = Pn[(NPOINTSYD / 2) - 2][j];
	//
	//	Pn[(NPOINTSYD / 2)][j] = Pn[(NPOINTSYD / 2) + 1][j];
	//
	//	//field[j + (((NPOINTSYD / 2) - 1)*NPOINTSXD)].P = field[j + (((NPOINTSYD / 2) - 2)*NPOINTSXD)].P;//pextrExtradorso;
	//
	//	//field[j + (((NPOINTSYD / 2))*NPOINTSXD)].P = field[j + (((NPOINTSYD / 2) + 1)*NPOINTSXD)].P;//pextrIntradorso;
	//}
	//
	// CASE 4 (outflow except inflow)


	concurrency::extent<1> computeDomainOutflow = concurrency::extent<1>(NPOINTSYD - 2);

	concurrency::parallel_for_each(computeDomainOutflow,
		[=](concurrency::index<1> idx) restrict(amp)
		{
			const int i = idx[0] + 1;

			un[i][NPOINTSXD - 1] = un[i][NPOINTSXD - 2];
			vn[i][NPOINTSXD - 1] = vn[i][NPOINTSXD - 2];
			Tn[i][NPOINTSXD - 1] = Tn[i][NPOINTSXD - 2];
		});

	//for (unsigned int i = 1; i < NPOINTSYD - 1; i++)
	//{
	//	
	//	un[i][NPOINTSXD - 1] = un[i][NPOINTSXD - 2];
	//	vn[i][NPOINTSXD - 1] = vn[i][NPOINTSXD - 2];
	//	Tn[i][NPOINTSXD - 1] = Tn[i][NPOINTSXD - 2];
	//
	//	//field[NPOINTSXD - 1 + (i*NPOINTSXD)].u = field[NPOINTSXD - 2 + (i*NPOINTSXD)].u;//uextr;
	//	//field[NPOINTSXD - 1 + (i*NPOINTSXD)].v = field[NPOINTSXD - 2 + (i*NPOINTSXD)].v;//vextr;
	//	//field[NPOINTSXD - 1 + (i*NPOINTSXD)].T = field[NPOINTSXD - 2 + (i*NPOINTSXD)].T;//textr;
	//}
}

float GetE1Predictor(int i, int j, concurrency::array_view<float, 2> rhon,
	concurrency::array_view<float, 2> un) restrict(amp)
{
	const float uij = un[i][j];
	const float rhoij = rhon[i][j];
	const float res = uij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);

	return res;
}

float GetE2Predictor(int i, int j, concurrency::array_view<float, 2> rhon,
	concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> Pn, concurrency::array_view<float, 2> mun) restrict(amp)
{
	const float uij = un[i][j];
	const float rhoij = rhon[i][j];
	const float Pij = Pn[i][j];
	float res = uij * uij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);

	res += Pij;
	res -= CalcTauXXPredictor(i, j, un, vn, mun);

	return res;
}


float GetE3Predictor(int i, int j, concurrency::array_view<float, 2> rhon,
	concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun)restrict(amp)
{
	const float vij = vn[i][j];
	const float uij = un[i][j];
	const float rhoij = rhon[i][j];
	float res = vij * uij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);

	res -= CalcTauXYPredictorE(i, j, un, vn, mun);

	return res;
}

float GetE4Predictor(int i, int j, concurrency::array_view<float, 2> Etn,
	concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun,
	concurrency::array_view<float, 2> Tn, concurrency::array_view<float, 2> Kn, concurrency::array_view<float, 2> Pn)restrict(amp)
{
	const float uij = un[i][j];
	const float vij = vn[i][j];
	const float Etij = Etn[i][j];
	const float Pij = Pn[i][j];

	const float Vquad = (uij* uij) + (vij* vij);
	float res = uij * Etij;
	res += Pij*uij;
	res -= uij * CalcTauXXPredictor(i, j, un,vn, mun);
	res -= vij * CalcTauXYPredictorE(i, j, un, vn, mun);
	res += CalcQxPredictor(i, j, Tn, Kn);

	return res;
}


float GetF1Predictor(int i, int j, concurrency::array_view<float, 2> rhon,
	concurrency::array_view<float, 2> vn)restrict(amp)
{
	const float vij = vn[i][j];
	const float rhoij = rhon[i][j];
	const float res = vij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);

	return res;
}

float GetF2Predictor(int i, int j, concurrency::array_view<float, 2> rhon,
	concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun) restrict(amp)
{
	const float uij = un[i][j];
	const float vij = vn[i][j];
	const float rhoij = rhon[i][j];
	float res = uij*vij * rhoij;// CalcDensity(fieldij.P, R, fieldij.T);

	res -= CalcTauXYPredictorF(i, j, un, vn, mun);

	return res;
}

float GetF3Predictor(int i, int j, concurrency::array_view<float, 2> rhon,
	concurrency::array_view<float, 2> un,
	concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun,
	concurrency::array_view<float, 2> Pn) restrict(amp) 
{
	const float vij = vn[i][j];
	const float rhoij = rhon[i][j];
	const float Pij = Pn[i][j];
	float res = vij*vij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);
	res += Pij;

	res -= CalcTauYYPredictor(i, j, un, vn, mun);

	return res;
}

float GetF4Predictor(int i, int j, concurrency::array_view<float, 2> Etn,
	concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun,
	concurrency::array_view<float, 2> Tn, concurrency::array_view<float, 2> Kn, concurrency::array_view<float, 2> Pn) restrict(amp) 
{
	const float uij = un[i][j];
	const float vij = vn[i][j];
	const float Etij = Etn[i][j];
	const float Pij = Pn[i][j];

	//const float Vquad = (uij*uij) + (vij*vij);

	float res = vij * Etij;
	res += Pij*vij;
	res -= uij * CalcTauXYPredictorF(i, j, un, vn, mun);
	res -= vij * CalcTauYYPredictor(i, j, un, vn, mun);
	res += CalcQyPredictor(i, j, Tn, Kn);

	return res;
}

float CalcTauXXPredictor(int i, int j, concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn,
	concurrency::array_view<float, 2> mun) restrict(amp)
{
	const float uij = un[i][j];
	const float uimj = un[i][j - 1];
	
	const float vijm = vn[i + 1][j];
	const float vijp = vn[i - 1][j];

	const float muij = mun[i][j];//CalcViscSutherland(muRef, tempRef, fieldij.T);

	const float res = (2.0f / 3.0f) * muij * (
		((2.0f / dxD)*(uij - uimj)) - (((vijp - vijm) / (2.0f*dyD)))
		);

	return res;
}

float CalcTauYYPredictor(int i, int j, concurrency::array_view<float, 2> un, concurrency::array_view<float, 2> vn,
	concurrency::array_view<float, 2> mun) restrict(amp)
{
	const float vij = vn[i][j];
	const float vijm = vn[i + 1][j];

	const float uipj = un[i][j + 1];
	const float uimj = un[i][j - 1];

	const float muij = mun[i][j];// = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const float res = (2.0f / 3.0f) * muij * (
		((2.0f / dyD)*(vij - vijm)) - (((uipj - uimj) / (2.0f*dxD)))
		);

	return res;
}

float CalcTauXYPredictorE(int i, int j, concurrency::array_view<float, 2> un,
	concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun) restrict(amp)
{
	const float vij = vn[i][j];
	const float vimj = vn[i][j - 1];

	const float uijp = un[i - 1][j];
	const float uijm = un[i + 1][j];

	const float muij = mun[i][j];// = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const float res = muij * (
		((vij - vimj) / dxD) + (((uijp - uijm) / (2.0f*dyD)))
		);

	return res;
}

float CalcTauXYPredictorF(int i, int j, concurrency::array_view<float, 2> un,
	concurrency::array_view<float, 2> vn, concurrency::array_view<float, 2> mun) restrict(amp)
{
	const float uij = un[i][j];
	const float uijm = un[i + 1][j];

	const float vipj = vn[i][j + 1];
	const float vimj = vn[i][j - 1];
	
	const float muij = mun[i][j];// = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const float res = muij * (
		((uij - uijm) / dyD) + (((vipj - vimj) / (2.0f*dxD)))
		);

	return res;
}

float CalcQxPredictor(int i, int j, concurrency::array_view<float, 2> Tn,
	concurrency::array_view<float, 2> Kn) restrict(amp)
{
	const float Tij = Tn[i][j];
	const float Timj = Tn[i][j - 1];

	//const float muij = mu[j + (i*NPOINTSX)];// = CalcViscSutherland(muRef, tempRef, fieldij.T);
	const float kij = Kn[i][j];// = Calck(mu);

	const float res = -kij * (Tij - Timj) / dxD;

	return res;
}

float CalcQyPredictor(int i, int j, concurrency::array_view<float, 2> Tn,
	concurrency::array_view<float, 2> Kn) restrict(amp)
{
	const float Tij = Tn[i][j];
	const float Tijm = Tn[i + 1][j];
	
	const float kij = Kn[i][j];// = Calck(mu);

	const float res = -kij * (Tij - Tijm) / dyD;

	return res;
}

float GetE1Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn,
	concurrency::array_view<float, 2> uPredn) restrict(amp)
{
	const float uij = uPredn[i][j];
	const float rhoij = rhoPredn[i][j];

	const float res = uij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);

	return res;
}

float GetE2Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn,
	concurrency::array_view<float, 2> uPredn, concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> PPredn, concurrency::array_view<float, 2> muPredn) restrict(amp)

{
	const float uij = uPredn[i][j];
	const float rhoij = rhoPredn[i][j];
	const float Pij = PPredn[i][j];

	float res = uij * uij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);

	res += Pij;
	res -= CalcTauXXCorrector(i, j, uPredn, vPredn, muPredn);

	return res;
}

float GetE3Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn,
	concurrency::array_view<float, 2> uPredn, concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn)restrict(amp)
{
	const float vij = vPredn[i][j];
	const float uij = uPredn[i][j];
	const float rhoij = rhoPredn[i][j];

	float res = vij * uij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);

	res -= CalcTauXYCorrectorE(i, j, uPredn, vPredn, muPredn);

	return res;
}

float GetE4Corrector(int i, int j, concurrency::array_view<float, 2> EtPredn,
	concurrency::array_view<float, 2> uPredn, concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn,
	concurrency::array_view<float, 2> TPredn, concurrency::array_view<float, 2> KPredn, concurrency::array_view<float, 2> PPredn)restrict(amp)
{
	const float uij = uPredn[i][j];
	const float vij = vPredn[i][j];
	const float Etij = EtPredn[i][j];
	const float Pij = PPredn[i][j];

	const float Vquad = (uij*uij) + (vij*vij);
	float res = uij * Etij;
	res += Pij*uij;
	res -= uij * CalcTauXXCorrector(i, j, uPredn, vPredn, muPredn);
	res -= vij * CalcTauXYCorrectorE(i, j, uPredn, vPredn, muPredn);
	res += CalcQxCorrector(i, j, TPredn, KPredn);

	return res;
}

float GetF1Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn,
	concurrency::array_view<float, 2> vPredn)restrict(amp)
{
	const float vij = vPredn[i][j];
	const float rhoij = rhoPredn[i][j];


	const float res = vij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);

	return res;
}

float GetF2Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn,
	concurrency::array_view<float, 2> uPredn, concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn) restrict(amp)
{
	const float uij = uPredn[i][j];
	const float vij = vPredn[i][j];
	const float rhoij = rhoPredn[i][j];

	float res = uij*vij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);

	res -= CalcTauXYCorrectorF(i, j, uPredn, vPredn, muPredn);

	return res;
}

float GetF3Corrector(int i, int j, concurrency::array_view<float, 2> rhoPredn, concurrency::array_view<float, 2> uPredn,
	concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn,
	concurrency::array_view<float, 2> PPredn) restrict(amp)
{
	const float vij = vPredn[i][j];
	const float rhoij = rhoPredn[i][j];
	const float Pij = PPredn[i][j];
	
	float res = vij*vij * rhoij;//CalcDensity(fieldij.P, R, fieldij.T);
	res += Pij;

	res -= CalcTauYYCorrector(i, j, uPredn, vPredn, muPredn);

	return res;
}

float GetF4Corrector(int i, int j, concurrency::array_view<float, 2> EtPredn,
	concurrency::array_view<float, 2> uPredn, concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn,
	concurrency::array_view<float, 2> TPredn, concurrency::array_view<float, 2> KPredn, concurrency::array_view<float, 2> PPredn) restrict(amp)
{
	const float uij = uPredn[i][j];
	const float vij = vPredn[i][j];
	const float Etij = EtPredn[i][j];

	const float Pij = PPredn[i][j];

	const float Vquad = (uij*uij) + (vij*vij);

	float res = vij * Etij;
	res += Pij*vij;
	res -= uij * CalcTauXYCorrectorF(i, j, uPredn, vPredn, muPredn);
	res -= vij * CalcTauYYCorrector(i, j, uPredn, vPredn, muPredn);
	res += CalcQyCorrector(i, j, TPredn, KPredn);

	return res;
}


float CalcTauXXCorrector(int i, int j, concurrency::array_view<float, 2> uPredn,
	concurrency::array_view<float, 2> vPredn,
	concurrency::array_view<float, 2> muPredn) restrict(amp)
{
	const float uipj = uPredn[i][j + 1];
	const float uij = uPredn[i][j];

	const float vijp = vPredn[i - 1][j];
	const float vijm = vPredn[i + 1][j];

	const float muij = muPredn[i][j];// = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const float res = (2.0f / 3.0f) * muij * (
		((2.0f / dxD)*(uipj - uij)) - (((vijp- vijm) / (2.0f*dyD)))
		);

	return res;
}

float CalcTauYYCorrector(int i, int j, concurrency::array_view<float, 2> uPredn,
	concurrency::array_view<float, 2> vPredn,
	concurrency::array_view<float, 2> muPredn) restrict(amp)
{
	const float vijp = vPredn[i - 1][j];
	const float vij = vPredn[i][j];

	const float uipj = uPredn[i][j + 1];
	const float uimj = uPredn[i][j - 1];

	const float muij = muPredn[i][j];// = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const float res = (2.0f / 3.0f) * muij * (
		((2.0f / dyD)*(vijp - vij)) - (((uipj - uimj) / (2.0f*dxD)))
		);

	return res;
}

float CalcTauXYCorrectorE(int i, int j, concurrency::array_view<float, 2> uPredn,
	concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn) restrict(amp)
{
	const float vipj = vPredn[i][j+1];
	const float vij = vPredn[i][j];

	const float uijp = uPredn[i - 1][j];
	const float uijm = uPredn[i + 1][j];

	const float muij = muPredn[i][j];// = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const float res = muij * (
		((vipj - vij) / dxD) + ((uijp- uijm) / (2.0f*dyD))
		);

	return res;
}

// i j in "array notation" (change in row and collumn respectively)
float CalcTauXYCorrectorF(int i, int j, concurrency::array_view<float, 2> uPredn,
	concurrency::array_view<float, 2> vPredn, concurrency::array_view<float, 2> muPredn) restrict(amp)
{
	const float uijp = uPredn[i - 1][j];
	const float uij = uPredn[i][j];

	const float vipj = vPredn[i][j + 1];
	const float vimj = vPredn[i][j - 1];

	const float muij = muPredn[i][j];// = CalcViscSutherland(muRef, tempRef, fieldij.T);

	const float res = muij * (
		((uijp - uij) / dyD) + ((vipj - vimj) / (2.0f*dxD))
		);

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
float CalcQxCorrector(int i, int j, concurrency::array_view<float, 2> Tn, concurrency::array_view<float, 2> Kn) restrict(amp)
{
	const float& Tnij =  Tn[i][j];
	const float& Tnipj = Tn[i][j + 1];

	//const float mu = CalcViscSutherland(muRef, tempRef, fieldij.T);
	const float kij = Kn[i][j];

	const float res = -kij * (Tnipj - Tnij) / dxD;

	return res;
}
// i j in "array notation" (change in row and collumn respectively)
float CalcQyCorrector(int i, int j, concurrency::array_view<float, 2> TPredn,
	concurrency::array_view<float, 2> KPredn) restrict(amp)
{
	const float Tij = TPredn[i][j];
	const float Tijp = TPredn[i - 1][j];

	//const float mu = CalcViscSutherland(muRef, tempRef, fieldij.T);
	const float kij = KPredn[i][j];// = Calck(mu);

	const float res = -kij * (Tijp - Tij) / dyD;

	return res;
}

StateData CalcState(const StateData& U)
{
	StateData estado;
	if (U.u == 0.0f)
	{
		assert(0);
	}
	estado.u = U.v / U.u;
	estado.v = U.P / U.u;

	const float e = (U.T / U.u) - (((estado.u*estado.u) + (estado.v * estado.v)) / 2.0f);
	if (e < 0.0f)
	{
		assert(0);

		MessageBox(NULL, L"negative Temperature :x", L"ERROR", 0);
	}
	estado.T = e / CvD;
	estado.P = U.u * RD * estado.T;

	return estado;
}

StateData CalcU(const StateData& estado)
{
	StateData U;
	const float Vquad = (estado.u*estado.u) + (estado.v * estado.v);
	if (estado.T == 0.0f)
	{
		assert(0);
	}
	U.u = estado.P / (RD*estado.T);
	U.v = U.u * estado.u;
	U.P = U.u* estado.v;
	U.T = U.u*((CvD * estado.T) + (Vquad / 2.0f));

	return U;
}


float CalcSrho(int i, int j, concurrency::array_view<float, 2> rhon,
	float xcoef, float ycoef) restrict(amp)
{
	const float rhoij = rhon[i][j];
	const float rhoipj = rhon[i][j + 1];
	const float rhoimj = rhon[i][j - 1];
	const float rhoijp = rhon[i - 1][j];
	const float rhoijm = rhon[i + 1][j];

	return (xcoef*(rhoipj - (2.0f*rhoij) + rhoimj)) +
		(ycoef*(rhoijp-(2.0f*rhoij)+rhoijm));
}


float CalcSrhou(int i, int j, concurrency::array_view<float, 2> U2n,
	float xcoef, float ycoef) restrict(amp) 
{

	const float& rhouij = U2n[i][j];
	const float& rhouipj = U2n[i][j + 1];
	const float& rhouimj = U2n[i][j - 1];
	const float& rhouijp = U2n[i - 1][j];
	const float& rhouijm = U2n[i + 1][j];

	return (xcoef*(rhouipj - (2.0f*rhouij) + rhouimj)) +
		(ycoef*(rhouijp - (2.0f*rhouij) + rhouijm));
}


float CalcSrhov(int i, int j, concurrency::array_view<float, 2> U3n,
	float xcoef, float ycoef) restrict(amp)
{
	const float& rhovij = U3n[i][j];
	const float& rhovipj = U3n[i][j + 1];
	const float& rhovimj = U3n[i][j - 1];
	const float& rhovijp = U3n[i - 1][j];
	const float& rhovijm = U3n[i + 1][j];

	return (xcoef*(rhovipj - (2.0f*rhovij) + rhovimj)) +
		(ycoef*(rhovijp - (2.0f*rhovij) + rhovijm));
}


float CalcSEt(int i, int j, concurrency::array_view<float, 2> Etn,
	float xcoef, float ycoef) restrict(amp)
{
	const float& Etij = Etn[i][j];
	const float& Etipj = Etn[i][j+1];
	const float& Etimj = Etn[i][j-1];
	const float& Etijp = Etn[i-1][j];
	const float& Etijm = Etn[i+1][j];

	return (xcoef*(Etipj - (2.0f*Etij) + Etimj)) +
		(ycoef*(Etijp - (2.0f*Etij) + Etijm));
}