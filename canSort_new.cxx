/* This code reads list-mode data files in the CANDLE format and writes spectra
 * in a ROOT Format. The first "N" spectra will be the RAW spectra,
 * corresponding to the number of parameters (N) in the experiment.
 *
 * If the calibration parameters are provided, then the next "N" spectra will be
 * calibrated to 0.5 keV/ch. If you don't have the calibration parameters, then
 * atleast a file with dummy parameters (0 1 0) should be provided.
 *
 * This is followed by generation of the aligned TDC and Add-Back spectra.
 *
 * RADWARE compatible matrix WITH or WITHOUT the TAC condition can be also made.
 *
 * A provision to make a dump file for CUBE sorting is also incorporated.
 * The dump file is used as an input to "incub8r_inga".
 *
 * Compile this using " g++ -std=c++11 -o canSort1 canSort1.cxx `root-config --cflags --glibs`" command.
 * Add "-D_ZERO_SUPPRESSION_LEVEL_=0" for a "non-zero-suppressed" tree.
 * Add "-D_ZERO_SUPPRESSION_LEVEL_=2" if you also want to suppress events with "clov_mult = 0".
 * Note: "-D_ZERO_SUPPRESSION_LEVEL_=1" is the default "zero-suppressed" tree.
 *
 * Written by Ajay Y. Deo, IIT Roorkee.
 * Date: August 25, 2016
 *
 * November 2016: 	Modified to generate PDCO matrices. (AYDeo)
 *
 *
 * December 1, 2016:	Modified to include "time condition" while generating
 *			the dump file for CUBE sorting. (AYDeo)
 *
 * September 2017:      Resolved a minor(?) issue in PDCO Matrices. (AYDeo)
 *
 * February 2018:	Some minor modifications
 *
 */


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <bitset>
#include <algorithm>    // std::all_of
#include <array>        // std::array
#include <cmath>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TSpectrum2.h"
#include "TSpectrum.h"

#include "COMMON.h"

#define	ONE_MINUTE	60000.0

using namespace std;

#if !defined(_ZERO_SUPPRESSION_LEVEL_)
// default setting: "variable size clover arrays"
#define _ZERO_SUPPRESSION_LEVEL_ 1
#endif /* !defined(_ZERO_SUPPRESSION_LEVEL_) */

TH1I *his[100];
TH1I *his_cal[100];
TH1I *his_AB[20];
TH1I *his_TDC[20];
TH1I *his_AB1[20];
TH1I *his_AB2[20];
TH1I *his_CVT[20];
TH1I *his_AFC[20];
TH1I *his_TAC;
TH2I *his2D;
TH2I *his_par;
TH2I *his_per;
TH2I *his_GT;
TH1 *bckgnd ;
TH1 *his2D_wb ;

TTree *t = 0;

TFile *f=NULL;
//FILE *MatFile = NULL;
//FILE *MatFile1 = NULL;
//FILE *MatFile2 = NULL;
FILE *DumpFile = NULL;

const int No_Clover = 16;

bitset<No_Clover> cl_mult;

char nm[30];

int NoPara, clov_mult, ZeroFold = 0, OneFold = 0, TwoFold = 0;
int ThreeFold = 0, FourFold = 0, FiveFold = 0, SixFold = 0, SixteenFold = 0;
int MatCase, Flag2D = 0, Flag3D = 0, FlagTAC = 0, TACL, TACH;

//unsigned short *Matrix, *Matrix_par, *Matrix_per, *Matrix_GT;
unsigned long x, y, z, MatDim = 4096, MatSize, Time_Diff;
unsigned long Eclab[16] = {0}, CLT[16] = {0}, energy_cal[500];
unsigned long E_par1[100], E_par2[100], E_par[100], E_per1[100], E_per2[100], E_per[100], EM_par[100], EM_per[100], CLT_M[20] = {0};

int BufSize, CubBufEvent = 0, BufNum = 1;
unsigned short CubBuf[8192] = {0};
unsigned short CubT[16] = {0};
static long sbuf[32768] = {0};
//static short sbuf[32768]; //AYDeo: "static short" was acting up; hence changed to "static long"
int valid_3F_events = 0;
bool  GammaTime = false;
/******************/
void AddBack();
void CountVsTime();
void MatGG();
void MatDCO();
void MatPDCO();
void ThreeD();
void MatED();
//void rembckgnd(TH1*);
//void rembckgnd(TH1* back);
/******************/

void expand (datapackptr inbuf, datapackptr outbuf);

int32		debug=0;
datapack	buf, readbuf;

int32 printHeader (dataheadptr h)
{
    printf ("+++++++++++++ HEADER +++++++++++++\n");
    printf ("\tBlock type      = %d\n", h->block);
    printf ("\tUnit Size       = %d\n", h->unitsize);
    printf ("\tNumber of Units = %d\n", h->number_of_units);
    printf ("\tComp Status     = %d\n", h->compstatus);
    printf ("\tSize In Bytes   = %d\n", h->size_in_bytes);
    printf ("++++++++++++++++++++++++++++++++++\n");
    return 0;
}

int32 process_event_block (datapackptr rdbuf)
{

    expand (rdbuf, &buf);

    if (debug) {
	unsigned	i, k;
	printf ("eventSize=%d. eventsPerBlock=%d compressedSize=%d words\n",
		rdbuf->header.unitsize, rdbuf->header.number_of_units,
		rdbuf->header.size_in_bytes / 2);

	for (k=0; k<buf.header.number_of_units; ++k) {
	    int32 offset = k * buf.header.unitsize;
	    fprintf (stdout, "%4d. ", k);
	    for (i=0; i<buf.header.unitsize; ++i)
		fprintf (stdout, "%4d ", buf.data[offset + i]);
	    	fprintf (stdout, "\n");
	}
	fprintf (stderr, "Unit#=%6d NumUnits=%6d\n",
	buf.header.unitsize, buf.header.number_of_units);
    }

    return buf.header.number_of_units;
}

int main (int argc, char *argv [])
{

cout << "\n/*************************************************/" << endl;
cout << "/**\t Written by Ajay Y. Deo \t**********/" << endl;
cout << "/**\t Dept. of Physics, IIT Roorkee \t**********/" << endl;
cout << "/*************************************************/" << endl;



time_t StartTimer, EndTimer;
double seconds;
time(&StartTimer);
FILE 		*ifp;
char		filename [512];
unsigned long	num_events = 0, num_ev_blocks = 0;
int32		prevBlock = -1, singID = 0, eventSize = 0, blockSize = 0;
ratepack	rate;
double		icr, dcr;

string calfile;
char ans[3], hisfile[512], dumpfilename[80];
bool first_event = true, first_file = true;
int parameters = 500, flag;
float coffset[500] = {0.0}, slope[500] = {0.5}, quad[500] = {0.0};

unsigned long energy[500];

    memset (&rate, 0, sizeof (ratepack));

    if (argc >= 3) {
	strcpy (filename, argv [1]);
	strcpy (hisfile, argv [2]);
	if (argc > 3) debug = 1;
    } else {
	fprintf (stderr, "Usage: canSort listmode_filename o/p_his.root_file_name\n");
	exit (1);
    }

/* 	for(int i=0;i<500;i++)
    {     snprintf(nm, sizeof(nm), "parameter%03i", i);
           his[i]=new TH1I(nm,nm,16384,0.,16384.);} */

// 101 = 100 * ROOT::kZLIB + 1, 201 = 100 * ROOT::kLZMA + 1
f = new TFile(hisfile, "RECREATE", "", 101);

t = new TTree("canSort", "canSort tree");
t->SetAutoFlush(1000000);
ROOT::TIOFeatures features;
features.Set(ROOT::Experimental::EIOFeatures::kGenerateOffsetMap);
t->SetIOFeatures(features);

// note: "_No_Clover_" is a "variable size" of some arrays so it MUST be an "Int_t"
Int_t _No_Clover_;

#if 1 /* 0 or 1 */

Int_t _clov_mult_, _Clover_[(No_Clover)];
ULong64_t _CLT_[(No_Clover)], _Eclab_[(No_Clover)];
t->Branch("clov_mult", &_clov_mult_, "clov_mult/I"); // a 32 bit signed integer (Int_t)
t->SetBasketSize("clov_mult", 16000);
t->Branch("No_Clover", &_No_Clover_, "No_Clover/I"); // a 32 bit signed integer (Int_t)
t->SetBasketSize("No_Clover", 16000);
t->Branch("Clover", _Clover_, "Clover[No_Clover]/I"); // 32 bit signed integers (Int_t)
t->SetBasketSize("Clover", 16000);
t->Branch("CLT", _CLT_, "CLT[No_Clover]/l"); // 64 bit unsigned integers (ULong64_t)
t->SetBasketSize("CLT", 16000);
t->Branch("Eclab", _Eclab_, "Eclab[No_Clover]/l"); // 64 bit unsigned integers (ULong64_t)
t->SetBasketSize("Eclab", 16000);

#elif 0 /* 0 or 1 */

Short_t _clov_mult_, _Clover_[(No_Clover)];
UShort_t _CLT_[(No_Clover)], _Eclab_[(No_Clover)];
t->Branch("clov_mult", &_clov_mult_, "clov_mult/S"); // a 16 bit signed integer (Short_t)
t->Branch("No_Clover", &_No_Clover_, "No_Clover/I"); // a 32 bit signed integer (Int_t)
t->Branch("Clover", _Clover_, "Clover[No_Clover]/S"); // 16 bit signed integers (Short_t)
t->Branch("CLT", _CLT_, "CLT[No_Clover]/s"); // 16 bit unsigned integers (UShort_t)
t->Branch("Eclab", _Eclab_, "Eclab[No_Clover]/s"); // 16 bit unsigned integers (UShort_t)

#else /* 0 or 1 */

Char_t _clov_mult_, _Clover_[(No_Clover)];
UShort_t _CLT_[(No_Clover)], _Eclab_[(No_Clover)];
t->Branch("clov_mult", &_clov_mult_, "clov_mult/B"); // an 8 bit signed integer (Char_t)
t->Branch("No_Clover", &_No_Clover_, "No_Clover/I"); // a 32 bit signed integer (Int_t)
t->Branch("Clover", _Clover_, "Clover[No_Clover]/B"); // 8 bit signed integers (Char_t)
t->Branch("CLT", _CLT_, "CLT[No_Clover]/s"); // 16 bit unsigned integers (UShort_t)
t->Branch("Eclab", _Eclab_, "Eclab[No_Clover]/s"); // 16 bit unsigned integers (UShort_t)

#endif /* 0 or 1 */

/*.......Matrix file name, memory allocation and initialization stuff begins here......*/
cout << "\nDo you want to generate Gamma-Gamma (Matrix)? (Y/N): ";
cin >> ans;
if (!strncmp(ans, "Y", 1)||!strncmp(ans, "y", 1)){
	bckgnd= new TH2I("bckgnd","bckgnd",MatDim,0.,MatDim,MatDim,0.,MatDim);
	his2D_wb= new TH2I("his2D_wb","his2D_wb",MatDim,0.,MatDim,MatDim,0.,MatDim);
	Flag2D=1; FlagTAC = 0;
	cout << "\nWhich type of the following matrix you want: \n";
		cout << "For [Gamma-Gamma(Enter \"1\")]  [DCO(Enter \"2\")]  [PDCO(Enter \"3\")] [Early_Delayed(Enter \"4\")] : ";
	cin >> MatCase;

if(MatCase != 3){
	his2D=new TH2I("his2D","his2D",MatDim,0.,MatDim,MatDim,0.,MatDim);
	//cout << "Enter file name for the matrix with .mat extension: ";
	//cin >> matfilnam;
	//MatFile = fopen(matfilnam, "wb");
	//MatSize = 2*(MatDim)*(MatDim);
	//Matrix = (unsigned short *)malloc(MatSize);
	//memset(Matrix,0,MatSize);
        }

if(MatCase == 4){
	cout << "Do you want Gamma-Time Matrix?" << endl;
	cin >> ans;
	if (!strncmp(ans, "Y", 1)||!strncmp(ans, "y", 1)) GammaTime = true;
	if(GammaTime == true){
	     his_GT=new TH2I("his_GT","his_GT",MatDim,0.,MatDim,MatDim,0.,MatDim);
        //cout << "Enter file name for the GT_matrix with .mat extension: ";
        //cin >> matfilnam1;
        //MatFile1 = fopen(matfilnam1, "wb");
        //MatSize = 2*(MatDim)*(MatDim);
        //Matrix_GT = (unsigned short *)malloc(MatSize);
        //memset(Matrix_GT,0,MatSize);
        }
	}


if(MatCase == 3){
	       his_par=new TH2I("his_par","his_par",MatDim,0.,MatDim,MatDim,0.,MatDim);
		   his_per=new TH2I("his_per","his_per",MatDim,0.,MatDim,MatDim,0.,MatDim);
        //cout << "Enter file name for Parallel vs. All matrix with .mat extension: ";
        //cin >> matfilnam_par;
        //cout << "Enter file name for Perpendicular vs. All matrix with .mat extension: ";
        //cin >> matfilnam_per;
        //MatFile1 = fopen(matfilnam_par, "wb");
        //MatFile2 = fopen(matfilnam_per, "wb");
        //MatSize = 2*(MatDim)*(MatDim);
        //Matrix_par = (unsigned short *)malloc(MatSize);
        //Matrix_per = (unsigned short *)malloc(MatSize);
        //memset(Matrix_par,0,MatSize);
        //memset(Matrix_per,0,MatSize);
        }

	cout << "\nDo you want to set a TAC gate? (Y/N): ";
	cin >> ans;
	if (!strncmp(ans, "Y", 1)||!strncmp(ans, "y", 1)){
	FlagTAC = 1;
	cout << "Enter the Lower and Upper Limit of the TAC: ";
	cin >> TACL >> TACH;
	cout << TACL << " " << TACH << endl;
		if(TACL > TACH){
		cout << "****ERROR****Lower Limit greater than the Upper Limit! **Exiting......**" << endl;
		exit(0);
		}
	}
}
/*.......Matrix file name and memory allocation stuff ends here......................*/

/****************** Dump file creation for CUBE sorting *****************/
cout << "\nDo you want to create a dump (.dmp) file for cube sorting? (Y/N): ";
cin >> ans;
if (!strncmp(ans, "Y", 1)||!strncmp(ans, "y", 1)){
	Flag3D = 1;
	cout << "Enter file name for the dump file with .dmp extension: ";
	cin >> dumpfilename;
	DumpFile = fopen(dumpfilename, "wb");
}
/**************************************************************************/
   do{
	if(first_file==false){
		first_event = true;
		printf("Enter the name of List Mode Data file: ");
		scanf ("%512s",filename);
	}

	if ((ifp = fopen (filename, "rb")) == (FILE*) NULL){
		fprintf (stderr, "File %s not found\n", filename);
		exit (1);
    	}



    while (! feof (ifp)) {
	if (fread (&readbuf.header, sizeof (dataheader), 1, ifp) != 1) {
	    fprintf (stderr, "Error reading header block\n");
	    exit (1);
	}
	if (readbuf.header.block == end_of_file) {
	    printf ("\n---END OF FILE---\n");
	    printf ("Total Number of Event Blocks = %lu\n", num_ev_blocks);
	    printf ("Num Events = BlockSize X EventSize = %d X %d = %lu\n\n",
		blockSize, eventSize, num_events);
	    break;
	}
	if ((readbuf.header.block < 0) || (readbuf.header.block > datarate)) {
	    struct stat	statBuf;

	    fstat (fileno (ifp), &statBuf);
	    fprintf (stderr, "\nUnknown header FILE-POS = %ld PREV = %d\n\n",
			ftell (ifp), prevBlock);
	    fprintf (stderr, "ByteSize = %ld BlockSize=%ld\n",
		statBuf.st_size, statBuf.st_blocks);
	    printHeader (&readbuf.header);
	    sleep (1);
//	    fseek (ifp, -sizeof (readbuf.header) + 1, SEEK_CUR);
	    continue;
	}
	if (fread (&readbuf.data, readbuf.header.size_in_bytes, 1, ifp) != 1) {
	    fprintf (stderr,
		"\n\n\nError reading data FILE-POS = %ld PREV = %d\n\n",
			ftell (ifp), prevBlock);
	    printHeader (&readbuf.header);
	    exit (1);
	}
	readbuf.size = sizeof (dataheader) + readbuf.header.size_in_bytes;
	expand (&readbuf, &buf);
	switch (buf.header.block) {
	    case datarate :
		memcpy (&rate, &readbuf, sizeof (ratepack));
		if (prevBlock != event) break;
		dcr = rate.evn_during_thisblock;
		if (rate.time_for_thisblock > 0) dcr /= rate.time_for_thisblock;
		icr = rate.evn_since_start;
		if (rate.time_since_start > 0) icr /= rate.time_since_start;
		/*printf ("%c[1;35m%8ld blocks in %.2lf sec \n",
			27, num_ev_blocks, rate.time_since_start*1e-3);
		printf (" DCR=%.2lf KE/s ICR=%.2lf KE/s%c[0;39m\r \n",
			dcr, icr, 27);*/
		fflush (stdout);
		break;

	    case event :

		if((first_event == true)){
		printHeader(&readbuf.header);
		NoPara = readbuf.header.unitsize;
		printf("No. of parameters are %d \n",NoPara);
		     if(first_file==true)
				    {   ifstream infile("Names.dat");
						for (unsigned int i=0; i<NoPara; ++i)
				        {  infile>>nm;
							/* int Cl_no,k1;
						 switch(i){
							 case 0:snprintf(nm, sizeof(nm), "Hit-1"); break;
							 case 1:snprintf(nm, sizeof(nm), "TLO-1");break;
							 case 2:snprintf(nm, sizeof(nm), "TMI-1");break;
							 case 3:snprintf(nm, sizeof(nm), "THI-1");break;
					         case 44:snprintf(nm, sizeof(nm), "Hit-2"); break;
							 case 45:snprintf(nm, sizeof(nm), "TLO-2");break;
							 case 46:snprintf(nm, sizeof(nm), "TMI-2");break;
							 case 47:snprintf(nm, sizeof(nm), "THI-2");break;
							 default:if(i>3 && i<44){
							 Cl_no=(i-4)/5+1;
							 k1=(i-4)%5;}
							 else{
							 Cl_no=(i-4)/5+1;k1=(i-3)%5;}
							 if(k1==0)
							 snprintf(nm, sizeof(nm), "CL%02i_T",Cl_no);
						     else snprintf(nm, sizeof(nm), "CL%02i_E%i",Cl_no,k1);
						           } */
                         his[i]=new TH1I(nm,nm,16384,0.,16384.);
						 strcat(nm,"_cal");
				         his_cal[i]=new TH1I(nm,nm,16384,0.,16384.);
						}
				        for (int i = 0; i < No_Clover; i++)
					    {
						 snprintf(nm, sizeof(nm), "TDCA_%02i", i);
				         his_TDC[i]=new TH1I(nm,nm,16384,0.,16384.);
						 snprintf(nm, sizeof(nm), "AB%02i", i);
				         his_AB[i]=new TH1I(nm,nm,16384,0.,16384.);
						 snprintf(nm, sizeof(nm), "AB_M1_%02i", i);
				         his_AB1[i]=new TH1I(nm,nm,16384,0.,16384.);
						 snprintf(nm, sizeof(nm), "AB_Mgt1_%02i", i);
				         his_AB2[i]=new TH1I(nm,nm,16384,0.,16384.);
				        						}
						 for (int i = 0; i < 8; i++)
						 {switch (i)
							 {case 0:snprintf(nm, sizeof(nm), "CL9_E_par");break;
							  case 1:snprintf(nm, sizeof(nm), "CL9_E_per");break;
							  case 2:snprintf(nm, sizeof(nm), "CL10_E_par");break;
							  case 3:snprintf(nm, sizeof(nm), "CL10_E_per");break;
							  case 4:snprintf(nm, sizeof(nm), "CL12_E_par");break;
							  case 5:snprintf(nm, sizeof(nm), "CL12_E_per");break;
							  case 6:snprintf(nm, sizeof(nm), "CL16_E_par");break;
							  case 7:snprintf(nm, sizeof(nm), "CL16_E_per");break;
							   }
						  his_AFC[i]=new TH1I(nm,nm,16384,0.,16384.);
						 }
						 his_CVT[0]=new TH1I("Time_Diff_i","Time_Diff_i",16384,0.,16384.);
						 his_CVT[1]=new TH1I("Time_Diff_p","Time_Diff_p",16384,0.,16384.);
						 his_CVT[2]=new TH1I("Time_Diff_bgi","Time_Diff_bgi",16384,0.,16384.);
						 his_CVT[3]=new TH1I("Time_Diff_bgp","Time_Diff_bgp",16384,0.,16384.);
						 his_CVT[4]=new TH1I("Time_Diff_ed","Time_Diff_ed",16384,0.,16384.);
						 his_CVT[5]=new TH1I("Time_Diff_bged","Time_Diff_bged",16384,0.,16384.);

                         his_TAC=new TH1I("TACspectra","TACspectra",16384,0.,16384.);

					}


		first_event=false;

		printf("Enter the name of file containing calibration parameters (e.g. calib.dat): ");
		cin  >> calfile;
		ifstream calfile2;
		calfile2.open(calfile.c_str());

		if(calfile2.fail())
			{
				cout << "ERROR*: Can't open calibration file:--> "<< calfile << endl;
				return(-2);
			}
		else{
			for(int l = 0; l < NoPara; l++){

				calfile2 >> coffset[l] >>  slope[l] >> quad[l];
				cout << l<<".\t" << coffset[l] << "\t" <<  slope[l] << "\t" << quad[l] << endl;	// Reading in calibration parameters
			}

		    }

		}
		++num_ev_blocks;
		num_events += process_event_block (&readbuf);



			for (unsigned int k=0; k<buf.header.number_of_units; ++k) {
				//cout << "buf.header.number_of_units  " << buf.header.number_of_units << endl;
			    int32 offset = k * buf.header.unitsize;

			    for (unsigned int i=0; i<buf.header.unitsize; ++i){ //Extracting RAW data and recalibrating to 0.5 kev/ch

				int data=0;
				energy[i] = 0.0;
		         	energy_cal[i] = 0.0;
				data = (int16)buf.data[offset + i];
				if(data > 0){
					energy[i] = data;		// Assigning data to each parameters - calling it as "energy"
					his[i]->Fill(energy[i]);
		             	   	//Putting RAW spectra in the his file (first "NoPara" spectra)
					/* Calibrating EVERY parameter to 0.5 keV/ch */
					energy_cal[i] = (unsigned short) (2*(coffset[i] + slope[i]*data + quad[i]*data*data)+((float) rand()/(float) RAND_MAX));
				     	his_cal[i]->Fill(energy_cal[i]);	// Putting calibrated spectra (next "NoPara" spectra)
			     		}
				}


				for(int l = 0; l < No_Clover; l++) CLT[l] = 0.0;
			/* Assigning & aligning TDC values to calculate TAC spectrum later*/
				if(energy[4] != 0.0)
					{CLT[0] = energy[4] - 372.0;}	// TDC of 1st Clover

				if(energy[9] != 0.0)
					{CLT[1] = energy[9] - 836.0;}	// ..

				if(energy[14] != 0.0)
					CLT[2] = energy[14] - 360.0;	// ..

				if(energy[19] != 0.0)
					CLT[3] = energy[19] - 316.0;	// ..

				if(energy[24] != 0.0)
					CLT[4] = energy[24] - 9.0;	// ..

				if(energy[29] != 0.0)
					CLT[5] = energy[29] - 372.0;	// ..

				if(energy[34] != 0.0)
					CLT[6] = energy[34] - 372.0;	// ..

				if(energy[39] != 0.0)
					CLT[7] = energy[39] - 299.0;	// ..

				if(energy[48] != 0.0)
					CLT[8] = energy[48] - 439.0;	// ..

				if(energy[53] != 0.0)
					CLT[9] = energy[53] - 432.0;	// ..

				if(energy[58] != 0.0)
					CLT[10] = energy[58] - 0.0;	// ..

				if(energy[63] != 0.0)
					CLT[11] = energy[63] - 381.0;	// ..

				if(energy[68] != 0.0)
					CLT[12] = energy[68] - 125.0;	// ..

				if(energy[73] != 0.0)
					CLT[13] = energy[73] - 451.0;	// ..

				if(energy[78] != 0.0)
					CLT[14] = energy[78] - 0.0;	// ..

				if(energy[83] != 0.0)
					CLT[15] = energy[83] - 0.0;	// TDC of 16th Clover

				for (int n = 0; n < No_Clover; n++){
					if (CLT[n] != 0.0) his_TDC[n]->Fill(CLT[n]);	// Putting aligned TDC spectra(Next 16 spectra)
				}
			/* Performing Addback here */
			AddBack();
				for (int n = 0; n < No_Clover; n++){
					if((Eclab[n] > 50) && (Eclab[n] < 8192) ){
					his_AB[n]->Fill(Eclab[n]);}	// Putting AddBack spectra (Next 16 spectra)
				}
			/* Addback done!*/
                       if(clov_mult == 1){
                       for (int n = 0; n < No_Clover; n++){
                                        if((Eclab[n] > 50) && (Eclab[n] < 8192) ){
                                        his_AB1[n]->Fill(Eclab[n]);}
                       }
                       }

                       if(clov_mult > 1){
                       for (int n = 0; n < No_Clover; n++){
                                        if((Eclab[n] > 50) && (Eclab[n] < 8192) ){
                                         his_AB2[n]->Fill(Eclab[n]);}
                       }
                       }

		       // copy current event data to all tree leaves
		       // (warning: if for tree leaves one uses more compact
		       // variables' types than for the event data, one should
		       // make sure here that the current values do not exceed
		       // limits of these more compact variables's types)
		       _clov_mult_ = clov_mult;
		       _No_Clover_ = 0;
		       for (int n = 0; n < No_Clover; n++) {
#if (_ZERO_SUPPRESSION_LEVEL_) > 0
			 // the "if" for the "zero-suppressed" tree case only
			 // (you can use any conditions which you want here)
			 if ((CLT[n] == 0) && (Eclab[n] == 0)) continue;
#endif /* (_ZERO_SUPPRESSION_LEVEL_) > 0 */
			 _Clover_[(_No_Clover_)] = n + 1;
			 _CLT_[(_No_Clover_)] = CLT[n];
			 _Eclab_[(_No_Clover_)] = Eclab[n];
			 _No_Clover_++;
		       }

		       // fill the tree with the current event data
		       // (you can add any conditions which you want here)
		       if (t
#if (_ZERO_SUPPRESSION_LEVEL_) > 1
			   && (_clov_mult_ != 0)
			   && (_No_Clover_ != 0)
#endif /* (_ZERO_SUPPRESSION_LEVEL_) > 1 */
			   ) t->Fill();

//Asymmetry factor Calculation only for 90 degree detectors and source data//
for(int i = 8; i <16; i++){
E_par1[i] = 0.0;
E_par2[i] = 0.0;
E_par[i] = 0.0;
E_per1[i] = 0.0;
E_per2[i] = 0.0;
E_per[i] = 0.0;
}
///*****Clover-9*****/
if ((energy_cal[49] > 10) && (energy_cal[52] > 10))
{
E_par1[8] = energy_cal[49] + energy_cal[52];
}
if ((energy_cal[50] > 10) && (energy_cal[51] > 10))
{
E_par2[8] = energy_cal[50] + energy_cal[51];
}
if ((energy_cal[49] > 10) && (energy_cal[50] > 10))
{
E_per1[8] = energy_cal[49] + energy_cal[50];
}
if ((energy_cal[51] > 10) && (energy_cal[52] > 10))
{
E_per2[8] = energy_cal[51] + energy_cal[52];
}
E_par[8] =  E_par1[8] + E_par2[8];
if (E_par[8] > 10) {his_AFC[0]->Fill(E_par[8]);}
E_per[8] =  E_per1[8] + E_per2[8];
if (E_per[8] >10) {his_AFC[1]->Fill(E_per[8]);}
///****************/

///*****Clover-10*****/
if ((energy_cal[54] > 10) && (energy_cal[57] > 10))
{
E_par1[9] = energy_cal[54] + energy_cal[57];
}
if ((energy_cal[55] > 10) && (energy_cal[56] > 10))
{
E_par2[9] = energy_cal[55] + energy_cal[56];
}
if ((energy_cal[54] > 10) && (energy_cal[55] > 10))
{
E_per1[9] = energy_cal[54] + energy_cal[55];
}
if ((energy_cal[56] > 10) && (energy_cal[57] > 10))
{
E_per2[9] = energy_cal[56] + energy_cal[57];
}
E_par[9] =  E_par1[9] + E_par2[9];
if (E_par[9] > 10) {his_AFC[2]->Fill(E_par[9]);}
E_per[9] =  E_per1[9] + E_per2[9];
if (E_per[9] > 10) {his_AFC[3]->Fill(E_per[9]);}
///*******************/

///*****Clover-12*****/         // 11th detector is not present in array
if ((energy_cal[64] > 10) && (energy_cal[67] > 10))
{
E_par1[11] = energy_cal[64] + energy_cal[67];
}
if ((energy_cal[65] > 10) && (energy_cal[66] > 10))
{
E_par2[11] = energy_cal[65] + energy_cal[66];
}
if ((energy_cal[64] > 10) && (energy_cal[65] > 10))
{
E_per1[11] = energy_cal[64] + energy_cal[65];
}
if ((energy_cal[66] > 10) && (energy_cal[67] > 10))
{
E_per2[11] = energy_cal[66] + energy_cal[67];
}
E_par[11] =  E_par1[11] + E_par2[11];
if (E_par[11] > 10) {his_AFC[4]->Fill(E_per[11]);}
E_per[11] =  E_per1[11] + E_per2[11];
if (E_per[11] > 10) {his_AFC[5]->Fill(E_per[11]);}
///***************/

///*****Clover-16*****/        // Not considering 13th and 14th detector because of bad resolution                             // 15th detector is not present
if ((energy_cal[84] > 10) && (energy_cal[85] > 10))
{
E_par1[15] = energy_cal[84] + energy_cal[85];
}
if ((energy_cal[86] > 10) && (energy_cal[87] > 10))
{
E_par2[15] = energy_cal[86] + energy_cal[87];
}
if ((energy_cal[85] > 10) && (energy_cal[86] > 10))
{
E_per1[15] = energy_cal[85] + energy_cal[86];
}
if ((energy_cal[84] > 10) && (energy_cal[87] > 10))
{
E_per2[15] = energy_cal[84] + energy_cal[87];
}
E_par[15] =  E_par1[15] + E_par2[15];
if (E_par[15] > 10) {his_AFC[6]->Fill(E_per[15]);}
E_per[15] =  E_per1[15] + E_per2[15];
if (E_per[15] > 10) {his_AFC[7]->Fill(E_per[15]);}
/*******************/

//Count_versus_time//
			CountVsTime();
//Gamma-Gamma Matrix//
			if((Flag2D == 1) && (MatCase == 1)) {MatGG();}
//DCO Matrix//
			if((Flag2D == 1) && (MatCase == 2)) {MatDCO();}
//PDCO Matrix//
			if((Flag2D == 1) && (MatCase == 3)) {MatPDCO();}
//ThreeD//
			if(Flag3D == 1) {ThreeD();}
//Early-Delayed Matrix//
			if((Flag2D == 1) && (MatCase == 4)) {MatED();}
//
		}



		if (! blockSize) {
		    blockSize = num_events;
		    eventSize = readbuf.header.unitsize;
		}

	        break;

	    case hgram :
		if (prevBlock != hgram) {
		    printf ("\n");
		    singID = 0;
		}
		printf ("\nMCA-%02d : Size = %d (TimeSinceRunStart=%d mins)\n",
			++singID, readbuf.header.unitsize,
			(int32) (rate.time_since_start/ONE_MINUTE));
	        break;

	    case scaler :
		printf ("\n--- SCALER BLOCK BEGIN ---\n");
		printf ("%s\n", (char *) readbuf.data);
		printf ("--- SCALER BLOCK END ---\n");
		break;

	    case start :
		printf ("\n--- START BLOCK BEGIN ---\n");
		printf ("%s", (char *) readbuf.data);
		printf ("--- START BLOCK END ---\n\n");
		break;

	    case stop :
		printf ("\n\n--- STOP BLOCK BEGIN ---\n");
		printf ("%s", (char *) readbuf.data);
		printf ("--- STOP BLOCK END ---\n\n");
		break;

	    case Pause :
		printf ("\n--- PAUSE BLOCK BEGIN ---\n");
		printf ("%s\n", (char *) readbuf.data);
		printf ("--- PAUSE BLOCK END ---\n");
		break;

	    case resume :
		printf ("\n--- RESUME BLOCK BEGIN START ---\n");
		printf ("%s\n",(char *) readbuf.data);
		printf ("--- RESUME BLOCK END ---\n");
		break;

	    case names :
		if(first_file==true){
		printf ("\n--- NAMES BLOCK BEGIN ---\n");
		printf ("%s", (char *) readbuf.data);
		printf ("--- NAMES BLOCK END ---\n");
		}
		break;

	    case user :
		if (! strncasecmp ((char*)readbuf.data, "$$FREE$$", 8)) {
		    printf ("\n--- USER BLOCK BEGIN ---\n");
		    printf ("%s", ((char *) readbuf.data) + 8);
		    printf ("--- USER BLOCK END ---\n");
		} else {
		    fprintf (stderr, "\n\n***Not a 'candle' format file\n\n");
		    exit (1);
		}
		break;

	    case error :
		printf ("\n--- ERROR BLOCK BEGIN ---\n");
		printf ("%s\n", (char *) readbuf.data);
		printf ("--- ERROR BLOCK END ---\n");
		break;

	    default :
		printf ("Unknown Blocktype %d\n", readbuf.header.block);
	}
	prevBlock = buf.header.block;
    }

    if (t) t->Write(); // "flush" tree entries from the current data file

			cout << "Zero Fold Data Events: "  << ZeroFold  << endl;
			cout << "One Fold Data Events: "   << OneFold   << endl;
			cout << "Two Fold Data Events: "   << TwoFold   << endl;
			cout << "Three Fold Data Events: " << ThreeFold << endl;
			cout << "Four Fold Data Events: "  << FourFold  << endl;
			cout << "Five Fold Data Events: "  << FiveFold  << endl;
			cout << "Six Fold Data Events: "   << SixFold   << endl;
			cout << "Sixteen Fold Data Events: "  << SixteenFold  << endl;

/* If the file is successfully sorted then ask if more data files are to be sorted */
   printf("\nDo you want to sort another file (Y/N)?");
   scanf("%s",ans);
   if (!strncmp(ans, "Y", 1)||!strncmp(ans, "y", 1)){
   	flag=1;
   	first_file=false;}
   else
   	flag=0;
  }while(flag==1);



/* if((Flag2D ==1)&& (MatCase != 3)){
    rembckgnd(his2D);
	//fwrite(Matrix, 1, MatSize, MatFile);
	//fclose(MatFile);
	//free(Matrix);
}*/
/*
if((Flag2D ==1)&& (MatCase == 4) && (GammaTime == true)){
	his_GT->Write();
	//fwrite(Matrix_GT, 1, MatSize, MatFile1);
	//fclose(MatFile1);
	//free(Matrix_GT);
}
if((Flag2D ==1)&& (MatCase == 3)){
	his_par->Write();
	//fwrite(Matrix_par, 1, MatSize, MatFile1);
	//fclose(MatFile1);
	//free(Matrix_par);
	his_per->Write();
	//fwrite(Matrix_per, 1, MatSize, MatFile2);
	//fclose(MatFile2);
	//free(Matrix_per);
}*/
f->cd();
f->Write();
delete f; // automatically deletes the tree and all histograms
if(Flag3D==1){
fclose(DumpFile);
			cout << "\n Valid 3 fold events: " << valid_3F_events << endl;
}
time(&EndTimer);
seconds = difftime(EndTimer,StartTimer);
cout << "Elapsed time " << seconds << " sec." << endl;

    return 0;
}


void expand (datapackptr inbuf, datapackptr outbuf)
{
int32		nwords = inbuf->header.number_of_units * inbuf->header.unitsize;
int16ptr 	wrtptr,rdptr,patptr;
int32		i, numpats;
unsigned short mask;

    if (inbuf->header.block != event) {
	memcpy (outbuf, inbuf, inbuf->size);
	return;
    }

    memset (outbuf, 0, sizeof (datapack));

    numpats = nwords / 16;
    if (nwords % 16) ++numpats;

    rdptr = inbuf->data + numpats;
    patptr = inbuf->data;
    wrtptr = outbuf->data;
    mask = 1;

    for (i=0; i<nwords; ++i) {
	if (mask & *patptr)
	    *wrtptr++ = *rdptr++;
	else
	    *wrtptr++ = 0;
	mask <<= 1;
	if (! mask) {
	    mask = 1;
	    ++patptr;
	}
    }

    outbuf->header.block  = inbuf->header.block;
    outbuf->header.number_of_units = inbuf->header.number_of_units;
    outbuf->header.unitsize = inbuf->header.unitsize;
    outbuf->header.compstatus = 0;
    outbuf->header.size_in_bytes = 2 * inbuf->header.number_of_units
		* inbuf->header.unitsize;

}
/*-------Performing Add-Back below-----------*/
void AddBack(){

cl_mult.reset();

for(int i = 0; i < No_Clover; i++){
Eclab[i]=0.0;}

/*****CLOVER-1******/
for(int k = 5; k < 9; k++){
Eclab[0]+=energy_cal[k];
Eclab[0] = Eclab[0] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[0] > 10) cl_mult.set(0);	//Eclab=> Energy Clover Addback
/******************/
/*****CLOVER-2******/
for(int k = 10; k < 14; k++){
Eclab[1]+=energy_cal[k];
Eclab[1] = Eclab[1] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[1] > 10) cl_mult.set(1);
/******************/
/*****CLOVER-3******/
for(int k = 15; k < 19; k++){
Eclab[2]+=energy_cal[k];
Eclab[2] = Eclab[2] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[2] > 10) cl_mult.set(2);
/******************/
/*****CLOVER-4******/
for(int k = 20; k < 24; k++){
Eclab[3]+=energy_cal[k];
Eclab[3] = Eclab[3] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[3] > 10) cl_mult.set(3);
/******************/
/*****CLOVER-5******/
for(int k = 25; k < 29; k++){
Eclab[4]+=energy_cal[k];
Eclab[4] = Eclab[4] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[4] > 10) cl_mult.set(4);
/******************/
/*****CLOVER-6******/
for(int k = 30; k < 34; k++){
Eclab[5]+=energy_cal[k];
Eclab[5] = Eclab[5] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[5] > 10) cl_mult.set(5);
/******************/
/*****CLOVER-7******/
for(int k = 35; k < 39; k++){
Eclab[6]+=energy_cal[k];
Eclab[6] = Eclab[6] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[6] > 10) cl_mult.set(6);
/******************/
/*****CLOVER-8******/
for(int k = 40; k < 44; k++){
Eclab[7]+=energy_cal[k];
Eclab[7] = Eclab[7] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[7] > 10) cl_mult.set(7);
/******************/
/*****CLOVER-9******/
for(int k = 49; k < 53; k++){
Eclab[8]+=energy_cal[k];
Eclab[8] = Eclab[8] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[8] > 10) cl_mult.set(8);
/******************/
/*****CLOVER-10*****/
for(int k = 54; k < 58; k++){
Eclab[9]+=energy_cal[k];
Eclab[9] = Eclab[9] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[9] > 10) cl_mult.set(9);
/******************/
/*****CLOVER-11*****/
for(int k = 59; k < 63; k++){
Eclab[10]+=energy_cal[k];
Eclab[10] = Eclab[10] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[10] >50) cl_mult.set(10);
/******************/
/*****CLOVER-12*****/
for(int k = 64; k < 68; k++){
Eclab[11]+=energy_cal[k];
Eclab[11] = Eclab[11] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[11] > 10) cl_mult.set(11);
/******************/
/*****CLOVER-13*****/
for(int k = 69; k < 73; k++){
Eclab[12]+=energy_cal[k];
Eclab[12] = Eclab[12] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[12] > 10) cl_mult.set(12);
/******************/
/*****CLOVER-14*****/
for(int k = 74; k < 78; k++){
Eclab[13]+=energy_cal[k];
Eclab[13] = Eclab[13] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[13] > 10) cl_mult.set(13);
/******************/
/*****CLOVER-15*****/
for(int k = 79; k < 83; k++){
Eclab[14]+=energy_cal[k];
Eclab[14] = Eclab[14] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[14] > 10) cl_mult.set(14);
/******************/
/*****CLOVER-16*****/
for(int k = 84; k < 88; k++){
Eclab[15]+=energy_cal[k];
Eclab[15] = Eclab[15] + (float) rand()/(float) RAND_MAX;
}
if(Eclab[15] > 10) cl_mult.set(15);
/******************/

clov_mult = cl_mult.count();
if(clov_mult == 0) ZeroFold++;
if(clov_mult == 1) OneFold++;
if(clov_mult == 2) TwoFold++;
if(clov_mult == 3) ThreeFold++;
if(clov_mult == 4) FourFold++;
if(clov_mult == 5) FiveFold++;
if(clov_mult == 6) SixFold++;
if(clov_mult == 16) SixteenFold++;
}
/**************************/
/**************************/
void CountVsTime(){

unsigned long Time_Differ_i, Time_Differ_p, Time_Differ_bgi, Time_Differ_bgp, Time_Differ_ed, Time_Differ_bged;

if(clov_mult > 1){      //Checking if atleast 2 Clovers fire
        for(int j = 0; j < 16; j++){    //Looping over the clovers
        for(int k = 0; k < 16; k++){    //Looping over the clovers
            if ((k != j) && (CLT[j] != 0.0) && (CLT[k] != 0.0)){
		if  ((((Eclab[j] >= 1088) && (Eclab[j] <= 1106)) || ((Eclab[j] >= 1050) && (Eclab[j] <= 1065))) && (Eclab[k] >= 918 && Eclab[k] <= 930)){
            	Time_Differ_i = 8000.0 + (CLT[k] - CLT[j]);
            	his_CVT[0]->Fill(Time_Differ_i);
	    	}
                if  ((((Eclab[k] >= 1088) && (Eclab[k] <= 1106)) || ((Eclab[k] >= 1050) && (Eclab[k] <= 1065))) && (Eclab[j] >= 918 && Eclab[j] <= 930)){
                Time_Differ_i = 8000.0 + (CLT[j] - CLT[k]);
                his_CVT[0]->Fill(Time_Differ_i);
                }
                if ((((Eclab[j] >= 821) && (Eclab[j] <= 834)) || ((Eclab[j] >= 876) && (Eclab[j] <= 890))) && (Eclab[k] >= 1048 && Eclab[k] <= 1062)){
                Time_Differ_p = 8000.0 + (CLT[k] - CLT[j]);
                his_CVT[1]->Fill(Time_Differ_p);
                }
                if ((((Eclab[j] >= 1070) && (Eclab[j] <= 1086)) || ((Eclab[j] >= 1032) && (Eclab[j] <= 1046))) && (Eclab[k] >= 931 && Eclab[k] <= 943)){
                Time_Differ_bgi = 8000.0 + (CLT[k] - CLT[j]);
                his_CVT[2]->Fill(Time_Differ_bgi);
                }
                if ((((Eclab[j] >= 835) && (Eclab[j] <= 848)) || ((Eclab[j] >= 860) && (Eclab[j] <= 874))) && (Eclab[k] >= 1032 && Eclab[k] <= 1046)){
                Time_Differ_bgp = 8000.0 + (CLT[k] - CLT[j]);
                his_CVT[3]->Fill(Time_Differ_bgp);
                }
                if ((((Eclab[j] >= 1088) && (Eclab[j] <= 1106)) || ((Eclab[j] >= 495) && (Eclab[j] <= 507)) || ((Eclab[j] >= 314) && (Eclab[j] <= 326))) && (Eclab[k] >= 260 && Eclab[k] <= 273)){
                Time_Differ_ed = 8000.0 + (CLT[k] - CLT[j]);
                his_CVT[4]->Fill(Time_Differ_ed);
                }
                if ((((Eclab[j] >= 1070) && (Eclab[j] <= 1086)) || ((Eclab[j] >= 508) && (Eclab[j] <= 520)) || ((Eclab[j] >= 300) && (Eclab[j] <= 312))) && (Eclab[k] >= 274 && Eclab[k] <= 287)){
                Time_Differ_bged = 8000.0 + (CLT[k] - CLT[j]);
                his_CVT[5]->Fill(Time_Differ_bged);
                }
              }
            }
          }
        }
     }
/*------------------*/
/*------------------*/
void MatGG(){
if(clov_mult > 1){	//Checking if atleast 2 Clovers fire
	for(int j = 0; j < 16; j++){	//Looping over the clovers	!X-Axis
	for(int k = 0; k < 16; k++){	//Looping over the clovers	!Y-Axis
		if((k != j) && (CLT[j] != 0.0) && (CLT[k] != 0.0)){	// Asymmetric Matrix AND Checking corresponding TDC signals
			Time_Diff = 8192.0 + (CLT[j] - CLT[k]);		// Calculating TAC spectra

		if(FlagTAC == 0){	// Gamma-Gamma Matrix WITHOUT a TAC condition
			his_TAC->Fill(Time_Diff);		// Putting TAC spectra
			x = 0,y = 0,z = 0;
			if ((Eclab[j] > 20) && (Eclab[j] < MatDim) && (Eclab[k] > 20) && (Eclab[k] < MatDim)){
			x = Eclab[j]; y = Eclab[k];
			his2D->Fill(x,y);
			//z = x+MatDim*y; Matrix[z]++;
			}
		}

		if(FlagTAC == 1){	// Gamma-Gamma Matrix WITH a TAC condition
			his_TAC->Fill(Time_Diff);		// Putting TAC spectra
			x = 0,y = 0,z = 0;
			if ((Time_Diff >= TACL) && (Time_Diff <= TACH) && ((Eclab[j] > 20) && (Eclab[j] < MatDim) && (Eclab[k] > 20) && (Eclab[k] < MatDim))){
			x = Eclab[j]; y = Eclab[k];
			his2D->Fill(x,y);
			//z = x+MatDim*y; Matrix[z]++;
			}
		}
		}
	}
	}
}
}
/**************************/
/**************************/
void MatDCO(){
//DCO Matrix
if(clov_mult > 1){
	for(int j = 0; j < 7; j++){	//Looping over detectors at angles other than 90 deg.	!X-Axis
	for(int k = 8; k < 16; k++){ //Looping over detectors at 90 deg.		!Y-Axis
		if ((CLT[j] != 0.0) && (CLT[k] != 0.0)){
			Time_Diff = 8192.0 + (CLT[j] - CLT[k]);}         // Calculating TAC spectra

		if(FlagTAC == 0){	// DCO Matrix WITHOUT a TAC condition
	  		his_TAC->Fill(Time_Diff);		// Putting TAC spectra
			x = 0,y = 0,z = 0;
			if ((Eclab[j] > 20) && (Eclab[j] < MatDim) && (Eclab[k] > 20) && (Eclab[k] < MatDim)){
			x = Eclab[j]; y = Eclab[k];
		    his2D->Fill(x,y);
			//z = x+MatDim*y; Matrix[z]++;

			}
		}

		if(FlagTAC == 1){	// DCO Matrix WITH a TAC condition
		his_TAC->Fill(Time_Diff);		// Putting TAC spectra
			x = 0,y = 0,z = 0;
			if ((Time_Diff >= TACL) && (Time_Diff <= TACH) && ((Eclab[j] > 20) && (Eclab[j] < MatDim) && (Eclab[k] > 20) && (Eclab[k] < MatDim))){
			x = Eclab[j]; y = Eclab[k];
			his2D->Fill(x,y);
			//z = x+MatDim*y; Matrix[z]++;
			}
		}
		}
	}
	}

}
/**************************/
/**************************/
void MatPDCO(){
//PDCO Matrix

for(int i = 0; i<6; i++){
EM_par[i] = 0.0;
EM_per[i] = 0.0;
}

EM_par[0] = E_par1[8];
EM_par[1] = E_par2[8];
EM_par[2] = E_par1[9];
EM_par[3] = E_par2[9];
EM_par[4] = E_par1[15];
EM_par[5] = E_par2[15];
EM_per[0] = E_per1[8];
EM_per[1] = E_per2[8];
EM_per[2] = E_per1[9];
EM_per[3] = E_per2[9];
EM_per[4] = E_per1[15];
EM_per[5] = E_per2[15];
CLT_M[0] = CLT[8];
CLT_M[1] = CLT[8];
CLT_M[2] = CLT[9];
CLT_M[3] = CLT[9];
CLT_M[4] = CLT[15];
CLT_M[5] = CLT[15];

if(clov_mult>1){
for(int j = 0; j<8; j++){//looping over detectors at angles other than 90 degree
for(int k = 0; k<6; k++){//looping over detectors only at 9th, 10th and 16th position of array
         if((CLT[j] != 0.0) && CLT_M[k]!= 0.0){
         Time_Diff = 8192.0 + (CLT[j] - CLT_M[k]);} // Calculating TAC Spectra

if (FlagTAC == 0){//PDCO matrix without TAC condition
         his_TAC->Fill(Time_Diff); // Writing TAC Spectra
         x = 0, y = 0, z = 0;
         if((Eclab[j] > 20) && (Eclab[j] < MatDim) && (EM_par[k] > 20) && (EM_par[k] < MatDim)){
         x = Eclab[j];
         y = EM_par[k];
		 his_par->Fill(x,y);
         //z = x+MatDim*y;
         //Matrix_par[z]++;
		 }

         x = 0, y = 0, z = 0;
         if((Eclab[j] > 20) && (Eclab[j] < MatDim) && (EM_per[k] > 20) && (EM_per[k] < MatDim)){
         x = Eclab[j];
         y = EM_per[k];
		 his_per->Fill(x,y);
         //z = x+MatDim*y;
         //Matrix_per[z]++;
		 }
         }

if (FlagTAC == 1){//PDCO matrix with TAC condition
         his_TAC->Fill(Time_Diff); // Writing TAC Spectra
         x = 0, y = 0, z = 0;
         if((Time_Diff >= TACL) && (Time_Diff <=TACH) && (Eclab[j] > 20) && (Eclab[j] < MatDim) && (EM_par[k] > 20) && (EM_par[k] < MatDim)){
         x = Eclab[j];
         y = EM_par[k];
		 his_par->Fill(x,y);
         //z = x+MatDim*y;
         //Matrix_par[z]++;
		 }

         x = 0, y = 0, z = 0;
         if((Time_Diff >= TACL) && (Time_Diff <=TACH) && (Eclab[j] > 20) && (Eclab[j] < MatDim) && (EM_per[k] > 20) && (EM_per[k] < MatDim)){
         x = Eclab[j];
         y = EM_per[k];
		 his_per->Fill(x,y);
         //z = x+MatDim*y;
         //Matrix_per[z]++;
		 }

}
}
}
}
}
/**************************/
/**************************/
/* Cube stuff begins here*/
void ThreeD(){

bool result;
unsigned long diff_CubT;
//array<unsigned long,6> CubT_Diff = {999999999,999999999,999999999,999999999,999999999,999999999};
vector<unsigned long> ClovT_Diff;

	for (int j=0; j<4; j++ ) CubBuf[j]=0;

	if(clov_mult > 2){
		int j = 0;
		int l = CubBufEvent;

		for (int i=0; i<16; i++){
       			if (Eclab[i]>50.){
					CubBuf[j]=Eclab[i];
					CubT[j] = CLT[i];
					j++;}
    		}

		int m = 0;
		for(int k = 0; k < j; k++){
			for(int l = 0; l < j; l++){
				diff_CubT = 999999999;
				if (l>k){
					diff_CubT = abs(CubT[k] - CubT[l]);
//					CubT_Diff[m] = diff_CubT;
					ClovT_Diff.push_back(diff_CubT);
					m++;}
			}
		}


//		result = all_of(CubT_Diff.begin(), CubT_Diff.end(), [](int i){return ((i<1001) || (i==999999999));});
		result = all_of(ClovT_Diff.begin(), ClovT_Diff.end(), [](int j){return (j<1000);}); //Here the number "1000" defines the prompt window
		ClovT_Diff.clear();

               	if ((j>2) && (j<5) && (result==true)){	//With time condition defined above. Also, ignoring events with folds >
//               	if ((j>2) && (j<5)){	//Ignoring events with folds > 4

               		for (int n=0; n<4; n++)
				{sbuf[4*l+4+n]=CubBuf[n];}
               			CubBufEvent++;
				valid_3F_events++;
		}

		if(CubBufEvent > 8191) {
			//cout << CubBufEvent << endl;
//======================Cube Buffer Update=================================
			BufSize = 8*CubBufEvent;
			cout << "Buffer # " << BufNum << endl;
			//write(TDF,&BufSize,4);
			fwrite(&BufSize, 1, 4, DumpFile);
			sbuf[0]=-8;
			//cout << "sbuf[0] is: " << sbuf[0] << endl;
			sbuf[1] = BufSize/2;
			//cout << "sbuf[1] is: " << sbuf[1] << endl;
			sbuf[2]=4;
			//cout << "sbuf[2] is: " << sbuf[2] << endl;
			sbuf[3]=4*CubBufEvent;
			//cout << "sbuf[3] is: " << sbuf[3] << endl;
			//for (int i=0;i<CubBufEvent;i++)
			//{
			//printf("Blocksize : %d %d %d %d\n",sbuf[4*i+4],sbuf[4*i+5],sbuf[4*i+6],sbuf[4*i+7]);
			//}
			//write(TDF,sbuf,BufSize);
			fwrite(sbuf, 1, BufSize, DumpFile);
			//
			//write(TDF,&BufSize,4);
			fwrite(&BufSize, 1, 4, DumpFile);
			//break;
			CubBufEvent = 0;
			BufNum++;
		}
//=======================End of Cube Buffer Update=========================
	}
}

/* Cube stuff ends here*/

/**************************/
/**************************/
/**Early Delayed Matrix
 *
 * NOT TESTED!
 *
 * **/
void MatED(){

unsigned long E_Early, E_Delayed, T_Diff, Delta_T;

if(clov_mult > 2){      //Checking if atleast 2 Clovers fire
        for(int j = 0; j < 16; j++){    //Looping over the clovers
        for(int k = 0; k < 16; k++){    //Looping over the clovers
                if((k > j) && (CLT[j] != 0.0) && (CLT[k] != 0.0)){      // Asymmetric Matrix AND Checking corresponding TDC signals
                        Time_Diff = 8192.0 + (CLT[j] - CLT[k]);         // Calculating TAC spectra
                        his_TAC->Fill(Time_Diff);          // Putting TAC spectra

			T_Diff = std::abs(static_cast<int>(CLT[j] - CLT[k]));
			if((T_Diff >= 1008) && (T_Diff <= 2008)){
				if(CLT[j]>CLT[k]){
					E_Early = Eclab[k];
					E_Delayed = Eclab[j];
                                        Delta_T = CLT[j] - CLT[k];
                                        }
					else if (CLT[j]<CLT[k]){
					E_Early = Eclab[j];
					E_Delayed = Eclab[k];
                                        Delta_T = CLT[k] - CLT[j];
					}
                        x = 0,y = 0,z = 0;
                        if ((Eclab[j] > 20) && (Eclab[j] < MatDim) && (Eclab[k] > 20) && (Eclab[k] < MatDim)){
                        x = E_Early; y = E_Delayed; his2D->Fill(x,y);//z = x+MatDim*y; Matrix[z]++;
						}
			}

  			if(GammaTime == true){
				if(CLT[j]>CLT[k]){
					E_Early = Eclab[k];
					E_Delayed = Eclab[j];
                                        Delta_T = 1000+(CLT[j] - CLT[k]);
                                        }
					else if (CLT[j]<CLT[k]){
					E_Early = Eclab[j];
					E_Delayed = Eclab[k];
                                        Delta_T = 1000+(CLT[k] - CLT[j]);
					}
                        x = 0,y = 0,z = 0;
                        //if (((E_Delayed > 100) && (E_Delayed < MatDim)) && (((E_Early > 706) && (E_Early < 732)) || ((E_Early > 747) && (E_Early < 757)))){ // for 213At_GDT
                        // if (((E_Delayed > 100) && (E_Delayed < MatDim)) && (((E_Early > 876) && (E_Early < 890)) || ((E_Early > 821) && (E_Early < 833)))){  // for Prompt_216Fr_GDT
                        if (((E_Delayed > 100) && (E_Delayed < MatDim)) && (((E_Early > 317) && (E_Early < 325)) || ((E_Early > 1092) && (E_Early < 1104)) || ((E_Early > 496) && (E_Early < 507)))){  // for 133_216Fr_GDT
                        //if (((E_Delayed > 100) && (E_Delayed < MatDim))){
                        x = Delta_T; y = E_Delayed; his_GT->Fill(x,y); //z = x+MatDim*y; Matrix_GT[z]++;
						}
			}

                }
        }
        }
}
}

/**Early Delayed Matrix Ends**/
/*
void rembckgnd(TH1* back)
{
 Int_t i, j;
   Int_t nbinsx =(Int_t) back->GetNbinsX();
   Int_t nbinsy =(Int_t) back->GetNbinsY();
   Float_t ** source = new Float_t*[nbinsx];
   for (i=0;i<nbinsx;i++)
      source[i]=new Float_t[nbinsy];
   TSpectrum2 *s = new TSpectrum2();
   for (i = 0; i < nbinsx; i++){
      for (j = 0; j < nbinsy; j++){
         source[i][j] = back->GetBinContent(i + 1,j + 1);
      }
   }
   s->Background(source,nbinsx,nbinsy,4,4,TSpectrum2::kBackIncreasingWindow,TSpectrum2::kBackOneStepFiltering);
   for (i = 0; i < nbinsx; i++){
      for (j = 0; j < nbinsy; j++)
         bckgnd->SetBinContent(i + 1,j + 1,(Int_t)source[i][j]);
   }
   his2D_wb->Add(back,bckgnd,1,-1);

}
*/
