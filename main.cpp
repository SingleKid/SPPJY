#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "JTime.h"
#include "MatC.h"

// USER LEVEL
#define INPUT_OBS "./RAW/Thu Jan 10 13_00_26 2019_.obs"
//#define INPUT_OBS "./RAW/UCAL00CAN_S_20190100000_01D_30S_MO.rnx"
#define INPUT_GPS_BRD "./RAW/brdc0100s.19n"
#define INPUT_GLO_BRD "./RAW/brdc0100.19g"
#define INPUT_INX "./RAW/emrg0100.19i"
#define INPUT_SP3 "./RAW/igv20354_18.sp3"
#define ELEV_THRESH 0.20
#define SAT_NUM_THRESH 4
#define LS_CONV_THRESH 1e-6
#define LS_MAX_ITER 30
#define OBS_SIG0 1

// DATA LEVEL
#define MAX_SAT_NUM 200
#define GLO_OFFSET 100

// RINEX LEVEL
#define GPS_FLAG 'G'
#define GLO_FLAG 'R'
#define MAX_LINE 400
#define END_OF_HEADER "END OF HEADER"

/// IONEX
#define CLK_AR "AR"
#define CLK_AS "AS"
#define INX_LAT_BINS 71
#define INX_LON_BINS 73

/// SP3
#define SP3_HEADER_LINE 22
#define SP3_EOF "EOF"

// CONSTS
#define LIGHT_SPEED 299792458.0
#define we 7.2921151467e-5
#define mu 3.986004415E14
#define R1 4.442807633e-10

const double station_blh[3] = { 0.891513720275876, -1.99201147820819, 1118.80358502921 }; // ucal
const double station_xyz[3] = { -1641945.2000, -3664804.1000, 4940009.3000 };             // ucal
const double center_earth[3] = { 0,0,0 };

const int ban_svn_amount = 2;
const int ban_svn[ban_svn_amount]{4, 27};


double s_sinB = sin(station_blh[0]);
double s_cosB = cos(station_blh[0]);
double s_sinL = sin(station_blh[1]);
double s_cosL = cos(station_blh[1]);

struct OBS_FRAME {
	//double C1C, L1C, D1C, S1C;
	double C1C, L1C, D1C, S1C, C2W, L2W, S2W, C2X, L2X, S2X, C5X, L5X, S5X;
};

struct SP3_FRAME {
	double X, Y, Z, dt;
};

struct LOC_FRAME {
	double X, Y, Z;
};

struct INX_FRAME {
	int map_index;
	int VTECs[INX_LAT_BINS][INX_LON_BINS];
};

struct BRD_FRAME_GPS {
	// PRN / EPOCH / SV CLK
	GPSTime toc;
	double sv_clock_bias;
	double sv_clock_drift;
	double sv_clock_drift_rate;

	// BROADCAST ORBIT - 1
	double idoe_issue_of_data;
	double crs;
	double delta_n;
	double m0;

	// BROADCAST ORBIT - 2
	double cuc;
	double eccentricity;
	double cus;
	double sqrt_a;

	// BROADCAST ORBIT - 3
	double toe;
	double cic;
	double OMEGA;
	double cis;

	// BROADCAST ORBIT - 4
	double i0;
	double crc;
	double omega;
	double omega_dot;

	// BROADCAST ORBIT - 5
	double idot;
	double codes_on_l2;
	double gpsweek;
	double l2_pdata_flag;

	// BROADCAST ORBIT - 6
	double sv_accuracy;
	double sv_health;
	double tgd;
	double iodc_issue_of_data;

	// BROADCAST ORBIT - 7
	double trans_time;
	double fit_interval;
};

OBS_FRAME OBS[MAX_SAT_NUM];
BRD_FRAME_GPS BRD_GPS[MAX_SAT_NUM];
SP3_FRAME SP3[MAX_SAT_NUM];
LOC_FRAME SAT[MAX_SAT_NUM];



OBS_FRAME * S_OBS[MAX_SAT_NUM];
BRD_FRAME_GPS * S_BRD_GPS[MAX_SAT_NUM];
SP3_FRAME * S_SP3[MAX_SAT_NUM];
LOC_FRAME * S_SAT[MAX_SAT_NUM];

int SVN_LIST[MAX_SAT_NUM] = { -1 };

FILE * ofp = fopen(INPUT_OBS, "r");
FILE * nfp = fopen(INPUT_GPS_BRD, "r");

bool obs_available[MAX_SAT_NUM] = { false };
bool brd_available[MAX_SAT_NUM] = { false };
bool solve_available[MAX_SAT_NUM] = { false };

char line_buffer[MAX_LINE] = "";
char SVN[4] = "";
char nav_buffer[20] = "";

UTC utc_obs = UTC::current_utc();
UTC utc_brd = UTC();
GPSTime gpst(0, 0);
int sat_num = 0;



double ion_param[8] = {0};
double sat_ek[MAX_SAT_NUM] = { 0 };
double sat_elevation[MAX_SAT_NUM] = { 0 };
double sat_azimuth[MAX_SAT_NUM] = { 0 };
double correction_values[MAX_SAT_NUM] = { 0 };
Matrix * r_hat = NULL;


int brd_prn = 0;

double solution[4] = { -1641945.2000, -3664804.1000, 4940009.3000, 0};

inline double _fastcall get_atan(double z, double y)
{
	double x = 0;
	if (z == 0)x = M_PI / 2;
	else if (y == 0)x = M_PI;
	else {
		x = atan(fabs(y / z));
		if ((y > 0) && (z < 0))x = M_PI - x;
		else if ((y < 0) && (z < 0))x = M_PI + x;
		else if ((y < 0) && (z > 0))x = 2 * M_PI - x;
	}
	return x;
}

bool skip_obs_header(FILE * fp)
{
	while (!feof(fp))
	{
		fgets(line_buffer, MAX_LINE, fp);
		if (strncmp(line_buffer + 60, END_OF_HEADER, 13) == 0)
			return true;
	}
	return false;
}

void fetch_obs(FILE * fp)
{
	int health = 0;

	fgets(line_buffer, MAX_LINE, fp);
	sscanf(line_buffer, "> %d %d %d %d %d %lf %d %d",
		&utc_obs.year, &utc_obs.month, &utc_obs.date, &utc_obs.hour, &utc_obs.minute, &utc_obs.sec, &health, &sat_num);
	utc_obs.year = utc_obs.year - 2000;

	for (int i = 0; i < sat_num; i++)
	{
		fgets(SVN, 4, fp);
		if (SVN[0] == GPS_FLAG)
		{
			int svn = atoi(SVN + 1);
			int index = 0;
			double * ptr = (double*)(OBS + svn - 1);
			while (fscanf(fp, "%lf", ptr + index) && !feof(fp)) index++;
			obs_available[svn - 1] = true;
		}
		else if (SVN[1] == GLO_FLAG)
		{
			int svn = atoi(SVN + 1) + GLO_OFFSET;
			int index = 0;
			double * ptr = (double*)(OBS + svn - 1);
			while (fscanf(fp, "%lf", ptr + index) && !feof(fp)) index++;
			obs_available[svn - 1] = true;
		}
		else {
			fgets(line_buffer, MAX_LINE, fp);
		}
	}
}

void brd_grab_time(FILE * fp)
{
	if (feof(fp)) throw - 1;
	double sec;
	fscanf(fp, "%d%d%d%d%d%d%lf",
		&brd_prn,
		&utc_brd.year,
		&utc_brd.month,
		&utc_brd.date,
		&utc_brd.hour,
		&utc_brd.minute,
		&sec
	);
	utc_brd.sec = (int)round(sec);
}

bool skip_brd_header(FILE * fp)
{
	for (int i = 0; i < 4; i++)
		fgets(line_buffer, MAX_LINE, fp);

	for (int i = 0; i < 60; i++)
		if (line_buffer[i] == 'D')line_buffer[i] = 'E';

	sscanf(line_buffer, "%lf%lf%lf%lf",
		ion_param, ion_param + 1, ion_param + 2, ion_param + 3);

	fgets(line_buffer, MAX_LINE, fp);

	for (int i = 0; i < 60; i++)
		if (line_buffer[i] == 'D')line_buffer[i] = 'E';

	sscanf(line_buffer, "%lf%lf%lf%lf",
		ion_param + 4, ion_param + 5, ion_param + 6, ion_param + 7);

	if (skip_obs_header(fp))
		brd_grab_time(fp);
	else return false;
	return true;
}

inline double _fastcall extract_double(const char * pointee)
{
	if (strlen(pointee) < 19) return 0;
	strncpy(nav_buffer, pointee, 19);
	nav_buffer[15] = 'E';
	return atof(nav_buffer);
}

void fetch_gps_brd(FILE * fp)
{
	while (utc_obs.larger_than(utc_brd) && !feof(fp))
	{
		brd_prn--;
		// PRN / EPOCH / SV CLK
		if (fgets(line_buffer, MAX_LINE, fp) == NULL) throw - 1;
		BRD_GPS[brd_prn].toc = GPSTime(utc_brd);
		BRD_GPS[brd_prn].sv_clock_bias = extract_double(line_buffer);
		BRD_GPS[brd_prn].sv_clock_drift = extract_double(line_buffer + 19);
		BRD_GPS[brd_prn].sv_clock_drift_rate = extract_double(line_buffer + 38);

		// BROADCAST ORBIT - 1
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD_GPS[brd_prn].idoe_issue_of_data = extract_double(line_buffer);
		BRD_GPS[brd_prn].crs = extract_double(line_buffer + 19);
		BRD_GPS[brd_prn].delta_n = extract_double(line_buffer + 38);
		BRD_GPS[brd_prn].m0 = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 2
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD_GPS[brd_prn].cuc = extract_double(line_buffer);
		BRD_GPS[brd_prn].eccentricity = extract_double(line_buffer + 19);
		BRD_GPS[brd_prn].cus = extract_double(line_buffer + 38);
		BRD_GPS[brd_prn].sqrt_a = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 3
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD_GPS[brd_prn].toe = extract_double(line_buffer);
		BRD_GPS[brd_prn].cic = extract_double(line_buffer + 19);
		BRD_GPS[brd_prn].OMEGA = extract_double(line_buffer + 38);
		BRD_GPS[brd_prn].cis = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 4
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD_GPS[brd_prn].i0 = extract_double(line_buffer);
		BRD_GPS[brd_prn].crc = extract_double(line_buffer + 19);
		BRD_GPS[brd_prn].omega = extract_double(line_buffer + 38);
		BRD_GPS[brd_prn].omega_dot = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 5
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD_GPS[brd_prn].idot = extract_double(line_buffer);
		BRD_GPS[brd_prn].codes_on_l2 = extract_double(line_buffer + 19);
		BRD_GPS[brd_prn].gpsweek = extract_double(line_buffer + 38);
		BRD_GPS[brd_prn].l2_pdata_flag = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 6
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD_GPS[brd_prn].sv_accuracy = extract_double(line_buffer);
		BRD_GPS[brd_prn].sv_health = extract_double(line_buffer + 19);
		BRD_GPS[brd_prn].tgd = extract_double(line_buffer + 38);
		BRD_GPS[brd_prn].iodc_issue_of_data = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 7
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD_GPS[brd_prn].trans_time = extract_double(line_buffer);
		BRD_GPS[brd_prn].fit_interval = extract_double(line_buffer + 19);

		brd_available[brd_prn] = true;
		brd_grab_time(fp);
	}
}

double distance(const double * p1, const double * p2, int dim = 3)
{
	double tot = 0;
	for (int i = 0; i < dim; i++)
	{
		tot += (p1[i] - p2[i]) * (p1[i] - p2[i]);
	}
	return sqrt(tot);
}

void elevation_and_azimuth(LOC_FRAME * satellite, int index)
{
	double dpos[3] = { 0 };
	double ori[3]{ 0,0,0 };
	dpos[0] = satellite->X - station_xyz[0];
	dpos[1] = satellite->Y - station_xyz[1];
	dpos[2] = satellite->Z - station_xyz[2];

	double user_distance_to_earth = distance(station_xyz, center_earth);
	double mod = sqrt(dpos[0] * dpos[0] + dpos[1] * dpos[1] + dpos[2] * dpos[2]);
	if (fabs(user_distance_to_earth * mod < 1.0)) {
		sat_elevation[index] = M_PI_2;
	}
	else {
		double m = dpos[0] * station_xyz[0] + dpos[1] * station_xyz[1] + dpos[2] * station_xyz[2];
		double n = m / (mod * user_distance_to_earth);
		sat_elevation[index] = M_PI_2 - acos(n);
	}

	double N = -s_sinB * s_cosL * dpos[0] - s_sinB * s_sinL * dpos[1] + s_cosB * dpos[2];
	double E = -s_sinL * dpos[0] + s_cosL * dpos[1];
	sat_azimuth[index] = atan2(E, N);
}

bool brdc_satell_gps(int index)
{
	BRD_FRAME_GPS * f = S_BRD_GPS[index];
	double n0 = sqrt(mu) / pow(f->sqrt_a, 3);
	double n = n0 + f->delta_n;

	double t = gpst.sec - S_OBS[index]->C1C / LIGHT_SPEED;
	double tk = t - f->toe;
	if (fabs(tk) > 7200)
		return false;
	if (tk > 302400)
		tk -= 604800;
	else if (tk < -302400)
		tk += 604800;

	double Mk = f->m0 + n * tk;

	sat_ek[index] = Mk;
	double Ek2 = 0;
	while (1)
	{
		Ek2 = Mk + f->eccentricity * sin(sat_ek[index]);
		if (fabs(sat_ek[index] - Ek2) <= 1.0e-12)break;
		sat_ek[index] = Ek2;
	}

	double sqt_1_e2 = sqrt(1 - pow(f->eccentricity, 2));
	double cosfk = (cos(sat_ek[index]) - f->eccentricity) / (1 - f->eccentricity * cos(sat_ek[index]));
	double sinfk = (sqt_1_e2 * sin(sat_ek[index])) / (1 - f->eccentricity * cos(sat_ek[index]));
	double fk = get_atan((cos(sat_ek[index]) - f->eccentricity), sqt_1_e2 * sin(sat_ek[index]));

	double faik = fk + f->omega;

	double cos_2_faik = cos(2 * faik);
	double sin_2_faik = sin(2 * faik);
	double su = f->cuc * cos_2_faik + f->cus * sin_2_faik;
	double sr = f->crc * cos_2_faik + f->crs * sin_2_faik;
	double si = f->cic * cos_2_faik + f->cis * sin_2_faik;

	double uk = faik + su;
	double rk = f->sqrt_a * f->sqrt_a * (1 - f->eccentricity * cos(sat_ek[index])) + sr;
	double ik = f->i0 + si + f->idot * tk;

	double xk = rk * cos(uk);
	double yk = rk * sin(uk);
	double zk = 0;

	double L = f->OMEGA + (f->omega_dot - we) * tk - we * f->toe;

	S_SAT[index]->X = xk * cos(L) - yk * cos(ik) * sin(L);
	S_SAT[index]->Y = xk * sin(L) + yk * cos(ik) * cos(L);
	S_SAT[index]->Z = yk * sin(ik);

	return true;
}

void corrections_gps(int index)
{
	double dtc = gpst.minus(&S_BRD_GPS[index]->toc) - S_OBS[index]->C1C / LIGHT_SPEED;
	double s = S_BRD_GPS[index]->sv_clock_bias
		+ S_BRD_GPS[index]->sv_clock_drift      * dtc
		+ S_BRD_GPS[index]->sv_clock_drift_rate * dtc * dtc;

	double r = (R1 * S_BRD_GPS[index]->eccentricity * sin(sat_ek[index]) * S_BRD_GPS[index]->sqrt_a);


	correction_values[SVN_LIST[index] - 1] -= r * LIGHT_SPEED;                                  //相对论
	correction_values[SVN_LIST[index] - 1] -= S_BRD_GPS[index]->tgd * LIGHT_SPEED;              //群延迟
	correction_values[SVN_LIST[index] - 1] += s * LIGHT_SPEED;									//星钟差

	S_OBS[index]->C1C += correction_values[SVN_LIST[index] - 1];

	double dA = we * (S_OBS[index]->C1C / LIGHT_SPEED - s + r);
	double Xc = S_SAT[index]->Y * dA;
	double Yc = S_SAT[index]->X * dA;

	S_SAT[index]->X += Xc;
	S_SAT[index]->Y -= Yc;
}

void overall_check()
{
	sat_num = 0;
	gpst = GPSTime(utc_obs);

	for (int i = 0; i < MAX_SAT_NUM; i++)
	{
		if (brd_available[i] && obs_available[i])
		{
			if (OBS[i].C1C == 0)continue;
			for (int j = 0; j < ban_svn_amount; j++)
				if (ban_svn[j] == i + 1)
					goto cont;

			S_OBS[sat_num] = OBS + i;
			S_BRD_GPS[sat_num] = BRD_GPS + i;
			S_SAT[sat_num] = SAT + i;
			

			if (!brdc_satell_gps(sat_num)) continue;

			elevation_and_azimuth(SAT + i, i);

			if (sat_elevation[i] <= ELEV_THRESH)
				continue;

			SVN_LIST[sat_num] = i + 1;
			corrections_gps(sat_num);
			solve_available[i] = true;
			sat_num++;
		}
	cont:;
	}
}

void erase()
{
	memset(OBS, 0, sizeof(OBS_FRAME) * MAX_SAT_NUM);
	memset(SP3, 0, sizeof(SP3_FRAME) * MAX_SAT_NUM);
	memset(SAT, 0, sizeof(LOC_FRAME) * MAX_SAT_NUM);
	memset(S_OBS, 0, sizeof(OBS_FRAME*) * MAX_SAT_NUM);
	memset(S_SP3, 0, sizeof(SP3_FRAME*) * MAX_SAT_NUM);
	memset(S_SAT, 0, sizeof(LOC_FRAME*) * MAX_SAT_NUM);
	memset(S_BRD_GPS, 0, sizeof(BRD_FRAME_GPS*) * MAX_SAT_NUM);
	memset(SVN_LIST, 0, sizeof(int) * MAX_SAT_NUM);
	memset(obs_available, 0, sizeof(bool) * MAX_SAT_NUM);
	memset(solve_available, 0, sizeof(bool) * MAX_SAT_NUM);
	memset(line_buffer, 0, sizeof(char) * MAX_LINE);
	memset(SVN, 0, sizeof(char) * 4);
	memset(nav_buffer, 0, sizeof(char) * 20);
	sat_num = 0;
	memset(sat_ek, 0, sizeof(double) * MAX_SAT_NUM);
	memset(sat_elevation, 0, sizeof(double) * MAX_SAT_NUM);
	memset(sat_azimuth, 0, sizeof(double) * MAX_SAT_NUM);
	memset(correction_values, 0, sizeof(double) * MAX_SAT_NUM);
}

bool solve()
{
	Matrix * L = malloc_mat(sat_num, 1);
	Matrix * A = malloc_mat(sat_num, 4);
	Matrix * Cl = malloc_mat(sat_num, sat_num);
	Matrix * Q = malloc_mat(4, 4);
	Matrix * δ = malloc_mat(4, 1);
	Matrix * Cx = malloc_mat(4, 4);

	double * DX0 = (double*)alloca(sat_num * sizeof(double));
	double * DY0 = (double*)alloca(sat_num * sizeof(double));
	double * DZ0 = (double*)alloca(sat_num * sizeof(double));
	double * S = (double*)alloca(sat_num * sizeof(double));

	double last_solution[4] = { 0,0,0,0 };
	//memcpy(solution, station_xyz, sizeof(double) * 4);


	for (int i = 0; i < LS_MAX_ITER; i++)
	{
		memcpy(last_solution, solution, sizeof(double) * 4);

		for (int j = 0; j < sat_num; j++)
		{
			DX0[j] = S_SAT[j]->X - solution[0];
			DY0[j] = S_SAT[j]->Y - solution[1];
			DZ0[j] = S_SAT[j]->Z - solution[2];
			S[j] = sqrt(DX0[j] * DX0[j] + DY0[j] * DY0[j] + DZ0[j] * DZ0[j]);
		}

		// get A, L, Cl matrices.
		for (int j = 0; j < sat_num; j++)
		{
			double sinval = sin(sat_elevation[SVN_LIST[j] - 1]);
			Cl->data[j][j] = OBS_SIG0 * OBS_SIG0 / (sinval * sinval);

			// for common spp
			L->data[j][0] = S_OBS[j]->C1C - S[j] - solution[3];

			// A matrix
			A->data[j][0] = -DX0[j] / S[j];  // for X
			A->data[j][1] = -DY0[j] / S[j];  // for Y
			A->data[j][2] = -DZ0[j] / S[j];  // for Z
			A->data[j][3] = 1;               // for cdt
		}

		LMS(L, A, Cl, δ, Q, r_hat, Cx);

		for (int j = 0; j < 4; j++)
			solution[j] += δ->data[j][0];

		// If converged
		if (distance(last_solution, solution, 4) <= LS_CONV_THRESH) {
			return true;
		}
	}

	return true;
}

FILE * out1 = fopen("out1.txt", "w");
FILE * out2 = fopen("out2.txt", "w");

void output()
{
	fprintf(out1, "%4d%15.2f%20.5lf%20.5lf%20.5lf%20.4lf%10d\n",
		gpst.week, gpst.sec, solution[0], solution[1], solution[2], solution[3], sat_num);

	for (int i = 0; i < sat_num; i++)
	{
		fprintf(out2, "%6.3lf ", r_hat->data[i][0]);
	}

	fprintf(out2, "\n");
}

int main()
{
	
	skip_obs_header(ofp);
	skip_brd_header(nfp);

	try {
		while (!feof(ofp)/* && !feof(nfp)*/)
		{
			fetch_obs(ofp);
			fetch_gps_brd(nfp);
			overall_check();

			if (fabs(gpst.sec - 430053.20) <= 0.01)
			{
				int i = 0;
			}
			if (sat_num >= SAT_NUM_THRESH)
			{
				if (solve())
				{
					output();
				}
			}
			erase();
		}
	}
	catch (...){}

	_fcloseall();
}
