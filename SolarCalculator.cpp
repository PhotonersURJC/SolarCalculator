
#include <stdio.h>
#include "spa.h"  //include the SPA header file
#include <iostream>
#include <vector>
#include <time.h>
#include <string>

# define M_PI           3.14159265358979323846

#define timeExp 1.3405  //time exponent of Hom's Law
#define k1 7.709199746 //s^(-timeExp)
#define k2 1.1437e+28
#define Ea1 7417.744722 //both are K-1  (Ea/R)
#define Ea2 22780.93448
#define alpha 1.60547561
#define rAM 707.8888889


//#CHANGEKINETICS: This function returns the kinetic constant of Hom's law as a function of Temp and TIR
double logKineticSODIS(double Temperature, double Radiation)
{
	return k1 * pow(Radiation, alpha)*exp(-Ea1 / (273.0 + Temperature) + k2 * exp(-Ea2 / (273.0 + Temperature)));
}


std::vector<double> getInstantTemp(double minTemp, double maxTemp, double timestep, int nSteps)
{
	double level = 0.25*timestep*nSteps;
	double xFactor = 7.27221E-05; //2*PI/día en segundos
	double yFactor = 0.5*(maxTemp - minTemp);
	double aveT = 0.5*(maxTemp + minTemp);
	std::vector<double> result;
	for (int cT = 0; cT < nSteps; cT++)
		result.push_back(aveT+yFactor*sin(xFactor*(cT*timestep-level)));
	return result;
}

double calculateDay(int day, int month, double latitude, double elevation, std::vector<double> &AM, std::vector<double> &cosZ, double &sunDistance, double temperature = 30.0)
{
	spa_data spa;  //declare the SPA structure
	int result;
	double min, sec;

	//enter required input values into SPA structure

	spa.year = 2018;
	spa.month = month;
	spa.day = day;
	spa.hour = 12;
	spa.minute = 00;
	spa.second = 00;
	spa.timezone = 0.0;
	spa.delta_ut1 = 0;
	spa.delta_t = 67;
	spa.longitude = 20.0;
	spa.latitude = latitude;
	spa.elevation = elevation;
	spa.pressure = 820;
	spa.temperature = temperature;
	spa.slope = 0.0;
	spa.azm_rotation = 0.0;
	spa.atmos_refract = 0.5667;
	spa.function = SPA_ALL;

	//call the SPA calculate function and pass the SPA structure
	result = spa_calculate(&spa);
	sunDistance = 1.0 / spa.r;
	min = 60.0*(spa.sunrise - (int)spa.sunrise);
	sec = 60.0*(min - (int)min);
	spa.hour = (int)spa.sunrise;
	spa.minute = (int)min;
	spa.second = sec;
	double timestep = (spa.sunset - spa.sunrise)*40.0;
	for (int i = 0; i < 90; i++)
	{
		result = spa_calculate(&spa);
		if (result)
			break;
		double a = cos(deg2rad(spa.zenith));
		double c = elevation / 9000.0;
		double b = sqrt((rAM + c)*(rAM + c)*a*a + (2 * rAM + 1 + c)) - (rAM + c)*a;
		//double b = (1.002432*a*a + 0.148386*a + 0.0096467) / (a*a*a + 0.149864*a*a + 0.0102963*a + 0.000303978);
		AM.push_back(b);
		cosZ.push_back(deg2rad(spa.zenith));
		spa.second += timestep;
		while (spa.second >= 60.0)
		{
			spa.minute++;
			spa.second -= 60.0;
		}
		while (spa.minute > 59)
		{
			spa.minute -= 60;
			spa.hour += 1;
		}
	}
	return timestep;
}

double calculateDay(int day, int month, double latitude, double elevation, double UTC, std::vector<double> &AM, std::vector<double> &cosZ, std::vector<double> &zenith,
	std::vector<double> &azimuth, double &sunDistance, double* hours, double temperature = 30.0)
{
	spa_data spa;  //declare the SPA structure
	int result;
	double min, sec;

	//enter required input values into SPA structure

	spa.year = 2018;
	spa.month = month;
	spa.day = day;
	spa.hour = 12;
	spa.minute = 00;
	spa.second = 00;
	spa.timezone = UTC;
	spa.delta_ut1 = 0;
	spa.delta_t = 67;
	spa.longitude = 20.0;
	spa.latitude = latitude;
	spa.elevation = elevation;
	spa.pressure = 820;
	spa.temperature = temperature;
	spa.slope = 0.0;
	spa.azm_rotation = 0.0;
	spa.atmos_refract = 0.5667;
	spa.function = SPA_ALL;

	//call the SPA calculate function and pass the SPA structure
	result = spa_calculate(&spa);
	sunDistance = 1.0 / spa.r;
	min = 60.0*(spa.sunrise - (int)spa.sunrise);
	sec = 60.0*(min - (int)min);
	spa.hour = (int)spa.sunrise;
	spa.minute = (int)min;
	spa.second = sec;
	double timestep = (spa.sunset - spa.sunrise)*40.0;
	for (int i = 0; i < 90; i++)
	{
		result = spa_calculate(&spa);
		if (result)
			break;
		zenith.push_back(spa.zenith);
		azimuth.push_back(spa.azimuth);
		double a = cos(deg2rad(spa.zenith));
		double c = elevation / 9000.0;
		double b = sqrt((rAM + c)*(rAM + c)*a*a + (2 * rAM + 1 + c)) - (rAM + c)*a;
		//double b = (1.002432*a*a + 0.148386*a + 0.0096467) / (a*a*a + 0.149864*a*a + 0.0102963*a + 0.000303978);
		AM.push_back(b);
		cosZ.push_back(deg2rad(spa.zenith));
		spa.second += timestep;
		while (spa.second >= 60.0)
		{
			spa.minute++;
			spa.second -= 60.0;
		}
		while (spa.minute > 59)
		{
			spa.minute -= 60;
			spa.hour += 1;
		}
	}
	*hours = spa.sunrise;
	return timestep;
}




double calculateDayRadTemp
(
	int day,
	int month,
	double latitude,
	double elevation,
	double cloud,
	double minTemp,
	double maxTemp,
	double UTC,
	std::vector<double> &TIR,
	std::vector<double> &instantT,
	std::vector<double> &zenith,
	std::vector<double> &azimuth,
	std::vector<double> &difFraction,
	double * WL,
	double * I0,
	double *kAtm,
	int lines,
	std::vector<double> &hours
)
{
	std::vector<double> daylyAM, cosZ;
	double sunDistance = 1.0;
	double sunrise = 0.0;
	double timestep = calculateDay(day, month, latitude, elevation, UTC, daylyAM, cosZ, zenith, azimuth, sunDistance,&sunrise, 0.5*(maxTemp + minTemp));
	double dayVal = 0.0;
	for (int cT = 0; cT < daylyAM.size(); cT++)
	{
		double currentAM = daylyAM[cT];
		double currentI = 0;
		for (int cWL = 20; cWL < lines - 10; cWL++)
			currentI += (WL[cWL] - WL[cWL - 1])*0.5*(I0[cWL] * exp(-currentAM * kAtm[cWL]) + I0[cWL - 1] * exp(-currentAM * kAtm[cWL - 1]));
		currentI *= sunDistance;
		//TIR value from MTI
		double DIC = 2.0*currentI*(0.182 - 9.26e-5*cloud + 6.77e-4*sin(0.06981317*latitude + 1.2));
		currentI *= (1.09090296 - 0.00536302*cloud + 0.00124499*sin(0.06981317*latitude + 1.2));
		difFraction.push_back(DIC / currentI);
		TIR.push_back(currentI);
		hours.push_back(sunrise+timestep*(double)cT/3600.0);
	}
	instantT = getInstantTemp(minTemp, maxTemp, timestep, daylyAM.size());
	return timestep;
}

int getCommonSpecter(double* &WL, double* &I0, double* &kAtm)
{
	char *p; char line_c[200]; int i = 0; FILE * myfile;
	fopen_s(&myfile, "kAM.txt","r");
	int lines = 1;
	while (fgets(line_c, 200, myfile))
		lines++;
	double* WL2 = new double[lines];
	double* I02 = new double[lines];
	double* kAtm2 = new double[lines];
	fclose(myfile);
	fopen_s(&myfile, "kAM.txt","r");
	while (fgets(line_c, 200, myfile))
	{
		std::string line(line_c);
		std::string WLstr = line.substr(0, line.find("nm"));
		WL2[i] = strtod(WLstr.c_str(), &p);
		std::string I0str = line.substr(line.find("nm") + 2, line.find("W"));
		I02[i] = strtod(I0str.c_str(), &p);
		std::string kAtmStr = line.substr(line.find("W") + 1);
		kAtm2[i] = strtod(kAtmStr.c_str(), &p);
		i++;
	}
	fclose(myfile);
	WL = WL2; I0 = I02; kAtm = kAtm2;
	return lines;
}

void getAllSodisData()
{
	char buffer[200];
	std::string s;
	FILE *monthBothAvg,*totalFile;
	fopen_s(&totalFile,"totalSodis", "w");
	for (int month = 0; month < 12; month++)
	{
		_itoa_s(month + 1, buffer, 10);
		s = "month" + std::string(buffer);
		fopen_s(&monthBothAvg, s.c_str(), "r");
		while (fgets(buffer, 200, monthBothAvg) != NULL)
			fprintf(totalFile, "%s", buffer);
	}

	fclose(monthBothAvg);
	fclose(totalFile);
}


double simpleInterp(double x, double x1, double x2, double y1, double y2)
{
	return y1 + (y2 - y1)*(x - x1)/(x2 - x1);
}

bool * readPositions(int & totalPositions)
{
	FILE * readFile;
	fopen_s(&readFile, "AvailableBoolHD", "rb");
	bool * positions = (bool *)malloc(sizeof(bool) * 720 * 360);
	fread(positions, sizeof(bool), 720 * 360, readFile);
	fclose(readFile);
	totalPositions = 0;
	for(int i = 0; i < 720*360; i++)
		if(positions[i])
			totalPositions++;
		return positions;
}
double * readElevations(int totalPositions)
{
	short* shortElev = (short*)malloc(sizeof(short)*totalPositions);
	FILE * readFile;
	fopen_s(&readFile, "elevHD", "rb");
	fread(shortElev, sizeof(short), totalPositions, readFile);
	fclose(readFile);
	double* elevations = (double*)malloc(sizeof(double)*totalPositions);
	for(int i = 0; i < totalPositions; i++)
		elevations[i]=shortElev[i]/100.0;
	free(shortElev);
	return elevations;
}



double * getAtemporalValues(double latitude, double longitude)
{
	//reading elevation and timezone
	//returns a double array of 13 values: 0-11 are deltaUTC values for each month. 12 is elevation
	//if nullpoint (sea), returns -15.0, to avoid extra computation
	FILE * readFile;
	fopen_s(&readFile, "AvailableIntHD", "rb");
	int * positions = (int *)malloc(sizeof(int) * 720 * 360);
	fread(positions, sizeof(int), 720 * 360, readFile);
	fclose(readFile);
	/* 0 --- 1
	   |     |
	   2-----3
	*/
	int latI = floor(2.0*latitude + 180.0);
	int longI = floor(2.0*longitude + 360.0);
	int pos[4];
	double lat1 = latI * 0.5 - 90.0;
	double lat2 = lat1 + 0.5;
	double long1 = longI * 0.5 - 180.0;
	double long2 = long1 + 0.5;
	pos[0] = positions[latI * 720 + longI];
	pos[1] = positions[latI * 720 + longI + 1];
	pos[2] = positions[latI * 720 + longI + 720];
	pos[3] = positions[latI * 720 + longI + 721];
	free(positions);

	double * result = (double *)malloc(sizeof(double) * 13);
	if (pos[0] == -1 && pos[1] == -1 && pos[2] == -1 && pos[3] == -1)
	{
		for (int i = 0; i < 13; i++)
			result[i] = -15.0;
		return result;
	}
	short * elevRead = (short*)malloc(sizeof(short) * (pos[3] + 1));
	short elev[4];
	fopen_s(&readFile, "elevHD", "rb");
	fread(elevRead, sizeof(short), (pos[3] + 1), readFile);
	elev[0] = elevRead[pos[0]];
	elev[1] = elevRead[pos[1]];
	elev[2] = elevRead[pos[2]];
	elev[3] = elevRead[pos[3]];
	fclose(readFile);

	double elev1 = 0.0;
	if (pos[0] == -1 && pos[1] == -1)
		elev1 = 0.0;
	else if (pos[0] == -1 || pos[1] == -1)
		elev1 = (pos[0] == -1) ? elev[1] : elev[0];
	else
		elev1 = simpleInterp(longitude, long1, long2, (double)elev[0] / 100.0, (double)elev[1] / 100.0);
	double elev2 = 0.0;
	if (pos[2] == -1 && pos[3] == -1)
		elev2 = 0.0;
	else if (pos[2] == -1 || pos[3] == -1)
		elev2 = (pos[2] == -1) ? elev[2] : elev[3];
	else
		elev2 = simpleInterp(longitude, long1, long2, (double)elev[2] / 100.0, (double)elev[3] / 100.0);

	if (pos[0] == -1 && pos[1] == -1)
		result[12] = elev2;
	else if (pos[2] == -1 && pos[3] == -1)
		result[12] = elev1;
	else
		result[12] = simpleInterp(latitude, lat1, lat2, elev1, elev2);

	char * UTCread = (char*)malloc(sizeof(char) * (pos[0] + 1));
	char UTCval;
	fopen_s(&readFile, "UTCmapDST", "rb");
	fread(UTCread, sizeof(char), (pos[0] + 1), readFile);
	UTCval = UTCread[pos[0]];

	//en horas
	double deltaUTC = 0.5* (double)(UTCval / 2 - 24);
	for (int i = 0; i < 12; i++)
		result[i] = deltaUTC;
	if (UTCval % 2)
	{
		if (latitude > 0.0)
			for (int i = 3; i < 10; i++)
				result[i] += 1.0;
		else
		{
			for (int i = 0; i < 3; i++)
				result[i] += 1.0;
			for (int i = 8; i < 12; i++)
				result[i] += 1.0;
		}
	}
	return result;
}

double * readMetheo(int totalPositions, const char * fileName)
{
	short * shortMet = (short*)malloc(sizeof(short)*totalPositions);
	FILE * readFile;
	fopen_s(&readFile, fileName, "rb");
	fread(shortMet, sizeof(short), totalPositions, readFile);
	fclose(readFile);
	double * metheoData = (double*)malloc(sizeof(double)*totalPositions);
	for(int i = 0; i < totalPositions; i++)
		metheoData[i]=shortMet[i]/100.0;
	free(shortMet);
	return metheoData;
}

double *  interpolateMonthData(double latitude, double longitude, int month)
{
	FILE * readFile;
	fopen_s(&readFile, "AvailableIntHD", "rb");
	int * positions = (int *)malloc(sizeof(int) * 720 * 360);
	fread(positions, sizeof(int), 720 * 360, readFile);
	fclose(readFile);
	/* 0 --- 1
	   |     |
	   2-----3
	*/
	int latI = floor(2.0*latitude + 180.0);
	int longI = floor(2.0*longitude + 360.0);
	int pos[4];
	double lat1 = latI * 0.5 - 90.0;
	double lat2 = lat1 + 0.5;
	double long1 = longI * 0.5 - 180.0;
	double long2 = long1 + 0.5;
	pos[0] = positions[latI * 720 + longI];
	pos[1] = positions[latI * 720 + longI + 1];
	pos[2] = positions[latI * 720 + longI + 720];
	pos[3] = positions[latI * 720 + longI + 721];
	free(positions);
	
	double * result = (double *)malloc(sizeof(double) * 3);
	if (pos[0] == -1 && pos[1] == -1 && pos[2] == -1 && pos[3] == -1)
	{
		result[0] = 0.0; result[1] = 0.0; result[2] = 0.0;
		return result;
	}

	short * monthRead = (short*)malloc(sizeof(short) * (pos[3]+1));
	fopen_s(&readFile, (std::string("minTHD") + std::to_string(month)).c_str(), "rb");
	fread(monthRead, sizeof(short), (pos[3]+1), readFile);
	short minT[4], maxT[4], cloud[4];
	minT[0] = monthRead[pos[0]];
	minT[1] = monthRead[pos[1]];
	minT[2] = monthRead[pos[2]];
	minT[3] = monthRead[pos[3]];
	fclose(readFile);
		
	double TLat1 = 0.0;
	if (pos[0] == -1 && pos[1] == -1)
		TLat1 = 0.0;
	else if (pos[0] == -1 || pos[1] == -1)
		TLat1 = (pos[0] == -1) ? minT[1] : minT[0];
	else
		TLat1 = simpleInterp(longitude, long1, long2, (double)minT[0] / 100.0, (double)minT[1] / 100.0);
	double TLat2 = 0.0;
	if (pos[2] == -1 && pos[3] == -1)
		TLat2 = 0.0;
	else if (pos[2] == -1 || pos[3] == -1)
		TLat2 = (pos[2] == -1) ? minT[2] : minT[3];
	else
		TLat2 = simpleInterp(longitude, long1, long2, (double)minT[2] / 100.0, (double)minT[3] / 100.0);
		
	if (pos[0] == -1 && pos[1] == -1)
		result[0] = TLat2;
	else if (pos[2] == -1 && pos[3] == -1)
		result[0] = TLat1;
	else
		result[0] = simpleInterp(latitude, lat1, lat2, TLat1, TLat2);

	fopen_s(&readFile, (std::string("maxTHD") + std::to_string(month)).c_str(), "rb");
	fread(monthRead, sizeof(short), (pos[3] + 1), readFile);
	maxT[0] = monthRead[pos[0]];
	maxT[1] = monthRead[pos[1]];
	maxT[2] = monthRead[pos[2]];
	maxT[3] = monthRead[pos[3]];
	fclose(readFile);

	if (pos[0] == -1 && pos[1] == -1)
		TLat1 = 0.0;
	else if (pos[0] == -1 || pos[1] == -1)
		TLat1 = (pos[0] == -1) ? maxT[1] : maxT[0];
	else
		TLat1 = simpleInterp(longitude, long1, long2, (double)maxT[0] / 100.0, (double)maxT[1] / 100.0);
	if (pos[2] == -1 && pos[3] == -1)
		TLat2 = 0.0;
	else if (pos[2] == -1 || pos[3] == -1)
		TLat2 = (pos[2] == -1) ? maxT[2] : maxT[3];
	else
		TLat2 = simpleInterp(longitude, long1, long2, (double)maxT[2] / 100.0, (double)maxT[3] / 100.0);

	if (pos[0] == -1 && pos[1] == -1)
		result[1] = TLat2;
	else if (pos[2] == -1 && pos[3] == -1)
		result[1] = TLat1;
	else
		result[1] = simpleInterp(latitude, lat1, lat2, TLat1, TLat2);

	fopen_s(&readFile, (std::string("cloudHD") + std::to_string(month)).c_str(), "rb");
	fread(monthRead, sizeof(short), (pos[3] + 1), readFile);
	cloud[0] = monthRead[pos[0]];
	cloud[1] = monthRead[pos[1]];
	cloud[2] = monthRead[pos[2]];
	cloud[3] = monthRead[pos[3]];
	fclose(readFile);
	free(monthRead);

	if (pos[0] == -1 && pos[1] == -1)
		TLat1 = 0.0;
	else if (pos[0] == -1 || pos[1] == -1)
		TLat1 = (pos[0] == -1) ? cloud[1] : cloud[0];
	else
		TLat1 = simpleInterp(longitude, long1, long2, (double)cloud[0] / 100.0, (double)cloud[1] / 100.0);
	if (pos[2] == -1 && pos[3] == -1)
		TLat2 = 0.0;
	else if (pos[2] == -1 || pos[3] == -1)
		TLat2 = (pos[2] == -1) ? cloud[2] : cloud[3];
	else
		TLat2 = simpleInterp(longitude, long1, long2, (double)cloud[2] / 100.0, (double)cloud[3] / 100.0);

	if (pos[0] == -1 && pos[1] == -1)
		result[2] = TLat2;
	else if (pos[2] == -1 && pos[3] == -1)
		result[2] = TLat1;
	else
		result[2] = simpleInterp(latitude, lat1, lat2, TLat1, TLat2);

	return result;
	
}


void generateTirChart(double latitude, double longitude, const char * fileName)
{
	//Obtiene las longitudes de onda y constantes de extincion para la atmosfera, en 1/AM
	double * WL, *I0, *kAtm;
	int lines = getCommonSpecter(WL, I0, kAtm);

	
	double * aTempValues = getAtemporalValues(latitude, longitude);
	if (aTempValues[0] == -15.0)
	{
		free(aTempValues);
		return;
	}
	
	//Obtengo los datos meteorólogicos medios de cada mes
	double ** monthValues = (double**)malloc(sizeof(double*) * 12);
	for (int i = 0; i < 12; i++)
		monthValues[i] = interpolateMonthData(latitude, longitude, i);
	char buffer[3];
	std::string s;
	FILE * resultsFile;
	fopen_s(&resultsFile, fileName, "w");
	int latI = 2.0*latitude + 180;
	int longI = 2.0*longitude + 360;
	int yearDay = 0;
	for (int month = 0; month < 12; month++)
	{
		//cálculo del número de días en el mes
		int daysInMonth = 31;
		if (month == 3 || month == 5 || month == 8 || month == 10)
			daysInMonth = 30;
		if (month == 1)
			daysInMonth = 28;
		for (int day = 1; day <= daysInMonth; day++)
		{
			yearDay++;
			//cálculo del std::vector de radiacion
			std::vector<double> daylyI, daylyT, daylyZenith, daylyAzimuth, hours, difFraction;
			double dayAccumulated = 0.0;
			double minT = monthValues[month][0];
			double maxT = monthValues[month][1];
			double cloudCoverage = monthValues[month][2];
			if (day < 15)
			{
				minT = simpleInterp(day, -15.0, 15.0, monthValues[(12+month - 1) % 12][0], minT);
				maxT = simpleInterp(day, -15.0, 15.0, monthValues[(12+month - 1) % 12][1], maxT);
				cloudCoverage = simpleInterp(day, -15.0, 15.0, monthValues[(12+month - 1) % 12][2], cloudCoverage);
			}
			if (day > 15)
			{
				minT = simpleInterp(day, 15.0, daysInMonth + 15, minT, monthValues[(month + 1) % 12][0]);
				maxT = simpleInterp(day, 15.0, daysInMonth + 15, maxT, monthValues[(month + 1) % 12][1]);
				cloudCoverage = simpleInterp(day, 15.0, daysInMonth + 15, cloudCoverage, monthValues[(month + 1) % 12][2]);
			}
			double timestep = calculateDayRadTemp
			(
				day,
				month + 1,
				latitude,
				aTempValues[12],
				cloudCoverage,
				minT,
				maxT,
				aTempValues[month],
				daylyI,
				daylyT,
				daylyZenith,
				daylyAzimuth,
				difFraction,
				WL,
				I0,
				kAtm,
				lines,
				hours
			);
			for (int i = 0; i < daylyI.size(); i++)
				fprintf(resultsFile, "%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", 
					month+1, day, yearDay, hours[i], daylyT[i], daylyI[i], difFraction[i], daylyZenith[i], daylyAzimuth[i]);
		}
	}
	fclose(resultsFile);
}

void disinfection()
{
	//Reading wavelength and atmospheric extinction coefficients (1/AM)
	double * WL, *I0, *kAtm;
	int lines = getCommonSpecter(WL, I0, kAtm);
	
	//reading terrain elevation
	int totalPositions;
	bool * positions = readPositions(totalPositions);
	double * elevations = readElevations(totalPositions);
	
	//reading metheorological data
	double ** minT = (double**)malloc(sizeof(double*)*12);
	double ** maxT = (double**)malloc(sizeof(double*)*12);
	double ** cloud = (double**)malloc(sizeof(double*)*12);
	for (int month = 0; month < 12; month++)
	{
		minT[month] = readMetheo(totalPositions, (std::string("minTHD") + std::to_string(month)).c_str());
		maxT[month] = readMetheo(totalPositions, (std::string("maxTHD") + std::to_string(month)).c_str());
		cloud[month] = readMetheo(totalPositions, (std::string("cloudTHD") + std::to_string(month)).c_str());
	}
	

	char buffer[3];
	std::string s;
	for (int month = 0; month < 12; month++)
	{
		//creating file to write results
		_itoa_s(month + 1, buffer, 10);
		FILE *monthBothAvg;
		s = "month" + std::string(buffer);
		fopen_s(&monthBothAvg, s.c_str(), "w");
		s = "15/" + std::string(buffer) + "/2019";
		
		int daysInMonth = 31;
		if (month == 3 || month == 5 || month == 8 || month == 10)
			daysInMonth = 30;
		if (month == 1)
			daysInMonth = 28;
		for (int latI = 68; latI < 314; latI++) //polar circles: from 47 to 314. 47 to 68 is empty
		{
			double latitude = 0.5*(double)(latI)-90.0;
			std::vector<double> minBatchTime(720, 0.0);
			std::vector<int> monthBatches(720, 0);
			for (int longI = 0; longI < 720; longI++)
			{
				//checking if point is available
				if (!positions[latI * 720 + longI])
					continue;
				if (maxTemp[month][latI*720+longI] < 4.0)
					continue;
			//loop along days in month
			for (int day = 1; day <= daysInMonth; day++)
			{
				std::vector<double> daylyI, daylyT;
				double timestep = calculateDayRadTemp
				(
					day,
					month + 1,
					latitude,
					elevations[latI*720+longI],
					(day<15)
						? simpleInterp(day, -15, 15, cloud[(month+11)%12][latI*720+longI], cloud[month][latI*720+longI])
						: simpleInterp(day, 15, daysInMonth+15, cloud[month][latI*720+longI], cloud[(month+1)%12][latI*720+longI]),
					(day<15)
						? simpleInterp(day, -15, 15, minT[(month+11)%12][latI*720+longI], minT[month][latI*720+longI])
						: simpleInterp(day, 15, daysInMonth+15, minT[month][latI*720+longI], minT[(month+1)%12][latI*720+longI]),
					(day<15)
						? simpleInterp(day, -15, 15, maxT[(month+11)%12][latI*720+longI], maxT[month][latI*720+longI])
						: simpleInterp(day, 15, daysInMonth+15, maxT[month][latI*720+longI], maxT[(month+1)%12][latI*720+longI]),
					daylyI,
					daylyT, 
					WL, 
					I0, 
					kAtm, 
					lines
				);

					//a vector of possible starting points of optimum batch to estimate minimum time
					std::vector<double> quickestBatch(daylyI.size(), 0.0);
					std::vector<double> qBtime(daylyI.size(), 100.0*daylyI.size()*timestep);
					//and a vector to store results
					double disinfection = 0.0;
					double elapsedTS = 0.0;
					double initPositiveTime = 0.0;
					for (int cT = 0; cT < daylyI.size(); cT++)
					{
						//at too low temperatures disinfection is neglected, and initTime is shifted towards the next step
						if (daylyT[cT] < 4.0)
						{
							initPositiveTime += timestep;
							continue;
						}
						
						//#CHANGEKINETICS: if temperature is between 4.0 ºC and 20.0 ºC, 20.0 ºC is used in the model
						//#CHANGEKINETICS: in the loop below, Hom's law is used to calculate the disinfection
						double instantKinetics = logKineticSODIS((daylyT[cT]>20.0?daylyT[cT]:20.0), daylyI[cT]);
						for (int cT2 = 0; cT2 <= cT; cT2++)
						{
							if (quickestBatch[cT2] == 4.0)
								continue;
							
							double addition = instantKinetics * (pow((cT - cT2 + 1.0)*timestep-initPositiveTime, timeExp)
								- pow((cT - cT2)*timestep- initPositiveTime, timeExp));
							if ((quickestBatch[cT2] + addition) >= 4.0)
							{
								qBtime[cT2] = pow(((4.0 - quickestBatch[cT2] + instantKinetics *
									pow((cT - cT2)*timestep- initPositiveTime, timeExp)) / instantKinetics), 1.0 / timeExp);
								quickestBatch[cT2] = 4.0;
								continue;
							}
							quickestBatch[cT2] += addition;
						}
						//#CHANGEKINETICS: again, Hom's law is used, in accumulated batches
						double addition = instantKinetics * (pow(elapsedTS + timestep, timeExp)
							- pow(elapsedTS, timeExp));
						if ((disinfection + addition) >= 4.0)
						{
							monthBatches[longI]++;
							double endTime = pow(((4.0 - disinfection + instantKinetics *
								pow(elapsedTS, timeExp)) / instantKinetics), 1.0 / timeExp);
							elapsedTS = elapsedTS + timestep - endTime;
							disinfection = instantKinetics * pow(elapsedTS, timeExp);
						}
						else
						{
							disinfection += addition;
							elapsedTS += timestep;
						}
					}

					double minTime = 100.0*daylyI.size()*timestep;
					for (int cT = 0; cT < daylyI.size(); cT++)
						if (qBtime[cT] < minTime)
							minTime = qBtime[cT];
					if (minTime > 90.0*daylyI.size()*timestep)
						minBatchTime[longI] += 0.0;
					else
						minBatchTime[longI] += (21600.0 / minTime);
				}
			}
			for (int longI = 0; longI < 720; longI++)
			{
				//compruebo si se descarta el dato
				if (!positions[latI * 720 + longI])
					continue;
				double currentBatches = (double)(monthBatches[longI]) / (double)(daysInMonth);
				double minTime = minBatchTime[longI] / (double)(daysInMonth);
				fprintf(monthBothAvg, "%f\t%f\t%s\t%f\t%f\n", latitude, 0.5*(double)(longI)-180.0, s.c_str(),  minTime, currentBatches);
			}
		}
		fclose(monthBothAvg);
	}
}


int main(int argc, char *argv[])
{
	size_t start = clock();
	//optionally, the model can be used to write the instant radiation along the full year at a given position
	//generateTirChart(-33.92584,18.42322,"Cape City");
	disinfection();
	getAllSodisData();

	float duration = (clock() - start) / (double)CLOCKS_PER_SEC;
	printf("%f\n", duration);
	return 0;
}

