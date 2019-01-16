#pragma once
#include <stdio.h>
#include <math.h>
#include <time.h>
struct MJDTime;
struct utc;
struct GPSTime;


struct JTime {
	static bool isBigMonth(int month)
	{
		int bigMonth[]{ 1,3,5,7,8,10,12 };
		for (int i = 0; i < 7; i++)
		{
			if (bigMonth[i] == month) return true;
		}
		return false;
	}
	static int date_amount_of_month(int year, int month)
	{
		if (month != 2)
		{
			return isBigMonth(month) ? 31 : 30;
		}
		else {
			if (year % 4 == 0 && year % 100 != 0 || year % 400 == 0)
			{
				return 29;
			}
			else {
				return 28;
			}
		}
	}
};

struct UTC {
	int year;
	int month;
	int date;
	int hour;
	int minute;
	double sec;

	bool equals(UTC * obj);
	int get_doy();
	static UTC current_lc();
	static UTC current_utc();
	UTC(int y, int m, int d, int h, int min, double s);
	void change_to_doy(int doy);
	UTC offset_hour(int h);
	UTC(MJDTime time);
	bool larger_than(UTC u1);
	UTC() = default;
	int minus(UTC * obj);
};
struct MJDTime
{
	int days;
	double frac_day;
	MJDTime(int d, double f);
	MJDTime(UTC time);
	MJDTime(GPSTime gps);
	MJDTime();
};

struct GPSTime
{
	int week;
	double sec;
	int date();
	bool valid();
	GPSTime(int w, double s);
	GPSTime(MJDTime time);
	GPSTime(UTC time);
	GPSTime();

	double minus(GPSTime * obj)
	{
		return (week - obj->week) * 7 * 24 * 3600 + sec - obj->sec;
	}
};