#include "JTime.h"

UTC UTC::current_lc()
{
	time_t t = time(0);
	tm * cur = localtime(&t);
	UTC total;
	total.sec = cur->tm_sec;
	total.minute = cur->tm_min;
	total.hour = cur->tm_hour;
	total.date = cur->tm_mday;
	total.month = cur->tm_mon + 1;
	total.year = cur->tm_year > 100 ? cur->tm_year - 100 : cur->tm_year;
	return total;
}
UTC UTC::current_utc()
{
	time_t t = time(0);
	tm * cur = gmtime(&t);
	UTC total;
	total.sec = cur->tm_sec;
	total.minute = cur->tm_min;
	total.hour = cur->tm_hour;
	total.date = cur->tm_mday;
	total.month = cur->tm_mon + 1;
	total.year = cur->tm_year > 100 ? cur->tm_year - 100 : cur->tm_year;
	return total;
}
int UTC::get_doy()
{
	int doy = 0;
	for (int i = 1; i < month; i++)
	{
		doy += JTime::date_amount_of_month(year, i);
	}
	doy += date;
	return doy;
}

UTC::UTC(int y, int m, int d, int h, int min, double s)
{
	year = y;
	month = m;
	date = d;
	hour = h;
	minute = min;
	sec = s;
}
void UTC::change_to_doy(int doy)
{
	int temp = doy;
	for (int i = 1; i <= 12; i++)
	{
		temp = JTime::date_amount_of_month(year, i);
		if (doy <= temp)
		{
			month = i;
			date = doy;
			return;
		}
		else
		{
			doy -= temp;
		}
	}
	throw "something wrong with the function change_to_doy";
}

bool UTC::equals(UTC * obj)
{
	if (sec == obj->sec)
		if (minute == obj->minute)
			if (hour == obj->hour)
				if (date == obj->date)
					if (month == obj->month)
						if (year == obj->year)
							return true;
	return false;
}

UTC UTC::offset_hour(int h)
{
	UTC total(year, month, date, hour, minute, sec);
	if (h > 0)
	{
		if (hour + h < 24) total.hour += h;
		else if (date < JTime::date_amount_of_month(total.year, total.month))
		{
			total.date += 1;
			total.hour = total.hour - 24 + h;
		}
		else if (month <= 11)
		{
			total.month += 1;
			total.date = 1;
			total.hour = total.hour - 24 + h;
		}
		else
		{
			total.year += 1;
			total.month = 1;
			total.date = 1;
			total.hour = total.hour - 24 + h;
		}
		return total;
	}
	else if (h < 0)
	{
		h = -h;
		if (hour >= h) total.hour -= h;
		else if (date > 1)
		{
			total.date -= 1;
			total.hour = total.hour + 24 - h;
		}
		else if (month > 1)
		{
			total.month -= 1;
			total.date = total.date - 1 + JTime::date_amount_of_month(total.year, total.month);
			total.hour = total.hour + 24 - h;
		}
		else
		{
			total.year -= 1;
			total.month += 11;
			total.date = total.date - 1 + JTime::date_amount_of_month(total.year, total.month);
			total.hour = total.hour + 24 - h;
		}
		return total;
	}
	return total;
}
int UTC::minus(UTC * obj)
{
	return (int)round(GPSTime(*this).minus(&GPSTime(*obj)));
}
UTC::UTC(MJDTime time)
{
	double rest = time.frac_day * 24;
	hour = (int)rest;
	rest = (rest - hour) * 60;
	minute = (int)rest;
	rest = (rest - minute) * 60;
	sec = (int)round(rest);

	year = (int)((time.days - 15078.2) / 365.25);
	month = (int)((time.days - 14956.1 - (int)(year * 365.25)) / 30.6001);
	date = time.days - 14956 - (int)(year * 365.25) - (int)(month * 30.6001);
	if (year > 100) year -= 100;
	month = month - 1;
}
bool UTC::larger_than(UTC u1)
{
	if (year == u1.year)
		if (month == u1.month)
			if (date == u1.date)
				if (hour == u1.hour)
					if (minute == u1.minute)
						if (sec == u1.sec)
							return true;
						else if (sec > u1.sec) return true;
						else return false;
					else if (minute > u1.minute) return true;
					else return false;
				else if (hour > u1.hour) return true;
				else return false;
			else if (date > u1.date) return true;
			else return false;
		else if (month > u1.month) return true;
		else return false;
	else if (year > u1.year) return true;
	else return false;
}

MJDTime::MJDTime(int d, double f)
{
	days = d;
	frac_day = f;
}
MJDTime::MJDTime(UTC time)
{
	int y, m, temp;
	y = time.year + (time.year < 80 ? 2000 : 1900);
	m = time.month;
	if (m <= 2)
	{
		y--;
		m += 12;
	}
	temp = (int)(365.25 * y);
	temp += (int)(30.6001 * (m + 1));
	temp += time.date;
	temp -= 679019;

	days = temp;
	frac_day = time.hour + time.minute / 60.0 + time.sec / 3600.0;
	frac_day /= 24.0;
}
MJDTime::MJDTime()
{

}
MJDTime::MJDTime(GPSTime gps)
{
	double sum = (gps.sec / 86400.0) + 44244 + 7 * gps.week;
	days = (int)sum;
	frac_day = sum - days;
}
int GPSTime::date()
{
	return sec / (24 * 3600);
}
bool GPSTime::valid()
{
	// for week
	if (week < 0 || week > 3000)return false;

	// for sec
	if (sec < 0 || sec > 7 * 24 * 3600) return false;

	return true;
}
GPSTime::GPSTime(int w, double s)
{
	week = w;
	sec = s;
}
GPSTime::GPSTime(MJDTime time)
{
	week = (int)((time.days - 44244) / 7);
	int remain = time.days - week * 7 - 44244;
	sec = (remain + time.frac_day) * 86400.0;
}
GPSTime::GPSTime(UTC time)
{
	MJDTime m_time(time);
	week = (int)((m_time.days - 44244) / 7.0);
	int remain = m_time.days - week * 7 - 44244;
	sec = (remain + m_time.frac_day) * 86400.0;
}
GPSTime::GPSTime()
{

}