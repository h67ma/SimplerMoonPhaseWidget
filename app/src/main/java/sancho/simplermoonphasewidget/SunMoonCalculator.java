package sancho.simplermoonphasewidget;

/**
 * A very simple Sun/Moon calculator without using JPARSEC library.
 * @author T. Alonso Albi - OAN (Spain), email t.alonso@oan.es
 * @version July 24, 2018 (new class to hold results, illumination phase, moon phases, equinoxes and solstices)
 * @version May 25, 2017 (fixed nutation correction and moon age, better accuracy in Moon)
 *
 * source: http://conga.oan.es/~alonso/doku.php?id=blog:sun_moon_position#demostration_program_for_java_desktop
 * modified code - removed few unused methods, refactored code a bit
 */
class SunMoonCalculator
{
	/** Radians to degrees. */
	private static final double RAD_TO_DEG = 180.0 / Math.PI;

	/** Degrees to radians. */
	static final double DEG_TO_RAD = 1.0 / RAD_TO_DEG;

	/** Radians to hours. */
	private static final double RAD_TO_HOUR = 180.0 / (15.0 * Math.PI);

	/** Radians to days. */
	private static final double RAD_TO_DAY = RAD_TO_HOUR / 24.0;

	/** Astronomical Unit in km. As defined by JPL. */
	private static final double AU = 149597870.691;

	/** Earth equatorial radius in km. IERS 2003 Conventions. */
	private static final double EARTH_RADIUS = 6378.1366;

	/** Two times Pi. */
	private static final double TWO_PI = 2.0 * Math.PI;

	/** The inverse of two times Pi. */
	private static final double TWO_PI_INVERSE = 1.0 / (2.0 * Math.PI);

	/** Four times Pi. */
	private static final double FOUR_PI = 4.0 * Math.PI;

	/** Pi divided by two. */
	private static final double PI_OVER_TWO = Math.PI / 2.0;

	/** Length of a sidereal day in days according to IERS Conventions. */
	private static final double SIDEREAL_DAY_LENGTH = 1.00273781191135448;

	/** Julian century conversion constant = 100 * days per year. */
	private static final double JULIAN_DAYS_PER_CENTURY = 36525.0;

	/** Seconds in one day. */
	private static final double SECONDS_PER_DAY = 86400;

	/** Our default epoch. <BR>
	 * The Julian Day which represents noon on 2000-01-01. */
	private static final double J2000 = 2451545.0;

	/**
	 * The set of twilights to calculate (types of rise/set events).
	 */
	public enum TWILIGHT {
		/**
		 * Event ID for calculation of rising and setting times for astronomical
		 * twilight. In this case, the calculated time will be the time when the
		 * center of the object is at -18 degrees of geometrical elevation below the
		 * astronomical horizon. At this time astronomical observations are possible
		 * because the sky is dark enough.
		 */
		TWILIGHT_ASTRONOMICAL,
		/**
		 * Event ID for calculation of rising and setting times for nautical
		 * twilight. In this case, the calculated time will be the time when the
		 * center of the object is at -12 degrees of geometric elevation below the
		 * astronomical horizon.
		 */
		TWILIGHT_NAUTICAL,
		/**
		 * Event ID for calculation of rising and setting times for civil twilight.
		 * In this case, the calculated time will be the time when the center of the
		 * object is at -6 degrees of geometric elevation below the astronomical
		 * horizon.
		 */
		TWILIGHT_CIVIL,
		/**
		 * The standard value of 34' for the refraction at the local horizon.
		 */
		HORIZON_34arcmin
	};

	/**
	 * The set of events to calculate (rise/set/transit events).
	 */
	public enum EVENT {
		/** Rise. */
		RISE,
		/** Set. */
		SET,
		/** Transit. */
		TRANSIT
	}

	/** Input values. */
	private double jd_UT = 0, t = 0, obsLon = 0, obsLat = 0, TTminusUT = 0;
	private TWILIGHT twilight = TWILIGHT.HORIZON_34arcmin;

	/**
	 * Class to hold the results of ephemerides.
	 * @author T. Alonso Albi - OAN (Spain)
	 */
	public class Ephemeris {
		private Ephemeris(double azi, double alt, double rise2, double set2,
		                  double transit2, double transit_alt, double ra, double dec,
		                  double dist, double eclLon, double eclLat) {
			azimuth = azi;
			elevation = alt;
			rise = rise2;
			set = set2;
			transit = transit2;
			transitElevation = transit_alt;
			rightAscension = ra;
			declination = dec;
			distance = dist;
			illuminationPhase = 100;
			eclipticLongitude = eclLon;
			eclipticLatitude = eclLat;
		}

		/** Values for azimuth, elevation, rise, set, and transit for the Sun. Angles in radians, rise ...
		 * as Julian days in UT. Distance in AU. */
		double azimuth, elevation, rise, set, transit, transitElevation, distance, rightAscension,
				declination, illuminationPhase, eclipticLongitude, eclipticLatitude;
	}

	/** Ephemeris for the Sun and Moon bodies. */
	private Ephemeris sun;

	/** Moon's age in days as an independent variable. */
	double moonAge;

	/**
	 * Main constructor for Sun/Moon calculations. Time should be given in
	 * Universal Time (UT), observer angles in radians.
	 * @param year The year.
	 * @param month The month.
	 * @param day The day.
	 * @param h The hour.
	 * @param m Minute.
	 * @param s Second.
	 * @param obsLon Longitude for the observer.
	 * @param obsLat Latitude for the observer.
	 * @throws Exception If the date does not exist.
	 */
	SunMoonCalculator(int year, int month, int day, int h, int m, int s, double obsLon, double obsLat) throws Exception {
		double jd = toJulianDay(year, month, day, h, m, s);

		this.TTminusUT = 0;
		if (year > -600 && year < 2200) {
			double x = year + (month - 1 + day / 30.0) / 12.0;
			double x2 = x * x, x3 = x2 * x, x4 = x3 * x;
			if (year < 1600) {
				this.TTminusUT = 10535.328003326353 - 9.995238627481024 * x + 0.003067307630020489 * x2 - 7.76340698361363E-6 * x3 + 3.1331045394223196E-9 * x4 +
						8.225530854405553E-12 * x2 * x3 - 7.486164715632051E-15 * x4 * x2 + 1.9362461549678834E-18 * x4 * x3 - 8.489224937827653E-23 * x4 * x4;
			} else {
				this.TTminusUT = -1027175.3477559977 + 2523.256625418965 * x - 1.885686849058459 * x2 + 5.869246227888417E-5 * x3 + 3.3379295816475025E-7 * x4 +
						1.7758961671447929E-10 * x2 * x3 - 2.7889902806153024E-13 * x2 * x4 + 1.0224295822336825E-16 * x3 * x4 - 1.2528102370680435E-20 * x4 * x4;
			}
		}
		this.obsLon = obsLon;
		this.obsLat = obsLat;
		setUTDate(jd);
	}

	private double toJulianDay(int year, int month, int day, int h, int m, int s) throws Exception {
		// The conversion formulas are from Meeus, chapter 7.
		boolean julian = false; // Use Gregorian calendar
		if (year < 1582 || (year == 1582 && month <= 10) || (year == 1582 && month == 10 && day < 15)) julian = true;
		int D = day;
		int M = month;
		int Y = year;
		if (M < 3)
		{
			Y--;
			M += 12;
		}
		int A = Y / 100;
		int B = julian ? 0 : 2 - A + A / 4;

		double dayFraction = (h + (m + (s / 60.0)) / 60.0) / 24.0;

		double jd = dayFraction + (int) (365.25D * (Y + 4716)) + (int) (30.6001 * (M + 1)) + D + B - 1524.5;

		if (jd < 2299160.0 && jd >= 2299150.0)
		{
			throw new Exception("invalid julian day " + jd + ". This date does not exist.");
		}
		return jd;
	}

	private void setUTDate(double jd) {
		this.jd_UT = jd;
		this.t = (jd + this.TTminusUT / SECONDS_PER_DAY  - J2000) / JULIAN_DAYS_PER_CENTURY;
	}

	/** Calculates everything for the Sun and the Moon. */
	void calcSunAndMoon() {
		double jd = this.jd_UT;

		// First the Sun
		sun = doCalc(getSun(), false);

		int niter = 3; // Number of iterations to get accurate rise/set/transit times
		sun.rise = obtainAccurateRiseSetTransit(sun.rise, EVENT.RISE, niter, true);
		sun.set = obtainAccurateRiseSetTransit(sun.set, EVENT.SET, niter, true);
		sun.transit = obtainAccurateRiseSetTransit(sun.transit, EVENT.TRANSIT, niter, true);
		if (sun.transit == -1) {
			sun.transitElevation = 0;
		} else {
			// Update Sun's maximum elevation
			setUTDate(sun.transit);
			sun.transitElevation = doCalc(getSun(), false).transitElevation;
		}

		// Now Moon
		setUTDate(jd);
		Ephemeris moon = doCalc(getMoon(), false);
		double ma = moonAge;

		niter = 5; // Number of iterations to get accurate rise/set/transit times
		moon.rise = obtainAccurateRiseSetTransit(moon.rise, EVENT.RISE, niter, false);
		moon.set = obtainAccurateRiseSetTransit(moon.set, EVENT.SET, niter, false);
		moon.transit = obtainAccurateRiseSetTransit(moon.transit, EVENT.TRANSIT, niter, false);
		if (moon.transit == -1) {
			moon.transitElevation = 0;
		} else {
			// Update Moon's maximum elevation
			setUTDate(moon.transit);
			getSun();
			moon.transitElevation = doCalc(getMoon(), false).transitElevation;
		}
		setUTDate(jd);
		moonAge = ma;

		// Compute illumination phase percentage for the Moon (do not use for other bodies!)
		double dlon = moon.rightAscension - sun.rightAscension;
		double elong = Math.acos(Math.sin(sun.declination) * Math.sin(moon.declination) +
				Math.cos(sun.declination) * Math.cos(moon.declination) * Math.cos(dlon));
		moon.illuminationPhase = 100 * (1.0 - Math.cos(elong)) * 0.5;
	}

	// Formulae here is a simplifcation of the expansion from
	// "Planetary Programs and Tables" by Pierre Bretagnon and
	// Jean-Louis Simon, Willman-Bell, 1986. This source also
	// have expansions for ephemerides of planets
	private static double sun_elements[][] = {
			new double[] { 403406.0, 0.0, 4.721964, 1.621043 },
			new double[] { 195207.0, -97597.0, 5.937458, 62830.348067 },
			new double[] { 119433.0, -59715.0, 1.115589, 62830.821524 },
			new double[] { 112392.0, -56188.0, 5.781616, 62829.634302 },
			new double[] { 3891.0, -1556.0, 5.5474, 125660.5691 },
			new double[] { 2819.0, -1126.0, 1.512, 125660.9845 },
			new double[] { 1721.0, -861.0, 4.1897, 62832.4766 },
			new double[] { 0.0, 941.0, 1.163, .813 },
			new double[] { 660.0, -264.0, 5.415, 125659.31 },
			new double[] { 350.0, -163.0, 4.315, 57533.85 },
			new double[] { 334.0, 0.0, 4.553, -33.931 },
			new double[] { 314.0, 309.0, 5.198, 777137.715 },
			new double[] { 268.0, -158.0, 5.989, 78604.191 },
			new double[] { 242.0, 0.0, 2.911, 5.412 },
			new double[] { 234.0, -54.0, 1.423, 39302.098 },
			new double[] { 158.0, 0.0, .061, -34.861 },
			new double[] { 132.0, -93.0, 2.317, 115067.698 },
			new double[] { 129.0, -20.0, 3.193, 15774.337 },
			new double[] { 114.0, 0.0, 2.828, 5296.67 },
			new double[] { 99.0, -47.0, .52, 58849.27 },
			new double[] { 93.0, 0.0, 4.65, 5296.11 },
			new double[] { 86.0, 0.0, 4.35, -3980.7 },
			new double[] { 78.0, -33.0, 2.75, 52237.69 },
			new double[] { 72.0, -32.0, 4.5, 55076.47 },
			new double[] { 68.0, 0.0, 3.23, 261.08 },
			new double[] { 64.0, -10.0, 1.22, 15773.85 }
	};

	private double[] getSun() {
		double L = 0.0, R = 0.0;
		double t2 = t * 0.01;
		for (double[] sun_element : sun_elements)
		{
			double v = sun_element[2] + sun_element[3] * t2;
			double u = normalizeRadians(v);
			L = L + sun_element[0] * Math.sin(u);
			R = R + sun_element[1] * Math.cos(u);
		}

		double lon = normalizeRadians(4.9353929 + normalizeRadians(62833.196168 * t2) + L / 10000000.0) * RAD_TO_DEG;
		double sdistance = 1.0001026 + R / 10000000.0;

		// Now, let calculate nutation
		double M1 = (124.90 - 1934.134 * t + 0.002063 * t * t) * DEG_TO_RAD;
		double M2 = (201.11 + 72001.5377 * t + 0.00057 * t * t) * DEG_TO_RAD;
		double d = -.00569 - .0047785 * Math.sin(M1) - .0003667 * Math.sin(M2);

		double slongitude = lon + d; // apparent longitude (error<0.001 deg)
		double slatitude = 0; // Sun's ecliptic latitude is always negligible

		return new double[] {slongitude, slatitude, sdistance, Math.atan(696000 / (AU * sdistance))};
	}

	private double[] getMoon() {
		// MOON PARAMETERS (Formulae from "Calendrical Calculations")
		double phase = normalizeRadians((297.8502042 + 445267.1115168 * t - 0.00163 * t * t + t * t * t / 538841 - t * t * t * t / 65194000) * DEG_TO_RAD);

		// Anomalistic phase
		double anomaly = (134.9634114 + 477198.8676313 * t + .008997 * t * t + t * t * t / 69699 - t * t * t * t / 14712000);
		anomaly = anomaly * DEG_TO_RAD;

		// Degrees from ascending node
		double node = (93.2720993 + 483202.0175273 * t - 0.0034029 * t * t - t * t * t / 3526000 + t * t * t * t / 863310000);
		node = node * DEG_TO_RAD;

		double E = 1.0 - (.002495 + 7.52E-06 * (t + 1.0)) * (t + 1.0);

		// Solar anomaly
		double sanomaly = (357.5291 + 35999.0503 * t - .0001559 * t * t - 4.8E-07 * t * t * t) * DEG_TO_RAD;

		// Now longitude, with the three main correcting terms of evection,
		// variation, and equation of year, plus other terms (error<0.01 deg)
		// P. Duffet's MOON program taken as reference
		double l = (218.31664563 + 481267.8811958 * t - .00146639 * t * t + t * t * t / 540135.03 - t * t * t * t / 65193770.4);
		l += 6.28875 * Math.sin(anomaly) + 1.274018 * Math.sin(2 * phase - anomaly) + .658309 * Math.sin(2 * phase);
		l +=  0.213616 * Math.sin(2 * anomaly) - E * .185596 * Math.sin(sanomaly) - 0.114336 * Math.sin(2 * node);
		l += .058793 * Math.sin(2 * phase - 2 * anomaly) + .057212 * E * Math.sin(2 * phase - anomaly - sanomaly) + .05332 * Math.sin(2 * phase + anomaly);
		l += .045874 * E * Math.sin(2 * phase - sanomaly) + .041024 * E * Math.sin(anomaly - sanomaly) - .034718 * Math.sin(phase) - E * .030465 * Math.sin(sanomaly + anomaly);
		l += .015326 * Math.sin(2 * (phase - node)) - .012528 * Math.sin(2 * node + anomaly) - .01098 * Math.sin(2 * node - anomaly) + .010674 * Math.sin(4 * phase - anomaly);
		l += .010034 * Math.sin(3 * anomaly) + .008548 * Math.sin(4 * phase - 2 * anomaly);
		l += -E * .00791 * Math.sin(sanomaly - anomaly + 2 * phase) - E * .006783 * Math.sin(2 * phase + sanomaly) + .005162 * Math.sin(anomaly - phase) + E * .005 * Math.sin(sanomaly + phase);
		l += .003862 * Math.sin(4 * phase) + E * .004049 * Math.sin(anomaly - sanomaly + 2 * phase) + .003996 * Math.sin(2 * (anomaly + phase)) + .003665 * Math.sin(2 * phase - 3 * anomaly);
		l += E * 2.695E-3 * Math.sin(2 * anomaly - sanomaly) + 2.602E-3 * Math.sin(anomaly - 2*(node+phase));
		l += E * 2.396E-3 * Math.sin(2*(phase - anomaly) - sanomaly) - 2.349E-3 * Math.sin(anomaly+phase);
		l += E * E * 2.249E-3 * Math.sin(2*(phase-sanomaly)) - E * 2.125E-3 * Math.sin(2*anomaly+sanomaly);
		l += -E * E * 2.079E-3 * Math.sin(2*sanomaly) + E * E * 2.059E-3 * Math.sin(2*(phase-sanomaly)-anomaly);
		l += -1.773E-3 * Math.sin(anomaly+2*(phase-node)) - 1.595E-3 * Math.sin(2*(node+phase));
		l += E * 1.22E-3 * Math.sin(4*phase-sanomaly-anomaly) - 1.11E-3 * Math.sin(2*(anomaly+node));
		double longitude = l;

		// Let's add nutation here also
		double M1 = (124.90 - 1934.134 * t + 0.002063 * t * t) * DEG_TO_RAD;
		double M2 = (201.11 + 72001.5377 * t + 0.00057 * t * t) * DEG_TO_RAD;
		double d = - .0047785 * Math.sin(M1) - .0003667 * Math.sin(M2);
		longitude += d;

		// Get accurate Moon age
		double Psin = 29.530588853;
		moonAge = normalizeRadians(longitude * DEG_TO_RAD - sun.eclipticLongitude) * Psin / TWO_PI;

		// Now Moon parallax
		double parallax = .950724 + .051818 * Math.cos(anomaly) + .009531 * Math.cos(2 * phase - anomaly);
		parallax += .007843 * Math.cos(2 * phase) + .002824 * Math.cos(2 * anomaly);
		parallax += 0.000857 * Math.cos(2 * phase + anomaly) + E * .000533 * Math.cos(2 * phase - sanomaly);
		parallax += E * .000401 * Math.cos(2 * phase - anomaly - sanomaly) + E * .00032 * Math.cos(anomaly - sanomaly) - .000271 * Math.cos(phase);
		parallax += -E * .000264 * Math.cos(sanomaly + anomaly) - .000198 * Math.cos(2 * node - anomaly);
		parallax += 1.73E-4 * Math.cos(3 * anomaly) + 1.67E-4 * Math.cos(4*phase-anomaly);

		// So Moon distance in Earth radii is, more or less,
		double distance = 1.0 / Math.sin(parallax * DEG_TO_RAD);

		// Ecliptic latitude with nodal phase (error<0.01 deg)
		l = 5.128189 * Math.sin(node) + 0.280606 * Math.sin(node + anomaly) + 0.277693 * Math.sin(anomaly - node);
		l += .173238 * Math.sin(2 * phase - node) + .055413 * Math.sin(2 * phase + node - anomaly);
		l += .046272 * Math.sin(2 * phase - node - anomaly) + .032573 * Math.sin(2 * phase + node);
		l += .017198 * Math.sin(2 * anomaly + node) + .009267 * Math.sin(2 * phase + anomaly - node);
		l += .008823 * Math.sin(2 * anomaly - node) + E * .008247 * Math.sin(2 * phase - sanomaly - node) + .004323 * Math.sin(2 * (phase - anomaly) - node);
		l += .0042 * Math.sin(2 * phase + node + anomaly) + E * .003372 * Math.sin(node - sanomaly - 2 * phase);
		l += E * 2.472E-3 * Math.sin(2 * phase + node - sanomaly - anomaly);
		l += E * 2.222E-3 * Math.sin(2 * phase + node - sanomaly);
		l += E * 2.072E-3 * Math.sin(2 * phase - node - sanomaly - anomaly);
		double latitude = l;

		return new double[] {longitude, latitude, distance * EARTH_RADIUS / AU, Math.atan(1737.4 / (distance * EARTH_RADIUS))};
	}

	private Ephemeris doCalc(double[] pos, boolean geocentric) {
		// Ecliptic to equatorial coordinates
		double t2 = this.t / 100.0;
		double tmp = t2 * (27.87 + t2 * (5.79 + t2 * 2.45));
		tmp = t2 * (-249.67 + t2 * (-39.05 + t2 * (7.12 + tmp)));
		tmp = t2 * (-1.55 + t2 * (1999.25 + t2 * (-51.38 + tmp)));
		tmp = (t2 * (-4680.93 + tmp)) / 3600.0;
		double angle = (23.4392911111111 + tmp) * DEG_TO_RAD; // obliquity

		// Add nutation in obliquity
		double M1 = (124.90 - 1934.134 * t + 0.002063 * t * t) * DEG_TO_RAD;
		double M2 = (201.11 + 72001.5377 * t + 0.00057 * t * t) * DEG_TO_RAD;
		double d = .002558 * Math.cos(M1) - .00015339 * Math.cos(M2);
		angle += d * DEG_TO_RAD;

		pos[0] *= DEG_TO_RAD;
		pos[1] *= DEG_TO_RAD;
		double cl = Math.cos(pos[1]);
		double x = pos[2] * Math.cos(pos[0]) * cl;
		double y = pos[2] * Math.sin(pos[0]) * cl;
		double z = pos[2] * Math.sin(pos[1]);
		tmp = y * Math.cos(angle) - z * Math.sin(angle);
		z = y * Math.sin(angle) + z * Math.cos(angle);
		y = tmp;

		if (geocentric) return new Ephemeris(0, 0, -1, -1, -1, -1, normalizeRadians(Math.atan2(y, x)),
				Math.atan2(z / Math.sqrt(x * x + y * y), 1.0), 	Math.sqrt(x * x + y * y + z * z), pos[0], pos[1]);

		// Obtain local apparent sidereal time
		double jd0 = Math.floor(jd_UT - 0.5) + 0.5;
		double T0 = (jd0 - J2000) / JULIAN_DAYS_PER_CENTURY;
		double secs = (jd_UT - jd0) * SECONDS_PER_DAY;
		double gmst = (((((-6.2e-6 * T0) + 9.3104e-2) * T0) + 8640184.812866) * T0) + 24110.54841;
		double msday = 1.0 + (((((-1.86e-5 * T0) + 0.186208) * T0) + 8640184.812866) / (SECONDS_PER_DAY * JULIAN_DAYS_PER_CENTURY));
		gmst = (gmst + msday * secs) * (15.0 / 3600.0) * DEG_TO_RAD;
		double lst = gmst + obsLon;

		// Obtain topocentric rectangular coordinates
		double radiusAU = EARTH_RADIUS / AU;
		double correction[] = new double[] {
				radiusAU * Math.cos(obsLat) * Math.cos(lst),
				radiusAU * Math.cos(obsLat) * Math.sin(lst),
				radiusAU * Math.sin(obsLat)};
		double xtopo = x - correction[0];
		double ytopo = y - correction[1];
		double ztopo = z - correction[2];

		// Obtain topocentric equatorial coordinates
		double ra = 0.0;
		double dec = PI_OVER_TWO;
		if (ztopo < 0.0)
			dec = -dec;
		if (ytopo != 0.0 || xtopo != 0.0)
		{
			ra = Math.atan2(ytopo, xtopo);
			dec = Math.atan2(ztopo / Math.sqrt(xtopo * xtopo + ytopo * ytopo), 1.0);
		}
		double dist = Math.sqrt(xtopo * xtopo + ytopo * ytopo + ztopo * ztopo);

		// Hour angle
		double angh = lst - ra;

		// Obtain azimuth and geometric alt
		double sinlat = Math.sin(obsLat);
		double coslat = Math.cos(obsLat);
		double sindec = Math.sin(dec), cosdec = Math.cos(dec);
		double h = sinlat * sindec + coslat * cosdec * Math.cos(angh);
		double alt = Math.asin(h);
		double azy = Math.sin(angh);
		double azx = Math.cos(angh) * sinlat - sindec * coslat / cosdec;
		double azi = Math.PI + Math.atan2(azy, azx); // 0 = north

		// Get apparent elevation
		if (alt > -3 * DEG_TO_RAD) {
			double r = 0.016667 * DEG_TO_RAD * Math.abs(Math.tan(PI_OVER_TWO - (alt * RAD_TO_DEG +  7.31 / (alt * RAD_TO_DEG + 4.4)) * DEG_TO_RAD));
			double refr = r * ( 0.28 * 1010 / (10 + 273.0)); // Assuming pressure of 1010 mb and T = 10 C
			alt = Math.min(alt + refr, PI_OVER_TWO); // This is not accurate, but acceptable
		}

		switch (twilight) {
			case HORIZON_34arcmin:
				// Rise, set, transit times, taking into account Sun/Moon angular radius (pos[3]).
				// The 34' factor is the standard refraction at horizon.
				// Removing angular radius will do calculations for the center of the disk instead
				// of the upper limb.
				tmp = -(34.0 / 60.0) * DEG_TO_RAD - pos[3];
				break;
			case TWILIGHT_CIVIL:
				tmp = -6 * DEG_TO_RAD;
				break;
			case TWILIGHT_NAUTICAL:
				tmp = -12 * DEG_TO_RAD;
				break;
			case TWILIGHT_ASTRONOMICAL:
				tmp = -18 * DEG_TO_RAD;
				break;
		}

		// Compute cosine of hour angle
		tmp = (Math.sin(tmp) - Math.sin(obsLat) * Math.sin(dec)) / (Math.cos(obsLat) * Math.cos(dec));
		double celestialHoursToEarthTime = RAD_TO_DAY / SIDEREAL_DAY_LENGTH;

		// Make calculations for the meridian
		double transit_time1 = celestialHoursToEarthTime * normalizeRadians(ra - lst);
		double transit_time2 = celestialHoursToEarthTime * (normalizeRadians(ra - lst) - TWO_PI);
		double transit_alt = Math.asin(Math.sin(dec) * Math.sin(obsLat) + Math.cos(dec) * Math.cos(obsLat));
		if (transit_alt > -3 * DEG_TO_RAD) {
			double r = 0.016667 * DEG_TO_RAD * Math.abs(Math.tan(PI_OVER_TWO - (transit_alt * RAD_TO_DEG +  7.31 / (transit_alt * RAD_TO_DEG + 4.4)) * DEG_TO_RAD));
			double refr = r * ( 0.28 * 1010 / (10 + 273.0)); // Assuming pressure of 1010 mb and T = 10 C
			transit_alt = Math.min(transit_alt + refr, PI_OVER_TWO); // This is not accurate, but acceptable
		}

		// Obtain the current event in time
		double transit_time = transit_time1;
		double jdToday = Math.floor(jd_UT - 0.5) + 0.5;
		double transitToday2 = Math.floor(jd_UT + transit_time2 - 0.5) + 0.5;
		// Obtain the transit time. Preference should be given to the closest event
		// in time to the current calculation time
		if (jdToday == transitToday2 && Math.abs(transit_time2) < Math.abs(transit_time1)) transit_time = transit_time2;
		double transit = jd_UT + transit_time;

		// Make calculations for rise and set
		double rise = -1, set = -1;
		if (Math.abs(tmp) <= 1.0)
		{
			double ang_hor = Math.abs(Math.acos(tmp));
			double rise_time1 = celestialHoursToEarthTime * normalizeRadians(ra - ang_hor - lst);
			double set_time1 = celestialHoursToEarthTime * normalizeRadians(ra + ang_hor - lst);
			double rise_time2 = celestialHoursToEarthTime * (normalizeRadians(ra - ang_hor - lst) - TWO_PI);
			double set_time2 = celestialHoursToEarthTime * (normalizeRadians(ra + ang_hor - lst) - TWO_PI);

			// Obtain the current events in time. Preference should be given to the closest event
			// in time to the current calculation time (so that iteration in other method will converge)
			double rise_time = rise_time1;
			double riseToday2 = Math.floor(jd_UT + rise_time2 - 0.5) + 0.5;
			if (jdToday == riseToday2 && Math.abs(rise_time2) < Math.abs(rise_time1)) rise_time = rise_time2;

			double set_time = set_time1;
			double setToday2 = Math.floor(jd_UT + set_time2 - 0.5) + 0.5;
			if (jdToday == setToday2 && Math.abs(set_time2) < Math.abs(set_time1)) set_time = set_time2;
			rise = jd_UT + rise_time;
			set = jd_UT + set_time;
		}

		return new Ephemeris(azi, alt, rise, set, transit, transit_alt,
				normalizeRadians(ra), dec, dist, pos[0], pos[1]);
	}

	/**
	 * Reduce an angle in radians to the range (0 - 2 Pi).
	 *
	 * @param r Value in radians.
	 * @return The reduced radians value.
	 */
	private static double normalizeRadians(double r)
	{
		if (r < 0 && r >= -TWO_PI) return r + TWO_PI;
		if (r >= TWO_PI && r < FOUR_PI) return r - TWO_PI;
		if (r >= 0 && r < TWO_PI) return r;

		r -= TWO_PI * Math.floor(r * TWO_PI_INVERSE);
		if (r < 0.) r += TWO_PI;

		return r;
	}

	private double obtainAccurateRiseSetTransit(double riseSetJD, EVENT index, int niter, boolean sun) {
		double step = -1;
		for (int i = 0; i< niter; i++) {
			if (riseSetJD == -1) return riseSetJD; // -1 means no rise/set from that location
			setUTDate(riseSetJD);
			Ephemeris out;
			if (sun) {
				out = doCalc(getSun(), false);
			} else {
				getSun();
				out = doCalc(getMoon(), false);
			}

			double val = out.rise;
			if (index == EVENT.SET) val = out.set;
			if (index == EVENT.TRANSIT) val = out.transit;
			step = Math.abs(riseSetJD - val);
			riseSetJD = val;
		}
		if (step > 1.0 / SECONDS_PER_DAY) return -1; // did not converge => without rise/set/transit in this date
		return riseSetJD;
	}
}
