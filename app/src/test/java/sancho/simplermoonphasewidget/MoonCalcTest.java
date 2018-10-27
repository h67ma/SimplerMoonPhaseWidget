package sancho.simplermoonphasewidget;

import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class MoonCalcTest
{
	private static final double ERROR_MARGIN = 0.1;
	private static final double MOON_CYCLE_LENGTH_DAYS = 29.530588853;

	// from http://www.giangrandi.ch/soft/mooncalc/mooncalc.shtml

	@Test
	public void NewTest()
	{
		double age = 0;
		try
		{
			SunMoonCalculator smc = new SunMoonCalculator(2018, 10, 9, 5, 47, 0, 0, 0);
			smc.calcSunAndMoon();
			age = smc.moonAge;
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}

		assertTrue("Age = " + String.valueOf(age) + ", expected " + 0,
				age < ERROR_MARGIN || Math.abs(MOON_CYCLE_LENGTH_DAYS - age) < ERROR_MARGIN);
	}

	@Test
	public void NewTest2()
	{
		double age = 0;
		try
		{
			SunMoonCalculator smc = new SunMoonCalculator(2018, 11, 7, 17, 2, 0, 0, 0);
			smc.calcSunAndMoon();
			age = smc.moonAge;
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}

		assertTrue("Age = " + String.valueOf(age) + ", expected " + 0,
				age < ERROR_MARGIN || Math.abs(MOON_CYCLE_LENGTH_DAYS - age) < ERROR_MARGIN);
	}

	@Test
	public void FirstQuarterTest()
	{
		double age = 0;
		try
		{
			SunMoonCalculator smc = new SunMoonCalculator(2018, 11, 15, 15, 34, 0, 0, 0);
			smc.calcSunAndMoon();
			age = smc.moonAge;
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}

		assertTrue("Age = " + String.valueOf(age) + ", expected " + MOON_CYCLE_LENGTH_DAYS/4,
				Math.abs(MOON_CYCLE_LENGTH_DAYS/4 - age) < ERROR_MARGIN);
	}

	@Test
	public void FullTest()
	{
		double age = 0;
		try
		{
			SunMoonCalculator smc = new SunMoonCalculator(2018, 11, 23, 6, 39, 0, 0, 0);
			smc.calcSunAndMoon();
			age = smc.moonAge;
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}

		assertTrue("Age = " + String.valueOf(age) + ", expected " + MOON_CYCLE_LENGTH_DAYS/2,
				Math.abs(MOON_CYCLE_LENGTH_DAYS/2 - age) < ERROR_MARGIN);
	}

	@Test
	public void LastQuarterTest()
	{
		double age = 0;
		try
		{
			SunMoonCalculator smc = new SunMoonCalculator(2018, 10, 31, 17, 40, 0, 0, 0);
			smc.calcSunAndMoon();
			age = smc.moonAge;
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}

		assertTrue("Age = " + String.valueOf(age) + ", expected " + MOON_CYCLE_LENGTH_DAYS*3/4,
				Math.abs(MOON_CYCLE_LENGTH_DAYS*3/4 - age) < ERROR_MARGIN);
	}

	// from http://cycletourist.com/moon/
	// age calculated from multiplying moon cycle by percent: e.g. 29.530588853 * 0.615 = 18.161312144595 (matches SunMoonCalculator's output)
	// differs from the one shown in days: 18 days 18 hours 31 minutes = 18.771527777777777777777777777778 by few percent
	// what gives??

	/*
	@Test
	public void RandomTest1()
	{
		double age = 0;
		try
		{
			SunMoonCalculator smc = new SunMoonCalculator(2018, 10, 28, 0, 18, 0, 0, 0);
			smc.calcSunAndMoon();
			age = smc.moonAge;
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}

		assertTrue("Age = " + String.valueOf(age) + ", expected " + 18.161312144595,
				Math.abs(18.161312144595 - age) < ERROR_MARGIN);
	}

	@Test
	public void RandomTest2()
	{
		double age = 0;
		try
		{
			SunMoonCalculator smc = new SunMoonCalculator(2018, 2, 13, 3, 9, 0, 0, 0);
			smc.calcSunAndMoon();
			age = smc.moonAge;
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}

		assertTrue("Age = " + String.valueOf(age) + ", expected " + 28.831,
				Math.abs(28.831 - age) < ERROR_MARGIN);
	}
	*/
}
