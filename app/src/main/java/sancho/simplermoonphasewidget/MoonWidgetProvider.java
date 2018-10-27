package sancho.simplermoonphasewidget;

import android.app.PendingIntent;
import android.appwidget.AppWidgetManager;
import android.appwidget.AppWidgetProvider;
import android.content.Context;
import android.content.Intent;
import android.content.SharedPreferences;
import android.preference.PreferenceManager;
import android.widget.RemoteViews;
import android.widget.Toast;

import java.text.DecimalFormat;
import java.util.Calendar;
import java.util.TimeZone;

public class MoonWidgetProvider extends AppWidgetProvider
{
	@Override
	public void onUpdate(Context context, AppWidgetManager appWidgetManager, int[] appWidgetIds)
	{
		updateWidgets(context, appWidgetManager, appWidgetIds);
	}

	private static String getPhaseName(Context context, double age)
	{
		int resId = R.string.phase_new;
		if(age < 1.05) resId = R.string.phase_new;
		else if(age < 6.33) resId = R.string.phase_waxing_crescent;
		else if(age < 8.44) resId = R.string.phase_first_quarter;
		else if(age < 13.71) resId = R.string.phase_waxing_gibbous;
		else if(age < 15.82) resId = R.string.phase_full;
		else if(age < 21.09) resId = R.string.phase_waning_gibbous;
		else if(age < 23.20) resId = R.string.phase_last_quarter;
		else if(age < 28.48) resId = R.string.phase_waning_crescent;

		return context.getResources().getString(resId);
	}

	private static int getMoonResource(double age, boolean northern)
	{
		if(!northern) age = 29.530588853 - age;
		if(age < 1.05) return R.drawable.m01; // new 1d
		else if(age < 2.11) return R.drawable.m02;
		else if(age < 3.16) return R.drawable.m03;
		else if(age < 3.16) return R.drawable.m03;
		else if(age < 4.22) return R.drawable.m04;
		else if(age < 5.27) return R.drawable.m05;
		else if(age < 6.33) return R.drawable.m06;
		else if(age < 8.44) return R.drawable.m07; // first q 2d
		else if(age < 9.49) return R.drawable.m08;
		else if(age < 10.55) return R.drawable.m09;
		else if(age < 11.60) return R.drawable.m10;
		else if(age < 12.66) return R.drawable.m11;
		else if(age < 13.71) return R.drawable.m12;
		else if(age < 15.82) return R.drawable.m13; // full 2d
		else if(age < 16.87) return R.drawable.m14;
		else if(age < 17.93) return R.drawable.m15;
		else if(age < 18.98) return R.drawable.m16;
		else if(age < 20.04) return R.drawable.m17;
		else if(age < 21.09) return R.drawable.m18;
		else if(age < 23.20) return R.drawable.m19; // last q 2d
		else if(age < 24.26) return R.drawable.m20;
		else if(age < 25.31) return R.drawable.m21;
		else if(age < 26.37) return R.drawable.m22;
		else if(age < 27.42) return R.drawable.m23;
		else if(age < 28.48) return R.drawable.m24;
		else return R.drawable.m01; // new 1d
	}

	public static void updateWidgets(Context context, AppWidgetManager appWidgetManager, int[] appWidgetIds)
	{
		for (int appWidgetId : appWidgetIds)
		{
			updateWidget(context, appWidgetManager, appWidgetId);
		}
	}

	public static void updateWidget(Context context, AppWidgetManager appWidgetManager, int appWidgetId)
	{
		double age = 0;

		Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
		int year = cal.get(Calendar.YEAR);
		int month = cal.get(Calendar.MONTH) + 1; // why is January 0????? JAVAAAAAAAAA!!!
		int day = cal.get(Calendar.DAY_OF_MONTH);
		int h = cal.get(Calendar.HOUR_OF_DAY);
		int m = cal.get(Calendar.MINUTE);
		int s = cal.get(Calendar.SECOND);

		try
		{
			SunMoonCalculator smc = new SunMoonCalculator(year, month, day, h, m, s, 0, 0);
			smc.calcSunAndMoon();
			age = smc.moonAge;
		}
		catch (Exception ex)
		{
			Toast.makeText(context, ex.getMessage(), Toast.LENGTH_LONG).show();
		}

		SharedPreferences prefs = PreferenceManager.getDefaultSharedPreferences(context);
		boolean showPhaseName = prefs.getBoolean(MoonWidgetConfig.CONFIG_PHASE_NAME, false);
		boolean showMoonAge = prefs.getBoolean(MoonWidgetConfig.CONFIG_MOON_AGE, false);
		boolean northern = prefs.getBoolean(MoonWidgetConfig.CONFIG_LOCATION_NORTHERN, true);

		String widgetText = "";

		if(showMoonAge)
		{
			widgetText += new DecimalFormat("#.#").format(age);
			if(showPhaseName) widgetText += "\n";
		}

		if(showPhaseName)
		{
			widgetText += getPhaseName(context, age);
		}

		// update layout
		RemoteViews views = new RemoteViews(context.getPackageName(), R.layout.widget_moon);

		views.setTextViewText(R.id.moon_text, widgetText);
		views.setImageViewResource(R.id.imageview_moon, getMoonResource(age, northern));

		Intent intent = new Intent(context, MoonWidgetConfig.class);
		intent.putExtra(AppWidgetManager.EXTRA_APPWIDGET_ID, appWidgetId);
		PendingIntent pendingIntent = PendingIntent.getActivity(context, 0, intent, PendingIntent.FLAG_UPDATE_CURRENT);
		views.setOnClickPendingIntent(R.id.imageview_moon, pendingIntent);

		appWidgetManager.updateAppWidget(appWidgetId, views);
	}
}
