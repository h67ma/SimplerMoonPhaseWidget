package sancho.simplermoonphasewidget;

import android.app.Activity;
import android.appwidget.AppWidgetManager;
import android.content.ComponentName;
import android.content.Intent;
import android.content.SharedPreferences;
import android.preference.PreferenceManager;
import android.os.Bundle;
import android.view.View;
import android.widget.CheckBox;
import android.widget.RadioGroup;

public class MoonWidgetConfig extends Activity
{
	public final static String CONFIG_MOON_AGE = "sancho.config.moonage";
	public final static String CONFIG_PHASE_NAME = "sancho.config.phasename";
	public final static String CONFIG_LOCATION_NORTHERN = "sancho.config.locationnorthern";

	private CheckBox _displayMoonAgeChk;
	private CheckBox _displayPhaseNameChk;
	private RadioGroup _locationRadioGroup;

	private int _widgetId = AppWidgetManager.INVALID_APPWIDGET_ID;

	@Override
	protected void onCreate(Bundle savedInstanceState)
	{
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_moon_widget_config);

		// get widget id
		Intent intent = getIntent();
		Bundle extras = intent.getExtras();
		if (extras != null)
			_widgetId = extras.getInt(AppWidgetManager.EXTRA_APPWIDGET_ID, AppWidgetManager.INVALID_APPWIDGET_ID);

		// id = 0 -> exit
		if (_widgetId == AppWidgetManager.INVALID_APPWIDGET_ID)	finish();

		_displayMoonAgeChk = findViewById(R.id.checkbox_display_moon_age);
		_displayPhaseNameChk = findViewById(R.id.checkbox_display_phase_name);
		_locationRadioGroup = findViewById(R.id.radiogroup_location);

		// default settings
		SharedPreferences prefs = PreferenceManager.getDefaultSharedPreferences(this);
		_displayMoonAgeChk.setChecked(prefs.getBoolean(MoonWidgetConfig.CONFIG_MOON_AGE, false));
		_displayPhaseNameChk.setChecked(prefs.getBoolean(MoonWidgetConfig.CONFIG_PHASE_NAME, false));
		if(prefs.getBoolean(MoonWidgetConfig.CONFIG_LOCATION_NORTHERN, true))
			_locationRadioGroup.check(R.id.radio_hemisphere_northern);
		else
			_locationRadioGroup.check(R.id.radio_hemisphere_southern);

		// save btn click
		findViewById(R.id.button_save).setOnClickListener(new View.OnClickListener()
		{
			@Override
			public void onClick(View v)
			{
				saveSettings();
			}
		});
	}

	private void saveSettings()
	{
		// save settings
		SharedPreferences.Editor prefs = PreferenceManager.getDefaultSharedPreferences(this).edit();
		prefs.putBoolean(CONFIG_MOON_AGE, _displayMoonAgeChk.isChecked());
		prefs.putBoolean(CONFIG_PHASE_NAME, _displayPhaseNameChk.isChecked());
		prefs.putBoolean(CONFIG_LOCATION_NORTHERN, _locationRadioGroup.getCheckedRadioButtonId() == R.id.radio_hemisphere_northern);
		prefs.apply();

		// update all widgets
		AppWidgetManager appWidgetManager = AppWidgetManager.getInstance(this);
		int[] ids = appWidgetManager.getAppWidgetIds(new ComponentName(this, MoonWidgetProvider.class));
		MoonWidgetProvider.updateWidgets(this, appWidgetManager, ids);

		// exit
		Intent resultValue = new Intent();
		resultValue.putExtra(AppWidgetManager.EXTRA_APPWIDGET_ID, _widgetId);
		setResult(RESULT_OK, resultValue);
		finish();
	}
}
