# SimplerMoonPhaseWidget
Simple widget that displays current moon phase.

Another one of my ["there are lots of apps that do this thing but they have far too many features and using them is a pain if I just want a single feature"](https://github.com/szycikm/LOSSimpleProfileSwitcher) apps. This one contains a single widget that displays current moon phase and age. There are hundreds of apps that provide this functionality, I know, but I couldn't find a single one that didn't have some kind of calendar, skymap, or other celestial bodies built in. Also all of them displayed pretty photos of moon, and I wanted just a simple schematic representation (a.k.a. material design). Also ads.

I guess another reason I made it is just for practise and fun.

So I present to you: an app with one widget and one configuration activity. Widget displays a simple vector graphic (all graphics made by me). Configuration activity lets you enable optional displaying of phase name and moon age, and setting location (Northern/Southern Hemisphere).

The whole thing is a whopping 124KB after install.

**I don't pretend to understand the algorithm for calculating moon age, it's taken from [here](http://conga.oan.es/~alonso/doku.php?id=blog:sun_moon_position#demostration_program_for_java_desktop).** Turns out that calculating moon age is not a simple task. Current solution calculates new/full/quarter moments perfectly, but somewhat "lags" in between those keypoints (and I don't have a clue why). It's still the best algorithm I could find (I also tried "simple", Conway, and some code from 1985, but all are inaccurate). After spending some time trying to fix this issue, I realized that for something as silly as a phone widget, current solution is already sufficient.
