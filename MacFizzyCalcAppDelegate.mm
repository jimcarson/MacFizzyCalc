/****************************************************************************
 FizzyCalcAppDelegate.m
 MacFizzyCalc
 
 Copyright (c) 2011 by James Carson (www.jimcarson.com).  
 Uses geocaching coordinate calculations Copyright (c) 2004 by David Knapp.
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

/* To do:
 * 08/30/2011:
 * Add preferences pane to change precision & suppress comma
 * After coordinates are parsed, set the precision and comma.
 * Add parsing where comma is not present (doable?  insane??)
 * 
 *   Support undo in text fields via %Z.
 *   Consider making a subroutine out of the integer versus string defaults values.
 *
 * Bugs fixed:
 *  8/28 If no distance value is entered in projection, you end up with NAN.
 *  8/28 If a user entered a bogus value in the combobox for the WAAS satellite, it would crash/hang.  
 Problem was NSComboBox was returning -1, which caused a bad reference into the struct
 definining satellites.  8/28/2011 - fixed
 *  8/28 Change the default distance value so "meters" makes sense.
 *  8/25 WAAS was using P_start instead of W_Coordinate.
 */

#import "MacFizzyCalcAppDelegate.h"
#import "Coord.h"


// set up the default units and satellite
// For most tabs
#define COORD_KEY               @"COORDINATE"
#define DEFAULT_COORDS_VALUE    @"N 47° 37.777 W 122° 06.310" // My First iPhone Hide (GC2XZKZ) :-)
#define COORD_KEY2              @"COORDINATE2"
#define DEFAULT_COORDS_VALUE2   @"N 47° 37.048 W 122° 07.907"
#define UNITS_KEY               @"UNITS"
#define DEFAULT_UNITS_VALUE     2 // meters
// For the Projection tab
#define DIST_KEY                @"DISTANCE"
#define DEFAULT_DIST_VALUE      @"75.5"
#define BEARING_KEY             @"BEARING"
#define DEFAULT_BEARING_VALUE   @"56"
// For the WAAS tab:
#define SAT_KEY                 @"SATELLITE"
#define DEFAULT_SAT_VALUE       0
// For the GC Number conversion
#define GC_KEY                  @"GCNUMBER"
#define DEFAULT_GC_VALUE        @"GC2XZKZ"

double units_to_m[UNITS] = { M_TO_FEET, M_TO_KILO, M_TO_M, M_TO_MILE, M_TO_NAUT, M_TO_YARD, M_TO_FIZZIES };
/*
 
 WAAS (and equivalent) satellite information was slightly out of date.  Information is current as of 8/8/2011 -- Jim Carson
 
 Sources: 
 http://www.kowoma.de/en/gps/waas_egnos.htm
 http://gpsinformation.net/exe/waas.html
 http://en.wikipedia.org/wiki/Wide_Area_Augmentation_System
 http://www.gpsworld.com/gnss-system/augmentation-assistance/news/waas-prn-135-resumes-normal-operation-11243
 http://en.wikipedia.org/wiki/European_Geostationary_Navigation_Overlay_Service
 http://en.wikipedia.org/wiki/Multi-functional_Satellite_Augmentation_System
 
 Real-time map of satellite positions and coverage:
 http://www.nstb.tc.faa.gov/RT_WaasSatelliteStatus.htm
 
 */
struct WAAS {
    NSString *name;
    NSString *shortname;
    NSString *system;
    int prn;
    int nmea;
    double longitude;
};
WAAS Satellites[9] = {
    {@"Anik F1R (USA - central)",@"Anik", @"WAAS",138,51,-107.3},
    {@"ARTEMIS (Africa)",@"Artemis", @"EGNOS",124,37,21.5},
    {@"Inmarsat 3F1 (Indian Ocean)",@"IOR",@"EGNOS",126,39,64.5},
    {@"Inmarsat 3F2 (Atlantic - east)",@"AOR-E",@"EGNOS",120,33,-15.5},
    {@"Inmarsat 4F2 (EMEA)",@"IND-W",@"EGNOS",131,44,25.0},
    {@"Inmarsat 4F3 (USA - east)",@"PAC-E",@"WAAS",133,46,-97.6},
    {@"Intelsat Galaxy 15 (USA - west)",@"G-15",@"WAAS",135,48,-133.1},
    {@"MTSat-1R/Himawari 6 (Japan)",@"M1R",@"MTSat",129,42,140.0},
    {@"MTSat-2/Himawari 7 (Japan)",@"M2",@"MTSat",137,50,145.0}
};

@implementation MacFizzyCalcAppDelegate

@synthesize window;

//
// Convert a coordinate string into various representations, including UTM
//
-(IBAction)doConvert:(id)sender {
    CLatLon *cs = new CLatLon;
    if (cs->ParseCoords([[CC_CoordString stringValue] UTF8String]) == TRUE) { 
        [CC_CoordString setTextColor:[NSColor blackColor]];
        [CC_d setStringValue:[NSString stringWithCString:cs->ToDDD().c_str() 
                                                encoding:[NSString defaultCStringEncoding]]];
        [CC_dm setStringValue:[NSString stringWithCString:cs->ToDMM().c_str() 
                                                 encoding:[NSString defaultCStringEncoding]]];
        [CC_dms setStringValue:[NSString stringWithCString:cs->ToDMS().c_str() 
                                                  encoding:[NSString defaultCStringEncoding]]];
        [CC_utm setStringValue:[NSString stringWithCString:cs->ToUTM().c_str() 
                                                  encoding:[NSString defaultCStringEncoding]]];
        //
        // Ease of use - preserve the values for the coordinates and units
        //
        [D_pt1 setStringValue:[CC_CoordString stringValue]];         
        [P_start setStringValue:[CC_CoordString stringValue]];        
        [W_CoordString setStringValue:[CC_CoordString stringValue]];        
        [CS_CoordString setStringValue:[CC_CoordString stringValue]];

        NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
        NSString *val = [CC_CoordString stringValue];
        if ([val isEqualToString:DEFAULT_COORDS_VALUE]) {
            [defaults removeObjectForKey:COORD_KEY];
        } else {
            [defaults setObject:val forKey:COORD_KEY];
        }
    } else {
        [[sender window] makeFirstResponder:CC_CoordString]; // send focus back to the field in error
        [CC_CoordString setTextColor:[NSColor ERRORCOLOR]];
        [CC_d setStringValue:@""];
        [CC_dm setStringValue:@""];
        [CC_dms setStringValue:@""];
        [CC_utm setStringValue:@""];
    }
    delete cs;
}

//
// Distance between two points.
//
-(IBAction)doDistance:(id)sender {
    CLatLon *pt1 = new CLatLon;
    CLatLon *pt2 = new CLatLon;
    BOOL pt1result, pt2result;
    NSInteger ival = [D_units indexOfSelectedItem];

    pt1result = pt1->ParseCoords([[D_pt1 stringValue] UTF8String]);
    pt2result = pt2->ParseCoords([[D_pt2 stringValue] UTF8String]);

    if (ival == -1) {
        //
        // You would think a pulldown would be straightforwad.  However, if the user inputs something bogus,
        // this will interefere with our cacluations.  If they do, select the default units and continue.
        //
        [D_units selectItemAtIndex:DEFAULT_UNITS_VALUE];
    }
    if ((pt1result == FALSE) || (pt2result == FALSE)) {
        // set the malformatted coordinate text to brown as an attention getter
        if (pt1result == FALSE ) {
            [D_pt1 setTextColor:[NSColor ERRORCOLOR]];
            [[sender window] makeFirstResponder:D_pt1]; // send focus back to the field in error
        } else {
            [D_pt1 setTextColor:[NSColor blackColor]];
        }
        if (pt2result == FALSE ) {
            [D_pt2 setTextColor:[NSColor ERRORCOLOR]];
            [[sender window] makeFirstResponder:D_pt2]; // send focus back to the field in error
        } else {
            [D_pt2 setTextColor:[NSColor blackColor]];        
        }
        [D_distance setStringValue:@""];
        [D_fwd setStringValue:@""];
        [D_rev setStringValue:@""];
    } else {
        double distance = 0.;
        double forwardAzimuth = 0.;
        double reverseAzimuth = 0.;
        // Copy the values to the other, similar attributes
        [CC_CoordString setStringValue:[D_pt1 stringValue]];         
        [P_start setStringValue:[D_pt1 stringValue]];        
        [W_CoordString setStringValue:[D_pt1 stringValue]];        
        [CS_CoordString setStringValue:[D_pt1 stringValue]];

        //
        // Ease of use - preserve the values for the coordinates and units
        //
        NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
        NSInteger ival;
        NSString *val = [D_pt1 stringValue];
        if ([val isEqualToString:DEFAULT_COORDS_VALUE]) {
            [defaults removeObjectForKey:COORD_KEY];
        } else {
            [defaults setObject:val forKey:COORD_KEY];
        }
        val = [D_pt2 stringValue];
        if ([val isEqualToString:DEFAULT_COORDS_VALUE2]) {
            [defaults removeObjectForKey:COORD_KEY2];
        } else {
            [defaults setObject:val forKey:COORD_KEY2];
        }
        ival = [D_units indexOfSelectedItem];
        if (ival == DEFAULT_UNITS_VALUE) {
            [defaults removeObjectForKey:UNITS_KEY];
        } else {
            [defaults setInteger:ival forKey:UNITS_KEY];
        }
        
        [D_pt1 setTextColor:[NSColor blackColor]];  
        [D_pt2 setTextColor:[NSColor blackColor]];  
        // I'm a little surprised I can't fetch the tag directly from NSMatrix versus doing this stuff.
        if ([D_accuracy selectedCell] == [D_accuracy cellWithTag:2]) {
            distance = pt1->RhumbDistance(*pt2, &forwardAzimuth, &reverseAzimuth);
        } else if ([D_accuracy selectedCell] == [D_accuracy cellWithTag:1]) {
            distance = pt1->SphericalDistance(*pt2, &forwardAzimuth,&reverseAzimuth);
        } else {
            distance = pt1->VincentyDistance(*pt2, &forwardAzimuth, &reverseAzimuth);
        }
        
        // Normalize to meters for future calculations
        distance *= units_to_m[[D_units indexOfSelectedItem]];
        
        [D_distance setStringValue:[NSString stringWithFormat:@"%9.6lf", distance]];
        [D_fwd setStringValue:[NSString stringWithFormat:@"%3.6lf", forwardAzimuth]];
        [D_rev setStringValue:[NSString stringWithFormat:@"%3.6lf", reverseAzimuth]];
    }
    delete pt1;
    delete pt2;
}

-(void) SetcDefaults: (NSString *) keyname: (NSString *) cvalue: (NSString  *) dcvalue {
   NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
   if (std::strcmp((const char *)dcvalue, (const char *)cvalue) == 0) {
      [defaults removeObjectForKey:keyname];
   } else {
      [defaults setObject:dcvalue forKey:keyname];    
   }
}
/*
+(void) SetiDefaults: (NSString *) keyname: (NSInteger) ivalue: (NSInteger) divalue {
   NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
   if (divalue == ivalue) {
      [defaults removeObjectForKey:keyname];
   } else {
      [defaults setInteger:ivalue forKey:keyname];    
   }    
}

+(void) setDefaultCoords: (NSString *) coordstring {
    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    NSString *val = coordstring;
    if ([val isEqualToString:DEFAULT_COORDS_VALUE]) {
        [defaults removeObjectForKey:COORD_KEY];
    } else {
        [defaults setObject:val forKey:COORD_KEY];
    }
}
*/
// 
// Project a point from a point.
//
// to do: error checking for the distance and azimuth fields?
//
-(IBAction)doProjection:(id)sender {
    CLatLon *startpoint = new CLatLon;    
    NSInteger ival = [P_units indexOfSelectedItem];
    if (ival == -1) {
        // user has given us some bogus units.
        [P_units selectItemAtIndex:DEFAULT_UNITS_VALUE];        
    }

    if (startpoint->ParseCoords([[P_start stringValue] UTF8String]) == FALSE) {
        // NSLog(@"Start point is not in a recognized format.\n");
        [[sender window] makeFirstResponder:P_start]; // send focus back to the field in error
        [P_start setTextColor:[NSColor ERRORCOLOR]];
        [P_result setStringValue:@""];
    } else {
        double distance = [P_distance doubleValue];
        double azimuth = [P_bearing doubleValue];
        // Update the default units if the user has set them.
        [CC_CoordString setStringValue:[P_start stringValue]];         
        [D_pt1 setStringValue:[P_start stringValue]];         
        [W_CoordString setStringValue:[P_start stringValue]];        
        [CS_CoordString setStringValue:[P_start stringValue]];
        NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
        NSInteger ival;
        NSString *val = [P_start stringValue];
        ival = [P_units indexOfSelectedItem];
        if (ival == DEFAULT_UNITS_VALUE) {
            [defaults removeObjectForKey:UNITS_KEY];
        } else {
            [defaults setInteger:ival forKey:UNITS_KEY];
        }
        if ([val isEqualToString:DEFAULT_COORDS_VALUE]) {
            [defaults removeObjectForKey:COORD_KEY];
        } else {
            [defaults setObject:val forKey:COORD_KEY];
        }
        [P_start setTextColor:[NSColor blackColor]];
        
        // Unlike the distance calculation where we render in the other units, the purpose here
        // is to standardize to meters.
        distance /= units_to_m[[P_units indexOfSelectedItem]];
        
        // I'm a little surprised I can't fetch the tag directly from NSMatrix versus doing this stuff.
        if ([P_accuracy selectedCell] == [P_accuracy cellWithTag:2]) {
            [P_result setStringValue:[NSString stringWithCString:startpoint->RhumbProjection(azimuth, distance).ToDMM().c_str() encoding:[NSString defaultCStringEncoding]]];
        } else if ([P_accuracy selectedCell] == [P_accuracy cellWithTag:1]) {
            [P_result setStringValue:[NSString stringWithCString:startpoint->SphericalProjection(azimuth, distance).ToDMM().c_str() encoding:[NSString defaultCStringEncoding]]];
        } else {
            [P_result setStringValue:[NSString stringWithCString:startpoint->VincentyProjection(azimuth, distance).ToDMM().c_str() encoding:[NSString defaultCStringEncoding]]]; 
        }
    }
    delete startpoint;
}

//
// Calculate the location of a WAAS satellite from the specific point.
//
-(IBAction)doWAAS:(id)sender {
    CLatLon *cs = new CLatLon;
    NSInteger ival = [W_satellite indexOfSelectedItem];
    
    if (ival == -1) {
        //
        // 08/28/2011 - jim - user selected a value outside of our array, perhaps through
        // typing bogus information in the field.  Select the default value in the dialog, 
        // but continue.
        //
        [W_satellite selectItemAtIndex:DEFAULT_SAT_VALUE];
    }

    if (cs->ParseCoords([[W_CoordString stringValue] UTF8String])== FALSE) {
        // NSLog(@"WAAS coordinate is not in a recognized format.\n");
        [[sender window] makeFirstResponder:W_CoordString]; // send focus back to the field in error
        [W_CoordString setTextColor:[NSColor ERRORCOLOR]];
    } else {
        [W_CoordString setTextColor:[NSColor blackColor]];
        double Azimuth, Elevation;
        double dDistance;
        double dSatLon = Satellites[[W_satellite indexOfSelectedItem]].longitude;
        NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
        NSInteger ival;
        NSString *val = [W_CoordString stringValue];
        // Update the default units if the user has set them.
        [CC_CoordString setStringValue:[W_CoordString stringValue]];         
        [D_pt1 setStringValue:[W_CoordString stringValue]];         
        [P_start setStringValue:[W_CoordString stringValue]];        
        [CS_CoordString setStringValue:[W_CoordString stringValue]];

        ival = [W_satellite indexOfSelectedItem];
        if (ival == DEFAULT_SAT_VALUE) {
            [defaults removeObjectForKey:SAT_KEY];
        } else {
            [defaults setInteger:ival forKey:SAT_KEY];
        } 
        if ([val isEqualToString:DEFAULT_COORDS_VALUE]) {
            [defaults removeObjectForKey:COORD_KEY];
        } else {
            [defaults setObject:val forKey:COORD_KEY];
        }
        
        if ([W_accuracy selectedCell] == [W_accuracy cellWithTag:1]) {
            dDistance = cs->GeoSatelliteAzElSpherical(dSatLon, 0., &Azimuth, &Elevation);
        } else {
            dDistance = cs->GeoSatelliteAzEl(dSatLon, 0., &Azimuth, &Elevation);
        }
        /**/
        
        [W_azimuth setStringValue:[NSString stringWithFormat:@"%3.2lf", Azimuth]];
        
        if (Elevation > 0.) {
            [W_elevation setStringValue:[NSString stringWithFormat:@"%3.2lf", Elevation]];
            [W_elevation setTextColor:[NSColor blackColor]];
        } else {
            [W_elevation setStringValue:@"Below horizon"];
            [W_elevation setTextColor:[NSColor ERRORCOLOR]];
        }
        
        [W_satname setStringValue:Satellites[[W_satellite indexOfSelectedItem]].shortname];
        [W_satnum setIntValue:Satellites[[W_satellite indexOfSelectedItem]].nmea];
        [W_system setStringValue:Satellites[[W_satellite indexOfSelectedItem]].system];
        [W_satlong setIntValue:dSatLon];
    }
    delete cs;
}

-(IBAction)doGCNumber:(id)sender {
    [GC_number setTextColor:[NSColor ERRORCOLOR]];
    [GC_number setTextColor:[NSColor blackColor]];
    
    unsigned long ulNumber = 0;
    CLatLon *cs = new CLatLon;
    std::string RetVal;
    ulNumber = cs->GetNumber([[GC_input stringValue] UTF8String]);
    if (ulNumber > 0) {
        [GC_number setStringValue:[NSString stringWithFormat:@"%u", ulNumber]];
    }
    RetVal = cs->MakeWaypoint(ulNumber);
    [GC_waypoint setStringValue:[NSString stringWithCString:RetVal.c_str() encoding:[NSString defaultCStringEncoding]]];
    delete cs;
}

-(IBAction)doChecksum:(id)sender {
    CLatLon *cs = new CLatLon;
    std::string tmp = [[CS_CoordString stringValue] UTF8String];
    if (cs->ParseCoords(tmp) == TRUE) { 
        // do we care to set the other coordinates?
        [CS_CoordString setTextColor:[NSColor blackColor]];
        std::string strLatString, strLonString;
        cs->SplitCoordString(tmp, strLatString, strLonString);
        [CS_lat_checksum setIntValue:cs->GetChecksum(strLatString)];
        [CS_lon_checksum setIntValue:cs->GetChecksum(strLonString)];
        [CS_entire_checksum setIntValue:(int) cs->GetChecksum(strLatString) + cs->GetChecksum(strLonString)];        
        [CS_lat_digiroot setIntValue:1+(([CS_lat_checksum intValue]-1)%9)];
        [CS_lon_digiroot setIntValue:1+(([CS_lon_checksum intValue]-1)%9)];
        [CS_entire_digiroot  setIntValue:1+(([CS_entire_checksum intValue]-1)%9)];  
    } else {
        // NSLog(@"Coordinate is not in a recognized format.\n");
        [[sender window] makeFirstResponder:CS_CoordString]; // send focus back to the field in error
        [CS_CoordString setTextColor:[NSColor ERRORCOLOR]];
        [CS_lat_checksum setStringValue:@"-"];
        [CS_lon_checksum setStringValue:@"-"];
        [CS_entire_checksum setStringValue:@"-"];
        [CS_lat_digiroot setStringValue:@"-"];
        [CS_lon_digiroot setStringValue:@"-"];
        [CS_entire_digiroot  setStringValue:@"-"]; 
    }
    delete cs;
}

-(BOOL) applicationShouldTerminateAfterLastWindowClosed:(NSApplication *) theApp; {return TRUE; }
+ (void)initialize{
    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    NSDictionary *appDefaults = [NSDictionary
                                 dictionaryWithObject:@"YES" forKey:@"D_units"];
    [defaults registerDefaults:appDefaults];
}


#define DEFAULT_UNITS_VALUE 2 // meters
#define UNITS_KEY @"UNITS"

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
    // Insert code here to initialize your application
    NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
    
    //
    // This is all about setting up the initial user dialog values.
    //
    NSUserDefaults *defaults = [NSUserDefaults standardUserDefaults];
    NSInteger   unitsval = [defaults integerForKey:UNITS_KEY];
    NSInteger   satval  = [defaults integerForKey:SAT_KEY];
    NSString    *coordsval = [defaults stringForKey:COORD_KEY];
    NSString    *coordsval2 =[defaults stringForKey:COORD_KEY2]; 
    // 10/06/2011 - these aren't used.  
    //if (unitsval == nil) unitsval = DEFAULT_UNITS_VALUE;
    //if (satval == nil) satval = DEFAULT_SAT_VALUE;
    if (coordsval == nil) coordsval = DEFAULT_COORDS_VALUE;
    if (coordsval2 == nil) coordsval2 = DEFAULT_COORDS_VALUE2;
    [CC_CoordString setStringValue:coordsval];
    [D_pt1 setStringValue:coordsval];
    [D_pt2 setStringValue:coordsval2];
    [D_units selectItemAtIndex:unitsval];
    [P_start setStringValue:coordsval];
    [P_units selectItemAtIndex:unitsval];
    [P_distance setStringValue:DEFAULT_DIST_VALUE];
    [P_bearing  setStringValue:DEFAULT_BEARING_VALUE];
    [W_CoordString setStringValue:coordsval];
    [W_satellite selectItemAtIndex:satval];
    [CS_CoordString setStringValue:coordsval];
    [GC_input setStringValue:DEFAULT_GC_VALUE];
    
    [pool drain];
}

@end
