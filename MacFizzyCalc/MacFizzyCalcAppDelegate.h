/****************************************************************************
 FizzyCalcAppDelegate.h
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
#import <Cocoa/Cocoa.h>
#import <Foundation/Foundation.h>
#import <stdio.h>
#import <stdlib.h>
#import <string.h>
#import <iostream>
#import <string>
#import <AppKit/NSColor.h>

#define ERRORCOLOR brownColor

@interface MacFizzyCalcAppDelegate : NSObject <NSApplicationDelegate> {
    NSWindow *window;
    // Convert tab
    IBOutlet NSTextField    *CC_CoordString;
    IBOutlet NSTextField    *CC_d;
    IBOutlet NSTextField    *CC_dm;
    IBOutlet NSTextField    *CC_dms;
    IBOutlet NSTextField    *CC_utm;
    IBOutlet NSButton       *CC_go;
    IBOutlet NSHelpManager  *CC_help;
    
    // Distance tab
    IBOutlet NSTextField    *D_pt1;
    IBOutlet NSTextField    *D_pt2;
    IBOutlet NSMatrix       *D_accuracy;
    IBOutlet NSTextField    *D_distance;
    IBOutlet NSComboBox     *D_units;
    IBOutlet NSTextField    *D_fwd;
    IBOutlet NSTextField    *D_rev;
    
    // Projection tab
    IBOutlet NSTextField    *P_start;
    IBOutlet NSMatrix       *P_accuracy;
    IBOutlet NSTextField    *P_distance;
    IBOutlet NSComboBox     *P_units;
    IBOutlet NSTextField    *P_bearing;
    IBOutlet NSTextField    *P_result;
    
    // WAAS tab
    IBOutlet NSTextField    *W_CoordString;
    IBOutlet NSMatrix       *W_accuracy;
    IBOutlet NSComboBox     *W_satellite;
    IBOutlet NSTextField    *W_satnum;
    IBOutlet NSTextField    *W_satlong;
    IBOutlet NSTextField    *W_satname;
    IBOutlet NSTextField    *W_azimuth;
    IBOutlet NSTextField    *W_system;
    IBOutlet NSTextField    *W_elevation;
    
    // GC number tab
    IBOutlet NSTextField    *GC_input;
    IBOutlet NSTextField    *GC_waypoint;
    IBOutlet NSTextField    *GC_number;
    
    // Checksum tab
    IBOutlet NSTextField    *CS_CoordString;
    IBOutlet NSTextField    *CS_lat_checksum;
    IBOutlet NSTextField    *CS_lon_checksum;
    IBOutlet NSTextField    *CS_entire_checksum;
    IBOutlet NSTextField    *CS_lat_digiroot;
    IBOutlet NSTextField    *CS_lon_digiroot;
    IBOutlet NSTextField    *CS_entire_digiroot;
}

@property (assign) IBOutlet NSWindow *window;
// process "go" button events
-(IBAction)doConvert:(id)sender;
-(IBAction)doDistance:(id)sender;
-(IBAction)doProjection:(id)sender;
-(IBAction)doWAAS:(id)sender;
-(IBAction)doGCNumber:(id)sender;
-(IBAction)doChecksum:(id)sender;
@end



