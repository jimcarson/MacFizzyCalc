// Description:  Coordinate calculations.

/* Coord.cpp
 * 
 * Copyright (C) 2003 by David Knapp
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *
 * 8/4/2011 - Jim Carson (jim.carson@gmail.com) - added functions GetNumber,
 * MakeWaypoint and moved several constants in the original MS dialogs into this file.  
 * Replaced utility functions available only in MS with C++ standards. 
 * 
 */

#include "Coord.h"

// Ellipsoid is initialized to WGS84 ellipsoid.
const CEllipsoid CLatLon::m_Ellipsoid = CEllipsoid(6378137.00, 298.257223563);
const double CLatLon::m_Radius = 6366707.01896486;// Earth's radius, in meters.
const double CLatLon::m_GeosyncRadius = 42164211.0;
const double CLatLon::m_Deg2Rad = 1.74532925199433E-02;

// UTM zone letters
const std::string CLatLon::m_strZoneLetters = "CDEFGHJKLMNPQRSTUVWX";

// used in GC number calculation.  
std::string Base31 = "0123456789ABCDEFGHJKMNPQRTVWXYZ";
#define ULTRANSITION    65535
#define ULOFFSET        411120


//
// Explanation of GC Numbers: http://support.groundspeak.com/index.php?pg=kb.page&id=221
//
unsigned long CLatLon::GetNumber(std::string strWaypoint)
{
    unsigned long ulNumber = 0;
    bool f_Numeric = true;
    bool f_Invalid = false;
    int i;
    std::string strID = strWaypoint;    
    std::transform(strID.begin(), strID.end(),strID.begin(), ::toupper);
    
    for (i=0; i < strID.length(); i++) {
        if (!isdigit(strID[i])) {
            f_Numeric = false;
            if (Base31.rfind(strID[i]) == std::string::npos) {
                f_Invalid = true;
            }
        }
    }
    if (f_Invalid) {        // If there are non-Base31 characters in it, discard the result
        ulNumber = 0;
    } else if (f_Numeric) { // If it's all digits, the user has specified a number to convert to waypoint.
        if (strID[0] == 'G' && strID[1] == 'C') {
            strID = strID.substr(2,strID.length());
        }
        ulNumber = std::strtoul(strID.c_str(),NULL,10);
    } else {
        //
        // If the first two characters are "GC" then remove them.
        //
        if (strID[0] == 'G' && strID[1] == 'C') {
            strID = strID.substr(2,strID.length());
        }
        while (strID.length() < 4) { // pad the GC
            strID = "0" + strID;
        }
        //
        // For values GC0 through GCFFFF, the number is in hex.  Afterwards, it's base 31.
        //
        if (strID.length() == 4 && strID.compare("FFFF") <= 0) {
            ulNumber = std::strtoul(strID.c_str(),NULL,16);
        } else {
            ulNumber = 0;
            for (i=0; i < strID.length(); i++) {
                ulNumber = ulNumber * 31 + Base31.rfind(strID[i]);
            }
            ulNumber -= ULOFFSET;
        }
    }
    return ulNumber;
}

std::string CLatLon::MakeWaypoint(unsigned long ulNumber)
{
    std::string RetVal = "GC";
    std::string tmp;
    char xxx[10];
    if (ulNumber <= ULTRANSITION) {
        sprintf(xxx,"%lx",ulNumber);
        RetVal.append(xxx);
    } else {
        unsigned long ulTmp = ulNumber + ULOFFSET;
        unsigned long iVal;
        while (ulTmp >= 1) {
            iVal = ulTmp % 31;
            tmp += Base31[iVal];
            ulTmp /= 31;
        }
        std::reverse(tmp.begin(), tmp.end());
        RetVal.append(tmp);
    }
    std::transform(RetVal.begin(), RetVal.end(),RetVal.begin(), ::toupper);
    return RetVal;
}


// -------------------------------------------------------------------------
// METHOD:  CLatLon::SphericalDistance()
/*! 
 \brief  Computes great-circle distance from this point to point P.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [double] - Distance between this point and P in meters.
 
 \param P [CLatLon&] - Point to which to compute distance.
 
 */
// -------------------------------------------------------------------------
double CLatLon::SphericalDistance(CLatLon& P)
{
    double dTmp1, dTmp2;
    return SphericalDistance(P, &dTmp1, &dTmp2);
}


// -------------------------------------------------------------------------
// METHOD:  CLatLon::SphericalDistance()
/*! 
 \brief  Computes great-circle distance from this point to point P.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [double] - Distance between this point and P in meters.
 
 \param P [CLatLon&] - Point to which to compute distance.
 \param pForwardAzimuth [double *] - Pointer to parameter to receive 
 forward azimuth in degrees.
 \param pReverseAzimuth [double *] - Pointer to parameter to receive
 reverse azimuth in degrees.
 */
// -------------------------------------------------------------------------
double CLatLon::SphericalDistance(CLatLon& P, double *pForwardAzimuth, double *pReverseAzimuth)
{
    double dDeltaLong = (m_Longitude - P.m_Longitude) * m_Deg2Rad;
    double dLat1 = m_Latitude * m_Deg2Rad;
    double dLat2 = P.m_Latitude * m_Deg2Rad;
    double dAngle = acos(sin(dLat1) * sin(dLat2) 
                         + cos(dLat1) * cos(dLat2) * cos(dDeltaLong));
    double dDistance = m_Radius * dAngle;
    
    *pForwardAzimuth = atan2(sin(-dDeltaLong)*cos(dLat2), cos(dLat1)*sin(dLat2)-sin(dLat1)*cos(dLat2)*cos(-dDeltaLong)) / m_Deg2Rad;
    if (*pForwardAzimuth < 0.) *pForwardAzimuth += 360.;
    
    *pReverseAzimuth = atan2(sin(dDeltaLong)*cos(dLat1), cos(dLat2)*sin(dLat1)-sin(dLat2)*cos(dLat1)*cos(dDeltaLong)) / m_Deg2Rad;
    if (*pReverseAzimuth < 0.) *pReverseAzimuth += 360.;
    
    // The above formula doesn't work well for small distances.
    // So if the distance is small, use simple linear approximation.
    if (dDistance < .01) {
        double dDeltaLat = dLat1 - dLat2;
        dDeltaLong *= cos(dLat2);
        dDistance = m_Radius * sqrt(dDeltaLat*dDeltaLat + dDeltaLong*dDeltaLong);
    }
    return dDistance;
}

double CLatLon::ApproxEllipsoidDistance(CLatLon& P) 
{
    double dDeltaLong = (m_Longitude - P.m_Longitude) * m_Deg2Rad;
    double dLat1 = m_Latitude * m_Deg2Rad;
    double dLat2 = P.m_Latitude * m_Deg2Rad;
    double dAngle = acos(sin(dLat1) * sin(dLat2) 
                         + cos(dLat1) * cos(dLat2) * cos(dDeltaLong));
    
    double a0 = m_Ellipsoid.m_a;
    double flat = 1./m_Ellipsoid.m_fInv;
    double r = 1. - flat;
    double b0 = a0 * r;
    
    double dGCLat1 = atan(r*r*tan(dLat1));
    double dGCLat2 = atan(r*r*tan(dLat2));
    double dAveLat = (dGCLat1 + dGCLat2) / 2.;
    double dAvgR = 1./ sqrt(cos(dAveLat)*cos(dAveLat)/(a0*a0) + sin(dAveLat)*sin(dAveLat)/(b0*b0));
    
    double dDistance = dAvgR * dAngle;
    
    return dDistance;
}


// -------------------------------------------------------------------------
// METHOD:  CLatLon::SphericalProjection()
/*! 
 \brief  Computes a projection along a great circle.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [CLatLon] - Projected coordinates.
 
 \param dAzimuth [double] - Forward azimuth in degrees.
 \param dDistance [double] - Distance in meters.
 */
// -------------------------------------------------------------------------
CLatLon CLatLon::SphericalProjection(double dAzimuth, double dDistance)
{
    CLatLon RetCoord;
    dAzimuth *= m_Deg2Rad;
    
    double dLat1 = m_Deg2Rad * m_Latitude;
    double dLong1 = m_Deg2Rad * m_Longitude;
    double s = dDistance / m_Radius;
    
    double dLat2 = asin(sin(dLat1)*cos(s) + cos(dLat1)*sin(s)*cos(dAzimuth));
    double dLong2 = dLong1 + atan2(sin(dAzimuth)*sin(s)*cos(dLat1), cos(s) - sin(dLat1)*sin(dLat2));
    RetCoord.m_Latitude = NormalizeLatitude(dLat2 / m_Deg2Rad);
    RetCoord.m_Longitude = NormalizeLongitude(dLong2 / m_Deg2Rad);
    
    return RetCoord;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::SphericalDistance()
/*! 
 \brief  Computes great-circle distance from this point to the great circle
 connecting P1 and P2.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [double] - Distance between this point and the GC in meters.
 If this point is not between the other points, the
 returned distance will be negative.
 
 \param P1 [CLatLon&] - First point defining great circle.
 \param P2 [CLatLon&] - Second point defining great circle.
 */
// -------------------------------------------------------------------------
double CLatLon::SphericalDistance(CLatLon& P1, CLatLon& P2)
{
    double dSign = 1.;
    if (!IsBetween(P1, P2)) {
        dSign = -1.;
    }
    CCartesianCoord P = ToSphericalCartesian();
    CCartesianCoord N = P1.SphericalCross(P2);
    N.Normalize();
    return dSign * fabs(m_Radius * (PI/2. - acos(N.dot(P))));
}


// -------------------------------------------------------------------------
// METHOD:  CLatLon::RhumbDistance()
/*! 
 \brief  Computes rhumb-line distance from this point to point P.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [double] - Distance between this point and P in meters.
 
 \param P [CLatLon&] - Point to which to compute distance.
 */
// -------------------------------------------------------------------------
double CLatLon::RhumbDistance(CLatLon& P)
{
    double dTmp1, dTmp2;
    return RhumbDistance(P, &dTmp1, &dTmp2);
}


// -------------------------------------------------------------------------
// METHOD:  CLatLon::RhumbDistance()
/*! 
 \brief  Computes rhumb-line distance from this point to point P.
 
 \author fizzymagic
 \date   7/3/2005
 
 \return  [double] - Distance between this point and P in meters.
 
 \param P [CLatLon&] - Point to which to compute distance.
 \param pForwardAzimuth [double *] - Pointer to parameter to receive 
 forward azimuth in degrees.
 \param pReverseAzimuth [double *] - Pointer to parameter to receive
 reverse azimuth in degrees.
 */
// -------------------------------------------------------------------------
double CLatLon::RhumbDistance(CLatLon& P, double *pForwardAzimuth, double *pReverseAzimuth)
{
    double dDeltaLong = (P.m_Longitude - m_Longitude) * m_Deg2Rad;
    double dLat1 = m_Latitude * m_Deg2Rad;
    double dLat2 = P.m_Latitude * m_Deg2Rad;
    double dDeltaLat = dLat2 - dLat1;
    double dQ, dPhi = 0.;
    
    if (fabs(dDeltaLat) < sqrt(EPSILON)) {
        dQ = cos(dLat1);
        *pForwardAzimuth = (dDeltaLong > 0)?90.:270.;
    } 
    else {
        dPhi = log(tan(dLat2/2. + PI/4.)/tan(dLat1/2. + PI/4.));
        dQ = (dDeltaLat)/dPhi;
        *pForwardAzimuth = atan2(dDeltaLong, dPhi) / m_Deg2Rad;
    }
    
    double dDistance = m_Radius * sqrt(dDeltaLat*dDeltaLat + dQ*dQ*dDeltaLong*dDeltaLong);
    
    if (*pForwardAzimuth < 0.) *pForwardAzimuth += 360.;
    
    *pReverseAzimuth = *pForwardAzimuth - 180.;
    if (*pReverseAzimuth < 0.) *pReverseAzimuth += 360.;
    
    return dDistance;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::RhumbProjection()
/*! 
 \brief  Computes a projection along a rhumb line.
 
 \author fizzymagic
 \date   7/3/2005
 
 \return  [CLatLon] - Projected coordinates.
 
 \param dAzimuth [double] - Forward azimuth in degrees.
 \param dDistance [double] - Distance in meters.
 */
// -------------------------------------------------------------------------
CLatLon CLatLon::RhumbProjection(double dAzimuth, double dDistance)
{
    CLatLon RetCoord;
    dAzimuth *= m_Deg2Rad;
    
    double dLat1 = m_Deg2Rad * m_Latitude;
    double dLong1 = m_Deg2Rad * m_Longitude;
    double s = dDistance / m_Radius;
    double dLat2 = dLat1 + s * cos(dAzimuth);
    double dDeltaLat = dLat2 - dLat1;
    double dQ, dPhi = 0.;
    
    if (fabs(dDeltaLat) < sqrt(EPSILON)) {
        dQ = cos(dLat1);
    } 
    else {
        dPhi = log(tan(dLat2/2. + PI/4.)/tan(dLat1/2. + PI/4.));
        dQ = (dDeltaLat)/dPhi;
    }
    
    double dLong2 = dLong1 + s*sin(dAzimuth)/dQ;
    RetCoord.m_Latitude = NormalizeLatitude(dLat2 / m_Deg2Rad);
    RetCoord.m_Longitude = NormalizeLongitude(dLong2 / m_Deg2Rad);
    
    return RetCoord;
}


// -------------------------------------------------------------------------
// METHOD:  CLatLon::IsBetween()
/*! 
 \brief  Approximate test to check if this point is between P1 and P2.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [bool] - true if point it between others; false otherwise.
 
 \param P1 [CLatLon&] - First point.
 \param P2 [CLatLon&] - Second point.
 */
// -------------------------------------------------------------------------
bool CLatLon::IsBetween(CLatLon& P1, CLatLon& P2)
{
    // Approximate scale factor for longitudes.
    double cosLat = cos(m_Latitude * m_Deg2Rad);
    double u = (m_Longitude - P1.m_Longitude) * (P2.m_Longitude - P1.m_Longitude) * cosLat * cosLat
    + (m_Latitude - P1.m_Latitude) * (P2.m_Latitude - P1.m_Latitude);
    u /= (P2.m_Longitude - P1.m_Longitude) * (P2.m_Longitude - P1.m_Longitude) * cosLat * cosLat
    + (P2.m_Latitude - P1.m_Latitude) * (P2.m_Latitude - P1.m_Latitude);
    
    // u is the relative position on the line connecting P1 and P2.
    return (u >= 0. && u < 1.);
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ParseCoords()
/*! 
 \brief  Parses single coordinate string in a variety of formats.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [bool] - true if coordinates were parsed without error.
 
 \param szCoordString [const char *] - Coordinate string.
 */
// -------------------------------------------------------------------------
bool CLatLon::ParseCoords(const char *szCoordString)
{
    if (ParseUTM(szCoordString)) {
        return true;
    }
    
    std::string strCoordString(szCoordString);
    CleanCoordString(strCoordString);
    std::string strLatString, strLongString;
    if (!SplitCoordString(strCoordString, strLatString, strLongString)) {
        return false;
    }
    return ParseCoords(strLatString.c_str(), strLongString.c_str());
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ParseCoords()
/*! 
 \brief  Parses separate lat/long coordinate strings in a variety of formats.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [bool] - true if coordinates were parsed without error.
 
 \param szLatString [const char *] - Coordinate string for latitude.
 \param szLongString [const char *] - Coordinate string for longitude.
 */
// -------------------------------------------------------------------------
bool CLatLon::ParseCoords(const char *szLatString, const char *szLongString)
{
    std::string strLatString(szLatString), strLongString(szLongString);
    CleanCoordString(strLatString);
    CleanCoordString(strLongString);
    
    int iLatSign = 1;
    std::string::iterator it = strLatString.begin();
    std::string strBuffer;
    double dLatitude, dLongitude;
    
    if (strLatString.find('N') != std::string::npos)  {
        findandreplace(strLatString, "N", "");
    }
    if (strLatString.find('S') != std::string::npos)  {
        iLatSign = -1;
        findandreplace(strLatString,"S", "");
    }
    while (it < strLatString.end() && ::isspace(*it)) it++;
    
    strBuffer.erase();
    strBuffer.append(it, strLatString.end());
    if (ParseDegreeString(strBuffer, &dLatitude)) {
        dLatitude *= iLatSign;
    }
    else return false;
    
    int iLongSign = 1;
    it = strLongString.begin();
    if (strLongString.find('E') != std::string::npos)  {
        findandreplace(strLongString, "E", "");
    }
    if (strLongString.find('W') != std::string::npos)  {
        iLongSign = -1;
        findandreplace(strLongString,"W", "");
    }
    while (it < strLongString.end() && ::isspace(*it)) it++;
    
    strBuffer.erase();
    strBuffer.append(it, strLongString.end());
    if (ParseDegreeString(strBuffer, &dLongitude)) {
        dLongitude *= iLongSign;
    }
    else return false;
    
    m_Latitude = dLatitude;
    m_Longitude = dLongitude;
    
    return true;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ParseDegreeString()
/*! 
 \brief  Parses a string containing degrees in various formats.
 
 \author fizzymagic
 \date   9/8/2003
 
 \return  [bool] - true if successful parse, false otherwise.
 
 \param szDegrees [char *] - Degree string to parse.
 \param dResult [double *] - Result of parsed string, in degrees.
 */
// -------------------------------------------------------------------------
bool CLatLon::ParseDegreeString(const char *szDegrees, double *dResult)
{
    std::string strDegrees(szDegrees);
    std::string::iterator it;
    double dTemp, dDegrees = 0.;
    int iSign = 1;
    std::string strBuffer;
    
    CleanCoordString(strDegrees);
    
    it = strDegrees.begin();
    if (it >= strDegrees.end() || !(::isdigit(*it) || *it == '-')) {
        *dResult = 0.;
        return false;
    }
    if (*it == '-') {
        iSign = -1;
        it++;
    }
    
    strBuffer.erase();
    while (it < strDegrees.end() && (::isdigit(*it) || *it == '.')) {
        if (it >= strDegrees.end()) return false;
        strBuffer.append(1, *it);
        it++;
    }
    dTemp = atof(strBuffer.c_str());
    if (dTemp < -180. || dTemp > 360.) {
        return false;
    }
    else {
        dDegrees = dTemp;
    }
    while (it < strDegrees.end() && ::isspace(*it)) it++;
    if (it < strDegrees.end()) {
        if (!::isdigit(*it)) return false;
        strBuffer.erase();
        while (it < strDegrees.end() && (::isdigit(*it) || *it == '.')) {
            strBuffer.append(1, *it);
            it++;
        }
        dTemp = atof(strBuffer.c_str());
        if (dTemp < 0. || dTemp >= 60.) {
            return false;
        }
        else {
            dDegrees += dTemp / 60.;
        }
        while (it < strDegrees.end() && ::isspace(*it)) it++;
        if (it < strDegrees.end() && *it != 0) {
            if (!::isdigit(*it)) return false;
            strBuffer.erase();
            while (it < strDegrees.end() && (::isdigit(*it) || *it == '.')) {
                strBuffer.append(1, *it);
                it++;
            }
            dTemp = atof(strBuffer.c_str());
            if (dTemp < 0. || dTemp >= 60.) {
                return false;
            }
            else {
                dDegrees += dTemp / 3600.;
            }
        }
    }
    *dResult = dDegrees * iSign;
    return true;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ParseUTM()
/*! 
 \brief  Parses UTM coordinate strings.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [bool] - true if the coordinate string is UTM and as parsed
 without error.
 
 \param szCoordString [const char *] - UTM coordinate string.
 */
// -------------------------------------------------------------------------
bool CLatLon::ParseUTM(const char *szCoordString)
{
    std::string strCoordString(szCoordString);
    std::transform(strCoordString.begin(), strCoordString.end(), strCoordString.begin(), ::toupper);
    CleanCoordString(strCoordString);
    
    std::string::iterator it = strCoordString.begin();
    std::string strBuffer;
    
    int iZone, iZoneLetter;
    double dEasting, dNorthing;
    
    if (!::isdigit(*it)) return false;
    strBuffer.erase();
    while (it < strCoordString.end() && ::isdigit(*it)) {
        strBuffer.append(1, *it);
        it++;
    }
    iZone = atoi(strBuffer.c_str());
    if (it >= strCoordString.end()) return false;
    if ((iZoneLetter = (int)m_strZoneLetters.find(*(it++))) == std::string::npos) return false;
    
    strBuffer.erase();
    while (it < strCoordString.end() && !::isdigit(*it)) it++;
    if (it >= strCoordString.end()) return false;
    while (it < strCoordString.end() && ::isdigit(*it)) {
        strBuffer.append(1, *it);
        it++;
    }
    if (it >= strCoordString.end()) return false;
    dEasting = atof(strBuffer.c_str());
    
    strBuffer.erase();
    while (it < strCoordString.end() &&  !::isdigit(*it)) it++;
    if (it >= strCoordString.end()) return false;
    while (it < strCoordString.end() && ::isdigit(*it)) {
        strBuffer.append(1, *it);
        it++;
    }
    dNorthing = atof(strBuffer.c_str());
    
    return ConvertUTM(iZone, m_strZoneLetters[iZoneLetter], dEasting, dNorthing);
}

double CLatLon::NormalizeLatitude(double dLatitude)
{
    if (dLatitude < -90.) dLatitude = -90.;
    else if (dLatitude > 90.) dLatitude = 90.;
    return dLatitude;
}

double CLatLon::NormalizeLongitude(double dLongitude)
{
    double dLongTemp = dLongitude + 180.;
    while (dLongTemp < 0.) {
        dLongTemp += 360.;
    }
    return fmod(dLongTemp, 360.) - 180.;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::LatToDDD()
/*! 
 \brief  Converts latitude coordinate into DD.DDDDDD format.
 
 \author fizzymagic7
 \date   9/6/2003
 
 \return  [std::string] - Latitude in DD.DDDDDD format.
 */
// -------------------------------------------------------------------------
std::string CLatLon::LatToDDD(void)
{
    std::ostringstream Stream;
    Stream << std::fixed << std::setprecision(6) 
    << m_Latitude;
    return Stream.str();
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::LatToDMM()
/*! 
 \brief  Converts latitude coordinate into DD MM.MMM format.
 
 \author fizzymagic7
 \date   9/6/2003
 
 \return  [std::string] - Latitude in DD MM.MMM format.
 */
// -------------------------------------------------------------------------
std::string CLatLon::LatToDMM(void)
{
    std::ostringstream Stream;
    int iLatDegrees = (int) floor(fabs(m_Latitude));
    double dLatMinutes = 60. * (fabs(m_Latitude) - floor(fabs(m_Latitude)));
    Stream << std::fixed << std::setprecision(3)
    << ((m_Latitude >= 0.)?'N':'S') << " "
    << iLatDegrees  << " "
    << std::setfill('0') << std::setw(6) << dLatMinutes;
    return Stream.str();
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::LatToDMS()
/*! 
 \brief  Converts latitude coordinate into DD MM SS.SS format.
 
 \author fizzymagic7
 \date   9/6/2003
 
 \return  [std::string] - Latitude in DD MM SS.SS format.
 */
// -------------------------------------------------------------------------
std::string CLatLon::LatToDMS(void)
{
    std::ostringstream Stream;
    int iLatDegrees = (int) floor(fabs(m_Latitude));
    int iLatMinutes = (int) floor(60. * (fabs(m_Latitude) - floor(fabs(m_Latitude))) + COORD_EPSILON);
    double dLatSeconds = 60. * (fabs(m_Latitude * 60.) - floor(fabs(m_Latitude * 60.)));
    Stream << std::fixed << std::setprecision(2)
    << ((m_Latitude >= 0.)?'N':'S') << " "
    << iLatDegrees  << " "
    << std::setfill('0') << std::setw(2) << iLatMinutes  << " "
    << std::setfill('0') << std::setw(5) << dLatSeconds;
    return Stream.str();
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::LongToDDD()
/*! 
 \brief  Converts longitude coordinate into DDD.DDDDDD format.
 
 \author fizzymagic7
 \date   9/6/2003
 
 \return  [std::string] - Longitude in DDD.DDDDDD format.
 */
// -------------------------------------------------------------------------
std::string CLatLon::LongToDDD(void)
{
    std::ostringstream Stream;
    Stream << std::fixed << std::setprecision(6) 
    << m_Longitude;
    return Stream.str();
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::LongToDMM()
/*! 
 \brief  Converts longitude coordinate into DDD MM.MMM format.
 
 \author fizzymagic7
 \date   9/6/2003
 
 \return  [std::string] - Longitude in DDD MM.MMM format.
 */
// -------------------------------------------------------------------------
std::string CLatLon::LongToDMM(void)
{
    std::ostringstream Stream;
    int iLongDegrees = (int) floor(fabs(m_Longitude));
    double dLongMinutes = 60. * (fabs(m_Longitude) - floor(fabs(m_Longitude)));
    Stream << std::fixed << std::setprecision(3)
    << ((m_Longitude >= 0.)?'E':'W') << " "
    << iLongDegrees  << " "
    << std::setfill('0') << std::setw(6) << dLongMinutes;
    return Stream.str();
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::LongToDMS()
/*! 
 \brief  Converts longitude coordinate into DDD MM SS.SS format.
 
 \author fizzymagic7
 \date   9/6/2003
 
 \return  [std::string] - Longitude in DDD MM SS.SS format.
 */
// -------------------------------------------------------------------------
std::string CLatLon::LongToDMS(void)
{
    std::ostringstream Stream;
    int iLongDegrees = (int) floor(fabs(m_Longitude));
    int iLongMinutes = (int) floor(60. * (fabs(m_Longitude) - floor(fabs(m_Longitude))) + COORD_EPSILON);
    double dLongSeconds = 60. * (fabs(m_Longitude * 60.) - floor(fabs(m_Longitude * 60.)));
    Stream << std::fixed << std::setprecision(2)
    << ((m_Longitude >= 0.)?'E':'W') << " "
    << iLongDegrees  << " "
    << std::setfill('0') << std::setw(2) << iLongMinutes  << " "
    << std::setfill('0') << std::setw(5) << dLongSeconds;
    return Stream.str();
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ToDDD()
/*! 
 \brief  Converts coordinates into DDD.DDDDDD format.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [std::string] - Output string containing both coordinates in
 DDD.DDDDDD format.
 
 This method converts both latitude and longitude into a single string
 in which both are formatted as decimal degrees.  The string will look
 like:
 
 DD.DDDDDD,-DDD.DDDDDD
 
 with the proper signs in front of both parts. No space is placed in 
 between the coordinates in order to facilitate importing into 
 spreadsheets.
 */
// -------------------------------------------------------------------------
std::string CLatLon::ToDDD(void)
{
    std::ostringstream Stream;
    Stream << LatToDDD() << "," << LongToDDD();
    return Stream.str();
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ToDMM()
/*! 
 \brief  Converts coordinates into DDD MM.MMM format.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [std::string] - Output string containing both coordinates in
 DDD MM.MMM format.
 
 This method converts both latitude and longitude into a single string
 in which both are formatted as degrees and decimal minutes.  The string 
 will look like:
 
 N DD MM.MMM, W DDD MM.MMM
 */
// -------------------------------------------------------------------------
std::string CLatLon::ToDMM(void)
{
    std::ostringstream Stream;
    Stream << LatToDMM() << ", " << LongToDMM();
    return Stream.str();
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ToDMS()
/*! 
 \brief  Converts coordinates into DDD MM SS.SS format.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [std::string] - Output string containing both coordinates in
 DDD MM SS.SS format.
 
 This method converts both latitude and longitude into a single string
 in which both are formatted as degrees, minutes and seconds.  The string 
 will look like:
 
 N DD MM SS.SS, W DDD MM SS.SS
 */
// -------------------------------------------------------------------------
std::string CLatLon::ToDMS(void)
{
    std::ostringstream Stream;
    Stream << LatToDMS() << ", " << LongToDMS();
    return Stream.str();
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ToUTM()
/*! 
 \brief  Converts coordinates into UTM coordinates.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [std::string] - Output string in which the coordinates are
 in standard UTM format.
 
 This method converts the coordinates into UTM coordinates using the
 WGS84 ellipsoid.  The output string will look like:
 
 ZZL E eeeeee N nnnnnnn 
 
 where ZZ is the zone number, L is the zone letter, eeeeee is the easting,
 and nnnnnnn is the northing.
 */
// -------------------------------------------------------------------------
std::string CLatLon::ToUTM(void)
{
    std::ostringstream Stream;
    int iZone = GetZone();
    char cZoneLetter = GetZoneLetter();
    const double k0 = 0.9996;
    const double dEastingOffset = 5.e5;
    const double dNorthingOffsetSouth = 1.e7;
    
    double a = m_Ellipsoid.m_a;
    double flat = 1./m_Ellipsoid.m_fInv;
    double ecc2 = 2. * flat - flat * flat;
    double ecc12 = (ecc2)/(1. - ecc2);
    
    
    double dLatitude = m_Latitude * m_Deg2Rad;
    double dLongitude = (fmod(m_Longitude + 180., 360.) - 180.) * m_Deg2Rad;
    double dCentralLongitude = ((iZone - 1) * 6. - 180. + 3.) * m_Deg2Rad;
    
    double N, T, C, A, M;
    double dUTMNorthing, dUTMEasting;
    
    N = a / sqrt(1. - ecc2*sin(dLatitude)*sin(dLatitude));
    T = tan(dLatitude)*tan(dLatitude);
    C = ecc12*cos(dLatitude)*cos(dLatitude);
    A = cos(dLatitude)*(dLongitude - dCentralLongitude);
    
    M = a * ((1. - ecc2/4. - 3.*ecc2*ecc2/64. - 5.*ecc2*ecc2*ecc2/256.)*dLatitude 
             - (3.*ecc2/8. + 3.*ecc2*ecc2/32. + 45.*ecc2*ecc2*ecc2/1024.)*sin(2.*dLatitude)
             + (15.*ecc2*ecc2/256. + 45.*ecc2*ecc2*ecc2/1024.)*sin(4.*dLatitude) 
             - (35.*ecc2*ecc2*ecc2/3072.)*sin(6.*dLatitude));
    
    dUTMEasting = (double)(k0*N*(A+(1. - T + C)*A*A*A/6.
                                 + (5. - 18.*T+T*T + 72.*C - 58.*ecc12)*A*A*A*A*A/120.)
                           + dEastingOffset);
    
    dUTMNorthing = (double)(k0*(M + N*tan(dLatitude)*(A*A/2. + (5. - T + 9.*C + 4.*C*C)*A*A*A*A/24.
                                                      + (61. - 58.*T+T*T + 600.*C - 330.*ecc12)*A*A*A*A*A*A/720.)));
    
    if (m_Latitude < 0.) {
        dUTMNorthing += dNorthingOffsetSouth;
    }
    
    dUTMEasting = floor(dUTMEasting + 0.5);
    dUTMNorthing = floor(dUTMNorthing + 0.5);
    
    Stream << std::fixed << std::setprecision(0)
    << iZone << cZoneLetter << " "
    << "E " << dUTMEasting << " " 
    << "N " << dUTMNorthing;
    return Stream.str();
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::VincentyDistance()
/*! 
 \brief  Calculates the distance between this point and P using the Vincenty
 method.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [double] - Distance between this point and P in meters.
 
 \param P [CLatLon&] - Point to which to compute distance.
 
 This method computes a high-accuracy distance between this point and P
 using the Vincenty method and the WGS84 ellipsoid.
 */
// -------------------------------------------------------------------------
double CLatLon::VincentyDistance(CLatLon& P)
{
    double Tmp1, Tmp2;
    return this->VincentyDistance(P, &Tmp1, &Tmp2);
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::VincentyDistance()
/*! 
 \brief  Calculates the distance and forward and reverse azimuths between 
 this point and P using the Vincenty method.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [double] - Distance between this point and P in meters.
 
 \param P [CLatLon&] - Point to which to compute distance.
 \param pForwardAzimuth [double *] - Pointer to parameter to receive 
 forward azimuth in degrees.
 \param pReverseAzimuth [double *] - Pointer to parameter to receive
 reverse azimuth in degrees.
 
 
 This method computes a high-accuracy distance between this point and P
 using the Vincenty method and the WGS84 ellipsoid.  It also calculates 
 and returns the forward and reverse azimuths.
 */
// -------------------------------------------------------------------------
double CLatLon::VincentyDistance(CLatLon& P, double *pForwardAzimuth, double *pReverseAzimuth)
{
    if (m_Latitude == P.m_Latitude
        && m_Longitude == P.m_Longitude) {
        *pForwardAzimuth = 0.;
        *pReverseAzimuth = 0.;
        return 0.;
    }
    
    // Check to see if either latitude is 90 degrees exactly.
    // In that case, the distance is a closed-form expression!
    
    double dLat1 = m_Deg2Rad * m_Latitude;
    double dLat2 = m_Deg2Rad * P.m_Latitude;
    double dLong1 = m_Deg2Rad * m_Longitude;
    double dLong2 = m_Deg2Rad * P.m_Longitude;
    
    
    double a0 = m_Ellipsoid.m_a;
    double flat = 1./m_Ellipsoid.m_fInv;
    double r = 1. - flat;
    double b0 = a0 * r;
    
    double tanu1 = r * tan(dLat1);
    double tanu2 = r * tan(dLat2);
    
    double dtmp;
    dtmp = atan(tanu1);
    if (fabs(m_Latitude) >= 90.) dtmp = dLat1;
    double cosu1 = cos(dtmp);
    double sinu1 = sin(dtmp);
    
    dtmp = atan(tanu2);
    if (fabs(P.m_Latitude) >= 90.) dtmp = dLat2;
    double cosu2 = cos(dtmp);
    double sinu2 = sin(dtmp);
    
    double omega = dLong2 - dLong1;
    
    double lambda = omega;
    
    double testlambda, ss1, ss2, ss, cs, sigma, sinalpha, cosalpha2, c2sm, c, dDeltaLambda;
    
    do {
        testlambda = lambda;
        ss1 = cosu2 * sin(lambda);
        ss2 = cosu1 * sinu2 - sinu1 * cosu2 * cos(lambda);
        ss = sqrt(ss1*ss1 + ss2 * ss2);
        cs = sinu1 * sinu2 + cosu1 * cosu2 * cos(lambda);
        sigma = atan2(ss, cs);
        sinalpha = cosu1 * cosu2 * sin(lambda) / ss;
        cosalpha2 = 1. - sinalpha * sinalpha; 
        c2sm = cs - 2.*sinu1*sinu2/cosalpha2;
        c = flat/16. * cosalpha2*(4. + flat*(4. - 3.*cosalpha2));
        lambda = omega + (1. - c)*flat*sinalpha*(sigma + c*ss*(c2sm + c*cs*(-1. + 2.*c2sm*c2sm)));
        dDeltaLambda = fabs(testlambda - lambda);
    } while (dDeltaLambda > EPSILON);
    
    double u2 = cosalpha2 * (a0*a0 - b0*b0)/(b0*b0);
    double a = 1. + (u2 / 16384.) * (4096. + u2 * (-768. + u2 * (320. - 175. * u2)));
    double b = (u2 / 1024.) * (256. + u2 * (-128. + u2 * (74. - 47. * u2)));
    
    double dsigma = b * ss * (c2sm + (b / 4.) * (cs * (-1. + 2. * c2sm*c2sm) 
                                                 - (b / 6.) * c2sm * (-3. + 4. * ss*ss) * (-3. + 4. * c2sm*c2sm)));
    
    double s = b0 * a * (sigma - dsigma);
    
    double alpha12 = atan2(cosu2 * sin(lambda), (cosu1 * sinu2 - sinu1 * cosu2 * cos(lambda)))/m_Deg2Rad;
    double alpha21 = atan2(cosu1 * sin(lambda), (-sinu1 * cosu2 + cosu1 * sinu2 * cos(lambda)))/m_Deg2Rad;
    
    *pForwardAzimuth = fmod(alpha12 + 360., 360.);
    *pReverseAzimuth = fmod(alpha21 + 180., 360.);
    
    return s;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::VincentyProjection()
/*! 
 \brief  Computes a projection along a geodesic using the Vincenty method.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [CLatLon] - Projected coordinates.
 
 \param dAzimuth [double] - Forward azimuth in degrees.
 \param dDistance [double] - Distance in meters.
 */
// -------------------------------------------------------------------------
CLatLon CLatLon::VincentyProjection(double dAzimuth, double dDistance)
{
    CLatLon RetCoord;
    
    double lastsigma, twosigmam, ss, cs, c2sm, deltasigma;
    double term1, term2, term3, term4;
    
    dAzimuth *= m_Deg2Rad;
    double dLat1 = m_Deg2Rad * m_Latitude;
    double dLong1 = m_Deg2Rad * m_Longitude;
    double s = dDistance;
    
    double a0 = m_Ellipsoid.m_a;
    double flat = 1./m_Ellipsoid.m_fInv;
    double r = 1. - flat;
    double b0 = a0 * r;
    double tanu1 = r * tan(dLat1);
    
    double tansigma1 = tanu1 / cos(dAzimuth);
    double u1 = atan(tanu1);
    double sinu1 = sin(u1);
    double cosu1 = cos(u1);
    
    double sinalpha = cosu1 * sin(dAzimuth);
    double cosalpha = sqrt(1. - sinalpha*sinalpha);
    
    double usqr = cosalpha*cosalpha * (a0*a0 - b0*b0) / (b0*b0);
    
    term1 = usqr / 16384.;
    term2 = 4096. + usqr * (-768. + usqr * (320. - 175. * usqr));
    double a = 1. + term1 * term2;
    double b = usqr / 1024. * (256. + usqr * (-128. + usqr * (74. - 47. * usqr)));
    
    double sigma = s / (b0 * a);
    double sigma1 = atan(tansigma1);
    
    do {
        lastsigma = sigma;
        twosigmam = 2.* sigma1 + sigma;
        ss = sin(sigma);
        cs = cos(sigma);
        c2sm = cos(twosigmam);
        
        deltasigma = b * ss * (c2sm + b / 4. * (cs * (-1. + 2. * c2sm*c2sm) 
                                                - b / 6. * c2sm * (-3. + 4. * ss*ss) * (-3. + 4. * c2sm*c2sm)));
        
        sigma = s / (b0 * a) + deltasigma;
        
    } while (fabs(sigma - lastsigma) > EPSILON);
    
    twosigmam = 2. * sigma1 + sigma;
    ss = sin(sigma);
    cs = cos(sigma);
    c2sm = cos(twosigmam);
    term1 = sinu1 * cs + cosu1 * ss * cos(dAzimuth);
    term4 = sinu1 * ss - cosu1 * cs * cos(dAzimuth);
    term2 = sinalpha*sinalpha + term4*term4;
    term3 = r * sqrt(term2);
    
    double dLat2 = atan2(term1, term3);
    
    term1 = ss * sin(dAzimuth);
    term2 = cosu1 * cs - sinu1 * ss * cos(dAzimuth);
    //double tanlambda = term1 / term2;
    double lambda = atan2(term1, term2);
    
    double c = flat / 16. * cosalpha*cosalpha * (4. + flat * (4. - 3. * cosalpha*cosalpha));
    
    double omega = lambda - (1. - c) * flat * sinalpha * (sigma + c * ss * (c2sm + c * cs * (-1. + 2. * c2sm*c2sm)));
    
    double dLong2 = dLong1 + omega;
    
   // term1 = -sinu1 * ss + cosu1 * cs * cos(dAzimuth);
    
    RetCoord.m_Latitude = NormalizeLatitude(dLat2 / m_Deg2Rad);
    RetCoord.m_Longitude = NormalizeLongitude(dLong2 / m_Deg2Rad);
    
    return RetCoord;
}

int CLatLon::GetChecksum(std::string& strCoordString) {
    int iSum = 0;
    std::string::iterator sit;
    for (sit = strCoordString.begin(); sit != strCoordString.end(); sit++) {
        if (::isdigit(*sit)) {
            iSum += (*sit - '0');
        }
    }
    return iSum;
}

int CLatLon::DigitalRoot(int n) {
    return 1 + ((n-1) % 9);
}


// -------------------------------------------------------------------------
// METHOD:  CLatLon::CleanCoordString()
/*! 
 \brief  Cleans coordinate strings.
 
 \author fizzymagic
 \date   9/6/2003
 
 \param strCoordString [std::string&] - Cleaned string.
 */
// -------------------------------------------------------------------------
void CLatLon::CleanCoordString(std::string& strCoordString)
{
    std::transform(strCoordString.begin(), strCoordString.end(), strCoordString.begin(), ::toupper);
    findandreplace(strCoordString, "LATITUDE", "");
    findandreplace(strCoordString, "LONGITUDE", "");
    findandreplace(strCoordString, "LAT", "");
    findandreplace(strCoordString, "LON", "");
    findandreplace(strCoordString, "NORTH", "N");
    findandreplace(strCoordString, "SOUTH", "S");
    findandreplace(strCoordString, "EAST", "E");
    findandreplace(strCoordString, "WEST", "W");
    std::string::iterator it;
    for (it = strCoordString.begin(); it != strCoordString.end(); it++) {
        if (!((*it >= 'A' && *it <= 'Z') || (*it >= '0' && *it <= '9') || *it == '-' || *it == '.' || *it == ',')) {
            *it = ' ';
        }
    }
    trim(strCoordString);
    findandreplace(strCoordString, "  ", " ");
}

bool CLatLon::SplitCoordString(std::string& strCoordString, std::string& strLatString, std::string& strLongString)
{
    std::string::size_type iSplitPoint = 0;
    int iLength = (int)strCoordString.size();
    
    CleanCoordString(strCoordString);
    
    // First, check to see if there is a comma anywhere in the string.
    if ((iSplitPoint = strCoordString.find(",")) == std::string::npos) {
        bool f_Digits = false;
        int i;
        for (i=0; i<iLength; i++) {
            char c = strCoordString[i];
            if (::isspace(c)) continue;
            else if (!f_Digits && (::isdigit(c) || c == '+' || c == '-')) {
                f_Digits = true;
            }
            else if (c == 'N' || c == 'S') {
                if (f_Digits) {
                    iSplitPoint = i+1;
                    break;
                }
            }
            else if (c == 'E' || c == 'W') {
                iSplitPoint = i-1;
                break;
            }
        }
    }
    strLatString.assign(strCoordString, 0, iSplitPoint);
    strLongString.assign(strCoordString, iSplitPoint+1, iLength - (iSplitPoint+1));
    findandreplace(strLatString, ",", "");
    findandreplace(strLongString, ",", "");
    if (strLongString[0] == ',') strLongString[0] = ' ';
    
    return true;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::GetZone()
/*! 
 \brief  Gets UTM zone of a point.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [int] - UTM zone for this point.
 */
// -------------------------------------------------------------------------
int CLatLon::GetZone(void)
{
    int iZone = ((int) floor((m_Longitude + 180.)/6.)) % 60  + 1;
    
    // Special zones for places in northern Europe.
    if (m_Latitude > 56. && m_Latitude <= 64. 
        && m_Longitude > 3. && m_Longitude <= 12.) {
        iZone = 32;
    }
    
    if (m_Latitude > 72. && m_Latitude < 84.) {
        if ((m_Longitude >= 0.) && (m_Longitude < 9.)) iZone = 31;
        else if ((m_Longitude >= 9.) && (m_Longitude < 21.)) iZone = 33;
        else if ((m_Longitude >= 21.) && (m_Longitude < 33.)) iZone = 35;
        else if ((m_Longitude >= 33.) && (m_Longitude < 42.)) iZone = 37;
    }
    
    return iZone;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::GetZoneLetter()
/*! 
 \brief  Gets UTM zone letter of a point.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return [char] - UTM zone letter for this point.
 */
// -------------------------------------------------------------------------
char CLatLon::GetZoneLetter(void)
{
    char cZoneLetter;
    if (m_Latitude >= 72. && m_Latitude <= 84.) {
        cZoneLetter = 'X';
    }
    else {
        int iZoneLetter = (int) floor(m_Latitude + 80.) / 8;
        if (iZoneLetter >= 0 && iZoneLetter < m_strZoneLetters.size()) {
            cZoneLetter = m_strZoneLetters[iZoneLetter];
        }
        else {
            cZoneLetter = 'Z';
        }
    }
    return cZoneLetter;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ConvertUTM()
/*! 
 \brief  Converts UTM Zone, Zone Letter, Easting, and Northing to a 
 CLatLon point.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [bool] - true if conversion was performed without error.
 
 \param iZone [int] - UTM Zone
 \param iZoneLetter [char] - UTM Zone Letter
 \param dEasting [double] - UTM Easting
 \param dNorthing [double] - UTM Northing
 
 */
// -------------------------------------------------------------------------
bool CLatLon::ConvertUTM(int iZone, char cZoneLetter, double dEasting, double dNorthing)
{
    const double k0 = 0.9996;
    const double dEastingOffset = 5.e5;
    const double dNorthingOffsetSouth = 1.e7;
    
    double dLatitude, dLongitude;
    
    double a = m_Ellipsoid.m_a;
    double flat = 1./m_Ellipsoid.m_fInv;
    double ecc2 = 2. * flat - flat * flat;
    double ecc12 = (ecc2)/(1. - ecc2);
    double e1 = (1. - sqrt(1. - ecc2))/(1. + sqrt(1. - ecc2));
    double N1, T1, C1, R1, D, M;
    double dCentralLongitude;
    double mu, phi1;
    double x, y;
    
    x = dEasting - dEastingOffset;
    y = dNorthing;
    cZoneLetter = toupper(cZoneLetter);
    if (m_strZoneLetters.find(cZoneLetter) == std::string::npos) return false;
    
    if ((cZoneLetter - 'N') < 0) {
        y -= dNorthingOffsetSouth;
    }
    
    dCentralLongitude = ((double) (iZone - 1) * 6. - 180. + 3.) * m_Deg2Rad;
    
    M = y / k0;
    mu = M / (a * (1. - ecc2/4. - 3.*ecc2*ecc2/64. - 5.*ecc2*ecc2*ecc2/256.));
    
    phi1 = mu + (3.*e1/2. - 27.*e1*e1*e1/32.) * sin(2.*mu) 
    + (21.*e1*e1/16. - 55.*e1*e1*e1*e1/32.) * sin(4.*mu)
    + (151.*e1*e1*e1/96.) * sin(6.*mu);
    
    N1 = a / sqrt(1. - ecc2 * sin(phi1)*sin(phi1));
    T1 = tan(phi1)*tan(phi1);
    C1 = ecc12 * cos(phi1)*cos(phi1);
    R1 = a * (1. - ecc2) / pow(1. - ecc2*sin(phi1)*sin(phi1), 1.5);
    D = x/(N1*k0);
    
    dLatitude = phi1 - (N1*tan(phi1)/R1) * (D*D/2. - (5. + 3.*T1 + 10.*C1 - 4.*C1*C1 - 9.*ecc12)*D*D*D*D/24.
                                            + (61. + 90.*T1 + 298.*C1 + 45.*T1*T1 - 252.*ecc12 - 3.*C1*C1)*D*D*D*D*D*D/720.);
    
    dLongitude = (D - (1.+ 2.*T1 + C1)*D*D*D/6. +(5.- 2.*C1 + 28.*T1 - 3.*C1*C1 + 8.*ecc12 + 24.*T1*T1)*D*D*D*D*D/120.)/cos(phi1);
    
    m_Latitude = dLatitude / m_Deg2Rad;
    m_Longitude = (dCentralLongitude + dLongitude) / m_Deg2Rad;
    
    return true;
}

double CLatLon::GeoSatelliteAzEl(double dSatLongitude, double dAltitude, double *pAzimuth, double *pElevation)
{
    CLatLon Sat(0., dSatLongitude);
    CCartesianCoord GroundCartesian(*this, dAltitude);
    CCartesianCoord SatCartesian(Sat, m_GeosyncRadius - m_Radius);
    CCartesianCoord DiffCartesian = GroundCartesian - SatCartesian;
    double Distance = DiffCartesian.Norm();
    
    double dLat = m_Deg2Rad * m_Latitude;
    double dDeltaLong = m_Deg2Rad * (dSatLongitude - m_Longitude);
    
    double Ra = sqrt(GroundCartesian.m_x*GroundCartesian.m_x + GroundCartesian.m_y*GroundCartesian.m_y);
    double Rz = GroundCartesian.m_z;
    
    double DeltaRx = m_GeosyncRadius * cos(dDeltaLong) - Ra;
    double DeltaRy = m_GeosyncRadius * sin(dDeltaLong);
    double DeltaRz = -Rz;
    
    double DeltaRNorth = DeltaRz * cos(dLat) - DeltaRx * sin(dLat);
    double DeltaRZenith = DeltaRx * cos(dLat) + DeltaRz * sin(dLat);
    
    // Use the spherical azimuth here because the Vincenty azimuth won't point right at the satellite.
    double Azimuth = atan2(sin(dDeltaLong), -sin(dLat)*cos(dDeltaLong)) / m_Deg2Rad;
    if (Azimuth < 0) Azimuth += 360.;
    double Elevation = atan(DeltaRZenith / sqrt(DeltaRNorth*DeltaRNorth + DeltaRy*DeltaRy)) / m_Deg2Rad;
    
    *pAzimuth = fmod(Azimuth, 360.);
    *pElevation = Elevation;
    
    return Distance;
}

double CLatLon::GeoSatelliteAzElSpherical(double dSatLongitude, double dAltitude, double *pAzimuth, double *pElevation)
{
    CLatLon Sat(0., dSatLongitude);
    
    double dLat = m_Deg2Rad * m_Latitude;
    double dDeltaLong = m_Deg2Rad * (dSatLongitude - m_Longitude);
    
    double R = m_Radius + dAltitude;
    double b = cos(dLat) * cos(dDeltaLong);
    double Azimuth = atan2(sin(dDeltaLong), -sin(dLat)*cos(dDeltaLong)) / m_Deg2Rad;
    if (Azimuth < 0) Azimuth += 360.;
    double Distance = sqrt(R*R + m_GeosyncRadius*m_GeosyncRadius - 2.*R*m_GeosyncRadius*b);
    double Elevation = atan((b - R/m_GeosyncRadius) / sqrt(1. - b*b)) / m_Deg2Rad;
    
    *pAzimuth = fmod(Azimuth, 360.);
    *pElevation = Elevation;
    
    return Distance;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ToCartesian()
/*! 
 \brief  Converts lat/lon into Cartesian coordinates using ellipsoid.
 
 \author fizzymagic
 \date   9/7/2003
 
 \return  [CCartesianCoord] - Cartesian coordinates in meters.
 */
// -------------------------------------------------------------------------
CCartesianCoord CLatLon::ToCartesian(double dElevation)
{
    double a = m_Ellipsoid.m_a;
    double flat = 1./m_Ellipsoid.m_fInv;
    double ecc2 = 2. * flat - flat * flat;
    double sinLat = sin(m_Latitude * m_Deg2Rad);
    double cosLat = cos(m_Latitude * m_Deg2Rad);
    double w = sqrt(1. - ecc2 * sinLat * sinLat);
    double r = a / w;
    
    double dLong = m_Longitude * m_Deg2Rad;
    double x = (r + dElevation) * cos(dLong) * cosLat;
    double y = (r + dElevation) * sin(dLong) * cosLat;
    double z = (r * (1. - ecc2) + dElevation) * sinLat;
    
    return CCartesianCoord(x, y, z);
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::FromCartesian()
/*! 
 \brief  Converts Cartesian coordinates into lat/lon using ellipsoid.
 
 \author fizzymagic
 \date   9/7/2003
 
 \param C [CCartesianCoord&] - Cartesian coordinates of point.
 */
// -------------------------------------------------------------------------
double CLatLon::FromCartesian(CCartesianCoord& C)
{
    double a = m_Ellipsoid.m_a;
    double flat = 1./m_Ellipsoid.m_fInv;
    double ecc2 = 2. * flat - flat * flat;
    double r = a * ecc2;
    
    double p = sqrt(C.m_x*C.m_x + C.m_y*C.m_y);
    
    double dTmp = C.m_z / (p * (1. - ecc2));
    double dLast;
    
    do {
        dLast = dTmp;
        dTmp = C.m_z / (p - r / sqrt(1. + (1. - ecc2) * dTmp * dTmp));
    } while (fabs(dLast - dTmp) > EPSILON);
    
    m_Latitude = atan(dTmp) / m_Deg2Rad;
    m_Longitude = atan2(C.m_y, C.m_x) / m_Deg2Rad;
    
    double cosLat = cos(m_Latitude * m_Deg2Rad);
    double dElevation = p / cosLat - a/sqrt(1. - ecc2 * (1. - cosLat * cosLat));
    return dElevation;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::ToSphericalCartesian()
/*! 
 \brief  Converts lat/lon into Cartesian coordinates using spheroid.
 
 \author fizzymagic
 \date   9/7/2003
 
 \return  [CCartesianCoord] - Cartesian coordinates in meters.
 */
// -------------------------------------------------------------------------
CCartesianCoord CLatLon::ToSphericalCartesian(double dElevation)
{
    double dLat = m_Latitude * m_Deg2Rad;
    double dLong = m_Longitude * m_Deg2Rad;
    return CCartesianCoord((m_Radius  + dElevation) * cos(dLong) * cos(dLat),
                           (m_Radius  + dElevation) * sin(dLong) * cos(dLat),
                           (m_Radius  + dElevation) * sin(dLat));
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::FromSphericalCartesian()
/*! 
 \brief  Converts Cartesian coordinates into lat/lon using spheroid.
 
 \author fizzymagic
 \date   9/7/2003
 
 \param C [CCartesianCoord&] - Cartesian coordinates of point.
 */
// -------------------------------------------------------------------------
void CLatLon::FromSphericalCartesian(CCartesianCoord& C)
{
    m_Latitude = atan2(C.m_z, sqrt(C.m_x*C.m_x + C.m_y*C.m_y)) / m_Deg2Rad;
    m_Longitude = atan2(C.m_y, C.m_x) / m_Deg2Rad;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::Molodensky()
/*! 
 \brief Performs a Molodensky transformation between datums
 
 \author fizzymagic
 \date   2/20/2010
 
 \return  [CLatLon] - Coordinates in the To datum
 
 \param From [CEllipsoid&] - Datum to transform from
 \param To [CEllipsoid&] - Datum to transform to
 
 */
// -------------------------------------------------------------------------
CLatLon CLatLon::Molodensky(CEllipsoid& From, CEllipsoid& To)
{
    CLatLon RetCoord;
    double dLat = m_Latitude * m_Deg2Rad;
    double dLong = m_Longitude * m_Deg2Rad;
    double sinLat = sin(dLat);
    double cosLat = cos(dLat);
    double sinLon = sin(dLong);
    double cosLon = cos(dLong);
    double sinLat2 = sinLat * sinLat;
    double a = From.m_a;
    double flat = 1./From.m_fInv;
    double ecc2 = 1. - (1. - flat)*(1. - flat);
    double ab = 1. / (1. - flat);
    double rn = a / sqrt(1.0 - ecc2 * sinLat2);
    double rm = a * (1. - ecc2) / pow((1.0 - ecc2*sinLat2), 1.5);
    double dx = To.m_dx - From.m_dx;
    double dy = To.m_dy - From.m_dy;
    double dz = To.m_dz - From.m_dz;
    double df = 1./To.m_fInv - 1./From.m_fInv;
    double da = To.m_a - From.m_a;
    double DeltaLat, DeltaLon;
    DeltaLat = (((((dx*sinLat*cosLon + dy*sinLat*sinLon) - dz*cosLat)
                  + (da*((rn*ecc2*sinLat*cosLat)/a)))
                 + (df*(rm*ab + rn/ab)*sinLat*cosLat))) / rm;
    DeltaLon = (dx*sinLon - dy*cosLon) / (rn*cosLat);
    
    RetCoord.m_Latitude = NormalizeLatitude((dLat + DeltaLat) / m_Deg2Rad);
    RetCoord.m_Longitude = NormalizeLongitude((dLong + DeltaLon) / m_Deg2Rad);
    
    return RetCoord;
}

// -------------------------------------------------------------------------
// METHOD:  CLatLon::SphericalCross()
/*! 
 \brief  Performs a cross-product between this and another CLatLon point
 using spheroidal approximation.
 
 \author fizzymagic
 \date   9/6/2003
 
 \return  [CCartesianCoord] - Cross-product vector (un-normalized).
 
 \param P [CLatLon&] - Point with which to perform cross-product.
 */
// -------------------------------------------------------------------------
CCartesianCoord CLatLon::SphericalCross(CLatLon& P)
{
    double dLat1 = m_Latitude * m_Deg2Rad;
    double dLong1 = m_Longitude * m_Deg2Rad;
    double dLat2 = P.m_Latitude * m_Deg2Rad;
    double dLong2 = P.m_Longitude * m_Deg2Rad;
    double dDeltaLat = dLat1 - dLat2;
    double dSumLat = dLat1 + dLat2;
    double dDeltaLong = (dLong1 - dLong2)/2.;
    double dAvgLong = (dLong1 + dLong2)/2.;
    
    return CCartesianCoord(sin(dSumLat) * cos(dAvgLong) * sin(dDeltaLong) 
                           - sin(dDeltaLat) * sin(dAvgLong) * cos(dDeltaLong),
                           sin(dDeltaLat) * cos(dAvgLong) * cos(dDeltaLong) 
                           + sin(dSumLat) * sin(dAvgLong) * sin(dDeltaLong),
                           cos(dLat1) * cos(dLat2) * sin(-2.*dDeltaLong));
}


void findandreplace(std::string& strSource, const char *szFind, const char *szReplace)
{
    std::string::size_type i;
    while ((i = strSource.find(szFind)) != std::string::npos) {
        strSource.replace(i, strlen(szFind), szReplace);
    }
}

void trim(std::string& strSource)
{
    std::reverse(strSource.begin(), strSource.end());
    std::string::iterator it = strSource.begin(), eit = strSource.end();
    while (it != eit && ::isspace(*it)) it++;
    strSource.erase(strSource.begin(), it);
    std::reverse(strSource.begin(), strSource.end());
    it = strSource.begin(), eit = strSource.end();
    while (it != eit && ::isspace(*it)) it++;
    strSource.erase(strSource.begin(), it);
}

