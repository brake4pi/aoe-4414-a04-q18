# eci_to_ecef.py
#
# Usage: python3 eci_to_ecef.py year month day hour minute second eci_x_km eci_y_km eci_z_km
# 
# Parameters:
#  year: what year is it (gregorian)
#  month: what month (1-12)
#  day: what day is it
#  hour: what hour is it
#  minute: what minute is it
#  second: what second is it
#  eci_x_km: x coordinate in ECI frame
#  eci_y_km: y coordinate in ECI frame
#  eci_z_km: z coordinate in ECI fram
# Output:
#  ecef_x_km: x coordinate in ECEF frame
#  ecef_y_km: y coordinate in ECEF frame
#  ecef_z_km: z coordinate in ECEF frame
#
# Written by Lee Wallenfang
# Other contributors: None
#

# import Python modules
import math # math module
import sys # argv
import numpy
# "constants"
w = 7.292115*(10**(-5))

# helper functions

## function description
def Julian(Y,M,D,H,Min,S):
    JD = D - 32075 + 1461*(Y+4800+(M-14)/12)/4 + 367*(M-2-(M-14)/12*12)/12 - 3*((Y+4900+(M-14)/12)/100)/4
    JDmid = JD-0.5
    DFrac = (S+60*(Min+60*H))/86400
    return JDmid+DFrac

def T_UT1(JulianDate):
    return (JulianDate-2451545.0)/36525

def Theta_GMST_sec(TUT1):
    return 67310.54841 + (876600*60*60+8640184.812866)*TUT1 + 0.093104*TUT1**2 + (-6.2*10**(-6))*TUT1**3
# initialize script arguments
year = ''
month = ''
day = ''
hour = ''
minute = ''
second = ''
eci_x_km = ''
eci_y_km = ''
eci_z_km = ''

# parse script arguments
if len(sys.argv)==10:
 year = float(sys.argv[1])
 month = float(sys.argv[2])
 day = float(sys.argv[3])
 hour = float(sys.argv[4])
 minute = float(sys.argv[5])
 second = float(sys.argv[6])
 eci_x_km = float(sys.argv[7])
 eci_y_km = float(sys.argv[8])
 eci_z_km = float(sys.argv[9])
else:
  print(\
   'Usage: '\
   'python3 arg1 arg2 ...'\
  )
  exit()

# write script below this line
FractionalJulian = Julian(year,month,day,hour,minute,second)
T_UT_1 = T_UT1(FractionalJulian)
Theta_GMST_second = Theta_GMST_sec(T_UT_1)
Theta_GMST_year = Theta_GMST_second/(86400*365)
ActualDay = 86400+3*60+56.555
Theta_GMST_start = Theta_GMST_year*2*math.pi*(math.pi/180)
Theta_GMST_remain = math.fmod(Theta_GMST_second,ActualDay)/(86400*365)*2*math.pi*math.pi/180
Theta_GMST_rad = 2*math.pi-(Theta_GMST_start - Theta_GMST_remain)
[ecef_x_km, ecef_y_km, ecef_z_km] = numpy.matmul([[math.cos(Theta_GMST_rad),-1*math.sin(Theta_GMST_rad),0],[math.sin(Theta_GMST_rad),math.cos(Theta_GMST_rad),0],[0,0,1]],[eci_x_km,eci_y_km,eci_z_km])
print(ecef_x_km)
print(ecef_y_km)
print(ecef_z_km)