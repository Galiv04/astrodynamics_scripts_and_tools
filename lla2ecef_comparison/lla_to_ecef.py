import sys
import json
from pyproj import Proj, transform

def lla_to_ecef(lat, lon, alt):
    # Define the projection for WGS84
    wgs84 = Proj(proj='latlong', datum='WGS84')
    
    # Define the projection for ECEF
    ecef = Proj(proj='geocent', datum='WGS84')
    
    # Transform coordinates from LLA to ECEF
    x, y, z = transform(wgs84, ecef, lon, lat, alt)
    
    return x, y, z

if len(sys.argv) != 4:
    print("Usage: python3 lla_to_ecef.py <latitude> <longitude> <altitude>")
    sys.exit(1)

lat = float(sys.argv[1])
lon = float(sys.argv[2])
alt = float(sys.argv[3])

x, y, z = lla_to_ecef(lat, lon, alt)

print(json.dumps({"x": x, "y": y, "z": z}, indent=2))