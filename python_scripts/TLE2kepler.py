import math
from datetime import datetime, timedelta

# Earth's gravitational constant (m^3/s^2)
MU = 3.986004418e14
# Earth's radius (km)
EARTH_RADIUS = 6378.137

def tle_to_keplerian(tle_line1, tle_line2):
    """
    Convert TLE (Two-Line Element) to Keplerian orbital elements.
    
    :param tle_line1: First line of the TLE
    :param tle_line2: Second line of the TLE
    :return: Dictionary containing Keplerian elements
    """
    # Parse TLE
    satellite_number = int(tle_line1[2:7])
    classification = tle_line1[7]
    
    # Parse epoch
    epoch_year = int(tle_line1[18:20])
    epoch_day = float(tle_line1[20:32])
    epoch = datetime(2000 + epoch_year, 1, 1) + timedelta(days=epoch_day - 1)
    
    inclination = float(tle_line2[8:16])
    raan = float(tle_line2[17:25])
    eccentricity = float("0." + tle_line2[26:33])
    argument_of_perigee = float(tle_line2[34:42])
    mean_anomaly = float(tle_line2[43:51])
    mean_motion = float(tle_line2[52:63])
    
    # Calculate semi-major axis
    n = mean_motion * 2 * math.pi / 86400  # Convert to radians per second
    a = (MU / (n * n)) ** (1/3) / 1000  # Convert to km
    
    # Calculate perigee and apogee
    perigee = a * (1 - eccentricity) - EARTH_RADIUS
    apogee = a * (1 + eccentricity) - EARTH_RADIUS
    
    # Calculate period
    period = 2 * math.pi / n
    
    return {
        "satellite_number": satellite_number,
        "classification": classification,
        "epoch": epoch,
        "inclination": inclination,
        "raan": raan,
        "eccentricity": eccentricity,
        "argument_of_perigee": argument_of_perigee,
        "mean_anomaly": mean_anomaly,
        "mean_motion": mean_motion,
        "semi_major_axis": a,
        "perigee": perigee,
        "apogee": apogee,
        "period": period
    }

# Example usage
if __name__ == "__main__":
    # Example TLE for ISS (ZARYA)
    tle_line1 = "1 25544U 98067A   23105.50000000  .00000000  00000-0  00000-0 0  9995"
    tle_line2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.50409622345678"
    
    keplerian_elements = tle_to_keplerian(tle_line1, tle_line2)
    
    print("Keplerian Orbital Elements:")
    for key, value in keplerian_elements.items():
        if isinstance(value, datetime):
            print(f"{key}: {value.strftime('%Y-%m-%d %H:%M:%S')}")
        elif isinstance(value, float):
            print(f"{key}: {value:.6f}")
        else:
            print(f"{key}: {value}")