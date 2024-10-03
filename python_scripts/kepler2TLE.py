import math
from datetime import datetime, timedelta

# Gravitational constant * Earth's mass (m^3/s^2)
GM = 3.986004418e14

def calculate_mean_motion(semi_major_axis):
    """
    Calculate mean motion based on the semi-major axis.
    
    :param semi_major_axis: Semi-major axis in kilometers
    :return: Mean motion in revolutions per day
    """
    # Convert semi-major axis to meters
    semi_major_axis *= 1000
    
    # Calculate orbital period in seconds
    period = 2 * math.pi * math.sqrt(semi_major_axis**3 / GM)
    
    # Convert period to mean motion in revolutions per day
    mean_motion = 86400 / period
    
    return mean_motion

def keplerian_to_tle(satellite_name, satellite_number, classification, epoch, inclination, raan, eccentricity, 
                     argument_of_perigee, mean_anomaly, mean_motion, revolution_number):
    """
    Convert Keplerian elements to TLE format.
    
    :param satellite_name: Name of the satellite
    :param satellite_number: NORAD catalog number
    :param classification: Classification (U, C, or S)
    :param epoch: Epoch time as a datetime object
    :param inclination: Inclination (degrees)
    :param raan: Right Ascension of the Ascending Node (degrees)
    :param eccentricity: Eccentricity
    :param argument_of_perigee: Argument of perigee (degrees)
    :param mean_anomaly: Mean anomaly (degrees)
    :param mean_motion: Mean motion (revolutions per day)
    :param revolution_number: Revolution number at epoch
    :return: TLE as a list of two strings
    """
    
    # Format epoch
    epoch_year = epoch.year % 100
    epoch_day = epoch.timetuple().tm_yday + epoch.hour/24.0 + epoch.minute/1440.0 + epoch.second/86400.0
    epoch_str = f"{epoch_year:02d}{epoch_day:012.8f}"
    
    # Format eccentricity
    eccentricity_str = f"{eccentricity:.7f}"[2:]  # Remove leading "0."
    
    # Prepare TLE lines
    line1 = f"1 {satellite_number:5d}{classification}  00000A  {epoch_str}  .00000000  00000-0  00000-0 0  9999"
    line2 = (f"2 {satellite_number:5d} {inclination:8.4f} {raan:8.4f} {eccentricity_str} "
             f"{argument_of_perigee:8.4f} {mean_anomaly:8.4f} {mean_motion:11.8f}{revolution_number:5d}")
    
    # Calculate checksums
    line1 = line1[:-1] + str(calculate_checksum(line1))
    line2 = line2[:-1] + str(calculate_checksum(line2))
    
    return [line1, line2]

def calculate_checksum(line):
    """Calculate the checksum for a TLE line."""
    return sum(int(c) if c.isdigit() else c == '-' for c in line[0:68]) % 10

# Example usage
if __name__ == "__main__":
    # Example Keplerian elements (these are approximate values for the ISS)
    satellite_name = "ISS (ZARYA)"
    satellite_number = 25544
    classification = "U"
    epoch = datetime(2023, 4, 15, 12, 0, 0)
    semi_major_axis = 6884.390540419565  # km
    inclination = 0.06872727105993418
    raan = 49.26928835280404
    eccentricity = 0.0013046051252213477
    argument_of_perigee = 22.404958451694593
    mean_anomaly = 22.466097815775502
    revolution_number = 28891
    
     # Print calculated mean motion
    mean_motion = calculate_mean_motion(semi_major_axis)
    print(f"\nCalculated Mean Motion: {mean_motion:.8f} revolutions per day")
    
    tle = keplerian_to_tle(satellite_name, satellite_number, classification, epoch, inclination, raan, 
                           eccentricity, argument_of_perigee, mean_anomaly, mean_motion, revolution_number)
    
    print(satellite_name)
    print(tle[0])
    print(tle[1])