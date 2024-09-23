import skyfield
from skyfield import api

# Use these functions from skyfield
load = api.load
EarthSatellite = api.EarthSatellite
from math import degrees, sqrt

def tle_to_orbital_elements(line1, line2):
    satellite = EarthSatellite(line1, line2)
    
    # Get the position and velocity vectors at the epoch
    t = satellite.epoch
    geocentric = satellite.at(t)
    position = geocentric.position.km
    velocity = geocentric.velocity.km_per_s

    # Calculate orbital elements
    mu = 398600.4418  # Earth's gravitational parameter (km^3/s^2)

    r = sqrt(sum(x*x for x in position))
    v = sqrt(sum(x*x for x in velocity))

    # Specific angular momentum
    h = [position[1]*velocity[2] - position[2]*velocity[1],
         position[2]*velocity[0] - position[0]*velocity[2],
         position[0]*velocity[1] - position[1]*velocity[0]]

    # Node vector
    n = [-h[1], h[0], 0]
    n_mag = sqrt(sum(x*x for x in n))

    # Eccentricity vector
    e = [(v*v - mu/r)*x/mu - sum(p*v for p, v in zip(position, velocity))*vx/mu for x, vx in zip(position, velocity)]
    e_mag = sqrt(sum(x*x for x in e))

    # Orbital elements
    a = 1 / (2/r - v*v/mu)  # Semi-major axis
    i = degrees(abs(sum(position[i]*velocity[i] for i in range(3)) / (r*v)))  # Inclination
    raan = degrees(abs(n[1] / n_mag))  # Right ascension of ascending node
    arg_pe = degrees(sum(n[i]*e[i] for i in range(3)) / (n_mag * e_mag))  # Argument of perigee
    M = degrees(abs(sum(e[i]*position[i] for i in range(3)) / (e_mag * r)))  # Mean anomaly

    return {
        "semi_major_axis": a,
        "eccentricity": e_mag,
        "inclination": i,
        "right_ascension_of_ascending_node": raan,
        "argument_of_perigee": arg_pe,
        "mean_anomaly": M
    }

# TLE data
chief_tle_line1 = "1 43797U 18099Q   23262.12345678 -.00000012  00000-0 -34562-5 0  9993"
chief_tle_line2 = "2 43797  97.7465 239.6012 0001247  76.9465 283.1880 15.21991542261234"
target_tle_line1 = "1 20044U 89067A   23261.50000000  .00000000  00000-0  00000-0 0  9999"
target_tle_line2 = "2 20044   0.0142 335.7784 0001641 216.5810 276.9542  1.00273272123456"

# Convert TLEs to orbital elements
chief_elements = tle_to_orbital_elements(chief_tle_line1, chief_tle_line2)
target_elements = tle_to_orbital_elements(target_tle_line1, target_tle_line2)

# Print results
print("Chief Satellite Orbital Elements:")
for key, value in chief_elements.items():
    print(f"{key}: {value}")

print("\nTarget Satellite Orbital Elements:")
for key, value in target_elements.items():
    print(f"{key}: {value}")