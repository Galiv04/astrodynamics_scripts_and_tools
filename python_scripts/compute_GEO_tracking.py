import numpy as np
import skyfield
from skyfield import api
import csv
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, FFMpegWriter

# Use these functions from skyfield
load = api.load
EarthSatellite = api.EarthSatellite
wgs84 = api.wgs84

# If these are available in your skyfield version, uncomment them:
# quaternion_to_matrix = skyfield.quaternion.quaternion_to_matrix
# matrix_to_quaternion = skyfield.quaternion.matrix_to_quaternion

# If the above are not available, keep the NumPy implementations:
def quaternion_to_matrix(q):
    w, x, y, z = q
    return np.array([
        [1 - 2*y*y - 2*z*z, 2*x*y - 2*z*w, 2*x*z + 2*y*w],
        [2*x*y + 2*z*w, 1 - 2*x*x - 2*z*z, 2*y*z - 2*x*w],
        [2*x*z - 2*y*w, 2*y*z + 2*x*w, 1 - 2*x*x - 2*y*y]
    ])

def matrix_to_quaternion(m):
    tr = np.trace(m)
    if tr > 0:
        S = np.sqrt(tr + 1.0) * 2
        w = 0.25 * S
        x = (m[2, 1] - m[1, 2]) / S
        y = (m[0, 2] - m[2, 0]) / S
        z = (m[1, 0] - m[0, 1]) / S
    elif m[0, 0] > m[1, 1] and m[0, 0] > m[2, 2]:
        S = np.sqrt(1.0 + m[0, 0] - m[1, 1] - m[2, 2]) * 2
        w = (m[2, 1] - m[1, 2]) / S
        x = 0.25 * S
        y = (m[0, 1] + m[1, 0]) / S
        z = (m[0, 2] + m[2, 0]) / S
    elif m[1, 1] > m[2, 2]:
        S = np.sqrt(1.0 + m[1, 1] - m[0, 0] - m[2, 2]) * 2
        w = (m[0, 2] - m[2, 0]) / S
        x = (m[0, 1] + m[1, 0]) / S
        y = 0.25 * S
        z = (m[1, 2] + m[2, 1]) / S
    else:
        S = np.sqrt(1.0 + m[2, 2] - m[0, 0] - m[1, 1]) * 2
        w = (m[1, 0] - m[0, 1]) / S
        x = (m[0, 2] + m[2, 0]) / S
        y = (m[1, 2] + m[2, 1]) / S
        z = 0.25 * S
    return np.array([w, x, y, z])

def tle_to_skyfield(line1, line2):
    return EarthSatellite(line1, line2)

def orbital_elements_to_skyfield(a, e, i, raan, argp, M, epoch):
    # Convert orbital elements to Skyfield EarthSatellite
    # a: semi-major axis (km)
    # e: eccentricity
    # i: inclination (degrees)
    # raan: right ascension of the ascending node (degrees)
    # argp: argument of perigee (degrees)
    # M: mean anomaly (degrees)
    # epoch: epoch time as a Time object from Skyfield

    # Convert to radians
    i_rad = np.radians(i)
    raan_rad = np.radians(raan)
    argp_rad = np.radians(argp)
    M_rad = np.radians(M)

    # Calculate position and velocity vectors
    mu = 398600.4418  # Earth's gravitational parameter (km^3/s^2)
    n = np.sqrt(mu / a**3)  # Mean motion

    E = M_rad  # Initial guess for eccentric anomaly
    for _ in range(10):  # Newton's method to solve Kepler's equation
        E = E - (E - e * np.sin(E) - M_rad) / (1 - e * np.cos(E))

    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E/2), np.sqrt(1 - e) * np.cos(E/2))  # True anomaly

    r = a * (1 - e * np.cos(E))  # Distance from central body
    h = np.sqrt(mu * a * (1 - e**2))  # Specific angular momentum

    # Position and velocity in orbital plane
    r_orbital = np.array([r * np.cos(nu), r * np.sin(nu), 0])
    v_orbital = np.array([-mu * np.sin(nu) / h, mu * (e + np.cos(nu)) / h, 0])

    # Rotation matrices
    R_argp = np.array([
        [np.cos(argp_rad), -np.sin(argp_rad), 0],
        [np.sin(argp_rad), np.cos(argp_rad), 0],
        [0, 0, 1]
    ])
    R_i = np.array([
        [1, 0, 0],
        [0, np.cos(i_rad), -np.sin(i_rad)],
        [0, np.sin(i_rad), np.cos(i_rad)]
    ])
    R_raan = np.array([
        [np.cos(raan_rad), -np.sin(raan_rad), 0],
        [np.sin(raan_rad), np.cos(raan_rad), 0],
        [0, 0, 1]
    ])

    # Transform to Earth-centered inertial frame
    R = R_raan @ R_i @ R_argp
    r_eci = R @ r_orbital
    v_eci = R @ v_orbital

    # Create Skyfield EarthSatellite object
    sat = EarthSatellite.from_vectors(wgs84, r_eci, v_eci, epoch)
    return sat

def compute_reference_quaternion(chief, target, t):
    # Compute positions and velocities of chief and target
    chief_gcrs = chief.at(t)
    target_gcrs = target.at(t)
    
    chief_pos = chief_gcrs.position.km
    chief_vel = chief_gcrs.velocity.km_per_s

    target_pos = target_gcrs.position.km

    # Compute the relative position vector (z-axis)
    z_axis = target_pos - chief_pos
    z_axis = z_axis / np.linalg.norm(z_axis)

    # Compute y-axis: cross product of z-axis and orbit normal
    y_axis = np.cross(z_axis, chief_vel)
    y_axis = y_axis / np.linalg.norm(y_axis)

    # Compute x-axis: cross product of y-axis and z-axis
    x_axis = np.cross(y_axis, z_axis)

    # Create the rotation matrix
    rotation_matrix = np.column_stack((x_axis, y_axis, z_axis))

    # Convert rotation matrix to quaternion
    quaternion = matrix_to_quaternion(rotation_matrix)

    # Ensure continuity of quaternion to avoid issues with numerical derivative
    quaternion = np.sign(quaternion[3])*quaternion 
    return quaternion, chief_pos, chief_vel, target_pos, x_axis, y_axis, z_axis

def ensure_quaternion_short_path(quaternions):
    quaternions_continuous = np.copy(quaternions)  # Copy the input array to preserve the original
    num_quaternions = quaternions.shape[0]
    
    for i in range(1, num_quaternions):
        # If the dot product of the current and previous quaternion is negative, flip the sign of the current quaternion
        if np.dot(quaternions_continuous[i], quaternions_continuous[i - 1]) < 0:
            quaternions_continuous[i] = -quaternions_continuous[i]
            # quaternions_continuous[i][3] = -quaternions_continuous[i][3]
    
    return quaternions_continuous

def compute_q_dot(quaternions, times):
    dt = (times[1] - times[0])*86400  # Assuming uniform time step / converting from days to seconds
    q_dot = np.zeros_like(quaternions)
    
    # Central finite difference
    q_dot[1:-1] = (quaternions[2:] - quaternions[:-2]) / (2 * dt)
    
    # Forward difference for the first point
    q_dot[0] = (quaternions[1] - quaternions[0]) / dt
    
    # Backward difference for the last point
    q_dot[-1] = (quaternions[-1] - quaternions[-2]) / dt
    
    return q_dot

def compute_angular_rates(q, q_dot):
    w = np.zeros((len(q), 3))
    for i in range(len(q)):
        q_matrix = np.array([
            [-q[i,1], -q[i,2], -q[i,3]],
            [q[i,0], -q[i,3], q[i,2]],
            [q[i,3], q[i,0], -q[i,1]],
            [-q[i,2], q[i,1], q[i,0]]
        ])
        w[i] = 2 * np.dot(q_matrix.T, q_dot[i])[:3]
    return w

def compute_angular_acceleration(angular_rates, times):
    # Ensure angular_rates is an Nx3 array and times is a 1D array of length N
    N = len(times)
    angular_acceleration = np.zeros_like(angular_rates)
    
    # Compute angular acceleration using central differences for inner points
    for i in range(1, N - 1):
        delta_t = times[i + 1] - times[i - 1]
        angular_acceleration[i] = (angular_rates[i + 1] - angular_rates[i - 1]) / delta_t
    
    # For the boundaries, use forward and backward differences
    angular_acceleration[0] = (angular_rates[1] - angular_rates[0]) / (times[1] - times[0])
    angular_acceleration[-1] = (angular_rates[-1] - angular_rates[-2]) / (times[-1] - times[-2])

    return angular_acceleration

def plot_satellites(chief_pos, chief_quaternions, target_pos, time_steps, args):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Set limits for the 3D plot
    ax.set_xlim([-30000, 30000])
    ax.set_ylim([-30000, 30000])
    ax.set_zlim([-30000, 30000])


    # Plot labels
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    # Initialize the chief satellite position and orientation
    chief_body, = ax.plot([], [], [], 'bo', markersize=8, label='Chief Satellite')

    # Initialize the target satellite position
    target_body, = ax.plot([], [], [], 'ro', markersize=8, label='Target Satellite')

    # Plot the orbits of both satellites
    ax.plot(chief_pos[:, 0], chief_pos[:, 1], chief_pos[:, 2], label="Chief Orbit", color='blue')
    ax.plot(target_pos[:, 0], target_pos[:, 1], target_pos[:, 2], label="Target Orbit", color='red')

    quivers = []  # To store quivers for later removal

    def update(frame):
        # Remove previous quivers if any
        for quiver in quivers:
            quiver.remove()
        quivers.clear()

        # Update the position of the chief satellite
        chief_body.set_data(chief_pos[frame][0], chief_pos[frame][1])
        chief_body.set_3d_properties(chief_pos[frame][2])

        # Update the position of the target satellite
        target_body.set_data(target_pos[frame][0], target_pos[frame][1])
        target_body.set_3d_properties(target_pos[frame][2])

        # Compute the orientation of the chief satellite using the quaternion
        rotation_matrix = quaternion_to_matrix(chief_quaternions[frame])

        # Draw the arrows representing the chief's orientation (using the quaternion)
        x_axis = rotation_matrix[:, 0]  # x-axis in the body frame
        y_axis = rotation_matrix[:, 1]  # y-axis in the body frame
        z_axis = rotation_matrix[:, 2]  # z-axis in the body frame

        quivers.append(ax.quiver(
            chief_pos[frame][0], chief_pos[frame][1], chief_pos[frame][2], 
            x_axis[0], x_axis[1], x_axis[2], color='red', length=50000
        ))
        quivers.append(ax.quiver(
            chief_pos[frame][0], chief_pos[frame][1], chief_pos[frame][2], 
            y_axis[0], y_axis[1], y_axis[2], color='green', length=50000
        ))
        quivers.append(ax.quiver(
            chief_pos[frame][0], chief_pos[frame][1], chief_pos[frame][2], 
            z_axis[0], z_axis[1], z_axis[2], color='blue', length=50000
        ))

        return chief_body, target_body

    ani = FuncAnimation(fig, update, frames=len(time_steps), interval=0.1, blit=False)
    
    plt.legend()
    plt.show()

    if args.save_ani:
         # Set up the writer
        writer = FFMpegWriter(fps=30)

        # Save the animation as an mp4 file
        ani.save('animation.mp4', writer=writer)
   
    

def main():

    # Create argument parser
    parser = argparse.ArgumentParser(description="Available options: --show_ani, --save_ani")
    
    # Add arguments
    parser.add_argument("--show_ani", action="store_true", help="show animation")
    parser.add_argument("--save_ani", action="store_true", help="save animation - you need to enable --show_ani")
    parser.add_argument("--full_telem", action="store_true", help="save csv with full telemetry")

    # Parse the arguments
    args = parser.parse_args()

    # Configuration
    use_tle = True  # Set to False to use orbital elements

    # Set time range and frequency
    # Load time scale
    ts = load.timescale()
    t_start = ts.tt_jd(2459000.5)  # Example start time
    t_end = ts.tt_jd(2459000.6)    # Example end time (1 day later)
    frequency_s = 0.5  # Hz - samples per second

    if use_tle:
        # Example TLEs (replace with actual TLEs)
        chief_tle_line1 = "1 43797U 18099Q   23262.12345678 -.00000012  00000-0 -34562-5 0  9993"
        chief_tle_line2 = "2 43797  97.7465 239.6012 0001247  76.9465 283.1880 15.21991542261234"
        target_tle_line1 = "1 20044U 89067A   23261.50000000  .00000000  00000-0  00000-0 0  9999"
        target_tle_line2 = "2 20044   0.0142 335.7784 0001641 216.5810 276.9542  1.00273272123456"

        chief = tle_to_skyfield(chief_tle_line1, chief_tle_line2)
        target = tle_to_skyfield(target_tle_line1, target_tle_line2)
    else:
        # Example orbital elements (replace with actual values)
        epoch = t_start # ts.tt_jd(2459000.5)
        chief = orbital_elements_to_skyfield(6778, 0.001, 51.6, 200, 0, 0, epoch)
        target = orbital_elements_to_skyfield(6778, 0.001, 51.6, 201, 0, 0, epoch)

    frequency = frequency_s*86400
    num_points = int((t_end.tt - t_start.tt) * frequency) + 1
    times = ts.tt_jd(np.linspace(t_start.tt, t_end.tt, num_points))

    # Compute reference quaternions
    results = [compute_reference_quaternion(chief, target, t) for t in times]
    quaternions = np.array([result[0] for result in results])
    
    # Ensure shortest path - avoid jumps in reference
    quaternions = ensure_quaternion_short_path(quaternions)

    chief_pos = np.array([result[1] for result in results])
    chief_vel = np.array([result[2] for result in results])
    target_pos = np.array([result[3] for result in results])
    x_vec = np.array([result[4] for result in results])
    y_vec = np.array([result[5] for result in results])
    z_vec = np.array([result[6] for result in results])


    # Compute q_dot using central finite differences
    q_dot = compute_q_dot(quaternions, times.tt)

    # Compute angular rates
    angular_rates = compute_angular_rates(quaternions, q_dot)
    angular_accelerations = compute_angular_acceleration(angular_rates, times.tt)

    # # Print results (you might want to save these to a file instead)
    # print("Time (JD)\tQuaternion\t\t\t\tAngular Rates (rad/s)")
    # for i in range(num_points):
    #     print(f"{times.tt[i]:.6f}\t{quaternions[i]}\t{angular_rates[i]}")

     # Save results to CSV file
    with open('satellite_tracking_data.csv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        
        if args.full_telem:
            # full telemetry
            csvwriter.writerow(['jdate', 'q_ref0', 'q_ref1', 'q_ref2', 'q_ref3', 'w_ref0', 'w_ref1', 'w_ref2', 'a_ref0', 'a_ref1', 'a_ref2', 'chief_pos0', 'chief_pos1', 'chief_pos2', 'chief_vel0', 'chief_vel1', 'chief_vel2', 'target_pos0', 'target_pos1', 'target_pos2', 'x_vec0', 'x_vec1', 'x_vec2', 'y_vec0', 'y_vec1', 'y_vec2', 'z_vec0', 'z_vec1', 'z_vec2',])

            # Write data
            for i in range(num_points):
                csvwriter.writerow([
                    times.tt[i],
                    quaternions[i][0], quaternions[i][1], quaternions[i][2], quaternions[i][3],
                    angular_rates[i][0], angular_rates[i][1], angular_rates[i][2],
                    angular_accelerations[i][0], angular_accelerations[i][1], angular_accelerations[i][2],

                    chief_pos[i][0], chief_pos[i][1], chief_pos[i][2],
                    chief_vel[i][0], chief_vel[i][1], chief_vel[i][2],
                    target_pos[i][0], target_pos[i][1], target_pos[i][2],
                    x_vec[i][0], x_vec[i][1], x_vec[i][2],
                    y_vec[i][0], y_vec[i][1], y_vec[i][2],
                    z_vec[i][0], z_vec[i][1], z_vec[i][2],
                ])

        else:
            # only controller data
            csvwriter.writerow(['jdate', 'q_ref0', 'q_ref1', 'q_ref2', 'q_ref3', 'w_ref0', 'w_ref1', 'w_ref2', 'a_ref0', 'a_ref1', 'a_ref2',])

            # Write data
            for i in range(num_points):
                csvwriter.writerow([
                    times.tt[i],
                    quaternions[i][0], quaternions[i][1], quaternions[i][2], quaternions[i][3],
                    angular_rates[i][0], angular_rates[i][1], angular_rates[i][2],
                    angular_accelerations[i][0], angular_accelerations[i][1], angular_accelerations[i][2],
                ])

    print(f"Data has been saved to 'satellite_tracking_data.csv'")

    if args.show_ani:
        # Call the function with the actual chief_pos, chief_quaternions, and target_pos data
        plot_satellites(chief_pos, quaternions, target_pos, times, args)

if __name__ == "__main__":
    main()
