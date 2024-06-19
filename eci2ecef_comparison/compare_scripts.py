import subprocess
import os
import numpy as np
from scipy.spatial.transform import Rotation as R

def compile_c_file(c_file, output_executable):
    try:
        subprocess.check_call(['gcc', c_file, '-o', output_executable])
        print(f"{c_file} compiled successfully to {output_executable}.")
    except subprocess.CalledProcessError as e:
        print(f"Error compiling {c_file}: {e}")
        return False
    return True

def run_executable(executable, JulianDate):
    try:
        result = subprocess.run([f'./{executable}', str(JulianDate)], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error running {executable}: {result.stderr}")
            return None
        return result.stdout.strip()
    except Exception as e:
        print(f"Error running {executable}: {e}")
        return None

def parse_quaternion(output):
    try:
        quaternion = [float(x) for x in output.split()]
        if len(quaternion) != 4:
            raise ValueError("Output is not a valid quaternion.")
        return quaternion
    except ValueError as e:
        print(f"Error parsing quaternion: {e}")
        return None
    
def compute_quaternion_error(q1, q2):
    q1 = R.from_quat(q1)
    q2 = R.from_quat(q2)
    q_error = q1.inv() * q2
    return q_error

def quaternion_to_euler(q):
    euler = q.as_euler('xyz', degrees=True)  # 'xyz' is a common sequence, adjust if needed
    return euler

def main():
    c_file1 = 'eci_to_ecef.c'
    c_file2 = 'eci_to_ecef_ref.c'
    exe1 = 'eci_to_ecef'
    exe2 = 'eci_to_ecef_ref'
    JulianDate = 2460484.5  # example input, replace as needed

    if not (compile_c_file(c_file1, exe1) and compile_c_file(c_file2, exe2)):
        return

    output1 = run_executable(exe1, JulianDate )
    output2 = run_executable(exe2, JulianDate )

    if output1 is None or output2 is None:
        return

    quaternion1 = parse_quaternion(output1)
    quaternion2 = parse_quaternion(output2)

    if quaternion1 is None or quaternion2 is None:
        return

    print("Quaternion from a3200 code:", quaternion1)
    print("Quaternion from ref code:", quaternion2)
    
    quaternion_error = compute_quaternion_error(quaternion1, quaternion2)
    euler_angles = quaternion_to_euler(quaternion_error)

    print("Quaternion error (as quaternion):", quaternion_error.as_quat())
    print("Euler angles of the quaternion error:")
    print("Roll (x):", euler_angles[0])
    print("Pitch (y):", euler_angles[1])
    print("Yaw (z):", euler_angles[2])

if __name__ == "__main__":
    main()
