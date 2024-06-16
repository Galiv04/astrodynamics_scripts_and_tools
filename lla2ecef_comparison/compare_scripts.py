import subprocess
import json
import math

def run_c_script(c_file, lat, lon, alt):
    # Ensure GCC is in the PATH or provide the full path to gcc
    gcc_path = 'gcc'  # Update this if necessary to the full path to gcc, e.g., 'C:\\MinGW\\bin\\gcc.exe'

    # Compile the C program
    compile_process = subprocess.run([gcc_path, c_file, '-o', 'c_program'], capture_output=True, text=True)
    if compile_process.returncode != 0:
        print(f"Compilation failed:\n{compile_process.stderr}")
        return None

    # Run the compiled C program with arguments
    run_process = subprocess.run(['./c_program', str(lat), str(lon), str(alt)], capture_output=True, text=True)
    if run_process.returncode != 0:
        print(f"Running C program failed:\n{run_process.stderr}")
        return None

    return run_process.stdout.strip()

def run_python_script(py_file, lat, lon, alt):
    run_process = subprocess.run(['python', py_file, str(lat), str(lon), str(alt)], capture_output=True, text=True)
    if run_process.returncode != 0:
        print(f"Running Python script failed:\n{run_process.stderr}")
        return None

    return run_process.stdout.strip()

def parse_output(output):
    try:
        return json.loads(output)
    except json.JSONDecodeError:
        print("Failed to parse output as JSON.")
        return None
        
def calculate_difference(c_output_parsed, py_output_parsed):
    diff_x = c_output_parsed["x"] - py_output_parsed["x"]
    diff_y = c_output_parsed["y"] - py_output_parsed["y"]
    diff_z = c_output_parsed["z"] - py_output_parsed["z"]
    return diff_x, diff_y, diff_z

def calculate_ground_distance(diff_x, diff_y, diff_z):
    return math.sqrt(diff_x**2 + diff_y**2 + diff_z**2)  # Convert meters to kilometers

def compare_outputs(c_output, py_output):
    c_output_parsed = parse_output(c_output)
    py_output_parsed = parse_output(py_output)
    
    if c_output_parsed is None or py_output_parsed is None:
        print("Comparison failed due to parsing errors.")
        return

    if c_output_parsed == py_output_parsed:
        print("The outputs are identical.")
    else:
        print("The outputs differ.")
        diff_x, diff_y, diff_z = calculate_difference(c_output_parsed, py_output_parsed)
        ground_distance_km = calculate_ground_distance(diff_x, diff_y, diff_z)

        print(f"C output:\n{c_output_parsed}")
        print(f"Python output:\n{py_output_parsed}")

        print(f"Differences in coordinates:")
        print(f"  X difference: {diff_x:.6f} meters")
        print(f"  Y difference: {diff_y:.6f} meters")
        print(f"  Z difference: {diff_z:.6f} meters")
        print(f"Ground distance difference: {ground_distance_km:.6f} meters")

def main():
    # Define the latitude, longitude, and altitude
    lat = 80
    lon = -160
    alt = 500

    c_output = run_c_script('lla_to_ecef.c', lat, lon, alt)
    if c_output is None:
        return

    py_output = run_python_script('lla_to_ecef.py', lat, lon, alt)
    if py_output is None:
        return

    compare_outputs(c_output, py_output)

if __name__ == "__main__":
    main()
