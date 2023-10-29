from datetime import *

def julian_day_number_to_gregorian(jdn):
    """Convert the Julian Day Number to the proleptic Gregorian Year, Month, Day."""
    L = jdn + 68569
    N = int(4 * L / 146_097)
    L = L - int((146097 * N + 3) / 4)
    I = int(4000 * (L + 1) / 1_461_001)
    L = L - int(1461 * I / 4) + 31
    J = int(80 * L / 2447)
    day = L - int(2447 * J / 80)
    L = int(J / 11)
    month = J + 2 - 12 * L
    year = 100 * (N - 49) + I + L
    return year, month, day

def julian_date_to_gregorian(jd):
    """Convert a decimal Julian Date to the equivalent proleptic Gregorian date and time."""
    jdn = int(jd)
    if jdn < 1_721_426:
        raise ValueError("Julian Day Numbers less than 1,721,426 are not supported, "
                         "because Python's date class cannot represent years before "
                         "AD 1.")
    year, month, day = julian_day_number_to_gregorian(jdn)
    offset = timedelta(days=(jd % 1), hours=+12)
    dt = datetime(year=year, month=month, day=day, tzinfo=timezone.utc)
    return dt + offset

gregorian = julian_date_to_gregorian(2460188.0842702985)
print(gregorian)