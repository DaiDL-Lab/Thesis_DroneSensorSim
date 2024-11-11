using UnityEngine;
using UnityEngine.Networking;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using Newtonsoft.Json;

// Need to include MathNet.Numerics.dll + Newtonsoft.Json exactly under Assets/Plugins

public class SatelliteSetup : MonoBehaviour
{
    #region classes

    [System.Serializable]
    public class Info
    {
        public int satid;
        public string satname;
        public int transactionscount;
    }

    [System.Serializable]
    public class SatelliteData
    {
        public Info info;
        public string tle;
    }

    [System.Serializable]
    public class SatelliteTLE
    {
        public string Line1;
        public string Line2;
    }

    public class TLE_Line2
    {
        public int norad_id;
        public double inc;   // inclination in °
        public double asc;   // right ascention of ascending node in °
        public double ecc;   // eccentricity
        public double peri;  // argument of perigee in °
        public double ano;   // mean anomaly in °
        public double mot;   // mean motion in revolutions per day
    }

    #endregion

    private readonly string apiKey = "XXXXXX-XXXXXX-XXXXXX-XXXX"; // personal API-Key from n2yo.com
    // User: duongledai2000 | PW: #GIMrequest2024
    private string apiURLTLE; // This will be initialized in the method
    public string TLEDataSet; // CSV filename (without extension)
    private string FilePath;
    public static DateTime UTC;
    public string LocationLLA; // Latitude, Longitude and Altitude of the location of Input (ILR: 52.514726,13.323701,42)
    public static (double, double, double) LocLLA = (0, 0, 0);    // Holder
    private List<SatelliteTLE> satelliteTLEs = new List<SatelliteTLE>();
    private List<TLE_Line2> TLE_L2 = new List<TLE_Line2>();
    private List<(string Line1, string Line2)> TLE_Data = new List<(string Line1, string Line2)>();
    public static List<((double, double, double), (double, double, double))> SatAboveECEF_ENU = new List<((double, double, double), (double, double, double))>();  // in ECEF, not blocked by Earth
    private List<int> norad_ID = new List<int> // Include all active GPS norad IDs
    {
        48859, 46826, 45854, 44506, 43873, 41328, 41019, 40730, 40534, 40294, 40105, 39741, 39533, 39166, 38833,
        37753, 36585, 35752, 32711, 32384, 32260, 29601, 29486, 28874, 28474, 28190, 27704, 27663, 26605, 26360, 24876
    };
    public delegate void Eventhandler();                  // Define the delegate type, Event handler
    public static event Eventhandler SetupDone;           // Define the event based on the delegate

    void Start()
    {
        // Convert LLA input into doubles
        string[] loc_LLA = LocationLLA.Split(',');
        LocLLA = (double.Parse(loc_LLA[0]), double.Parse(loc_LLA[1]), double.Parse(loc_LLA[2]));

        //// Run function in coroutine for Convenient Timing Control
        StartCoroutine(HandleTLEData());
    }

    private IEnumerator HandleTLEData()
    {
        // (1) If no TLE Data is entered, retrieve new TLE Data >> export+import; (2) If an existing TLE Data is entered, import
        if (string.IsNullOrWhiteSpace(TLEDataSet))
        {
            // Start the coroutine to fetch satellite data
            Debug.Log("[Setup] TLE data is requested from n2yo.com.");
            yield return StartCoroutine(GetTLE()); // function needs to run completely

            // Set current UTC as Name for the Data Set, ISO 8601 would be "yyyy-MM-ddTHH:mm:ssZ", but can't have ":" in filenames
            TLEDataSet = DateTime.UtcNow.ToString("yyyy-MM-ddTHH-mm-ssZ");

            // Determine the file path
            FilePath = Path.Combine(Application.dataPath, TLEDataSet + ".csv");

            ExportTLE(FilePath);
            TLE_Data = ImportTLE(FilePath);
        }
        else
        {
            // Import CSV with TLE lines
            Debug.Log($"[Setup] Import existing TLE file {TLEDataSet}.csv");
            TLE_Data = ImportTLE(Path.Combine(Application.dataPath, TLEDataSet + ".csv"));
        }

        // Extracting all values from TLE Line 2 for each TLE
        foreach (var TLE in TLE_Data)
        {
            ExtractTLELine2Value(TLE);
        }

        // Convert TLE to ECEF and ECEF to ENU and store satellites
        foreach (var set in TLE_L2)
        {
            var ECI = TLE2ECI(set);
            var ECEF = ECI2ECEF(ECI, TLEDataSet);
            if (Above(ECEF))
            {
                // Satellites coord in local frame
                var ENU = ECEF2ENU(ECEF);
                SatAboveECEF_ENU.Add((ECEF, ENU));
            }
        }

        // Notify subscribers (other codes) that setup is complete
        Debug.Log("[Setup] Setup successfully completed.");
        SetupDone?.Invoke();
    }

    private IEnumerator GetTLE()
    {
        foreach (int id in norad_ID)
        {
            apiURLTLE = $"https://api.n2yo.com/rest/v1/satellite/tle/{id}?apiKey={apiKey}";

            using (UnityWebRequest webRequest = UnityWebRequest.Get(apiURLTLE))
            {
                yield return webRequest.SendWebRequest();

                string response = webRequest.downloadHandler.text;
                Debug.Log("[Setup] Received JSON: " + response);  // Log the JSON for debugging

                // Deserialize using Json.NET
                SatelliteData data = JsonConvert.DeserializeObject<SatelliteData>(response);

                string[] lines = data.tle.Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);

                string TLE_line1 = lines[0].Trim();
                string TLE_line2 = lines[1].Trim();

                SatelliteTLE newTLE = new SatelliteTLE
                {
                    Line1 = TLE_line1,
                    Line2 = TLE_line2
                };

                satelliteTLEs.Add(newTLE);
            }
        }
        Debug.Log("[Setup] TLE data received.");
    }

    private void ExportTLE(string FilePath)
    {
        using (StreamWriter writer = new StreamWriter(FilePath))
        {
            foreach (var tle in satelliteTLEs)
            {
                string combinedTLE = $"{tle.Line1.Replace("\r", "").Replace("\n", " ")},{tle.Line2.Replace("\r", "").Replace("\n", " ")}";
                writer.WriteLine(combinedTLE);
            }
        }
        Debug.Log($"[Setup] TLE data ({FilePath}) exported successfully.");
    }

    private List<(string Line1, string Line2)> ImportTLE(string FilePath)
    {
        List<(string Line1, string Line2)> TLE_List = new List<(string Line1, string Line2)>();

        using (StreamReader reader = new StreamReader(FilePath))
        {
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                // Split line by comma
                string[] parts = line.Split(',');

                // Add the tuple to the list
                TLE_List.Add((parts[0].Trim(), parts[1].Trim()));
            }
        }

        // Print all imported TLE for testing
        foreach (var tle in TLE_List)
        {
            //Debug.Log($"{tle.Line1}\n{tle.Line2}");
        }
        Debug.Log("[Setup] TLE data imported successfully.");
        return TLE_List;
    }

    private void ExtractTLELine2Value((string, string) TLE)
    {

        // Seperate the values from line 2
        var TLE_Line2 = TLE.Item2.Split(new[] { " " }, StringSplitOptions.RemoveEmptyEntries);

        //Debug.Log(TLE.Item2);

        TLE_L2.Add(new TLE_Line2
        {
            norad_id = int.Parse(TLE_Line2[1]),
            inc = double.Parse(TLE_Line2[2]),
            asc = double.Parse(TLE_Line2[3]),
            ecc = double.Parse("0." + TLE_Line2[4]),
            peri = double.Parse(TLE_Line2[5]),
            ano = double.Parse(TLE_Line2[6]),
            mot = double.Parse(TLE_Line2[7])
        });
    }

    public static (double, double, double) TLE2ECI(TLE_Line2 Line2)
    {
        // Calculations in double due to plugin MathNet.Numerics plugins
        // Convert angles to radians
        double peri_rad = Line2.peri * Math.PI / 180.0;
        double inc_rad = Line2.inc * Math.PI / 180.0;
        double asc_rad = Line2.asc * Math.PI / 180.0;
        double ano_rad = Line2.ano * Math.PI / 180.0;       // converge with function

        // Rotation matrices
        Matrix<double> R1 = DenseMatrix.OfArray(new double[,]
        {
        {Math.Cos(peri_rad), -Math.Sin(peri_rad), 0},
        {Math.Sin(peri_rad), Math.Cos(peri_rad), 0},
        {0, 0, 1}
        });

        Matrix<double> R2 = DenseMatrix.OfArray(new double[,]
        {
        {1, 0, 0},
        {0, Math.Cos(inc_rad), -Math.Sin(inc_rad)},
        {0, Math.Sin(inc_rad), Math.Cos(inc_rad)}
        });

        Matrix<double> R3 = DenseMatrix.OfArray(new double[,]
        {
        {Math.Cos(asc_rad), -Math.Sin(asc_rad), 0},
        {Math.Sin(asc_rad), Math.Cos(asc_rad), 0},
        {0, 0, 1}
        });

        double mu = 3.986005e14;    // Earth's gravitational parameter m^3/s^2

        // Calculation of the semi-major axis
        double n = Line2.mot * 2 * Math.PI / 86400;  // Convert rev/d to rad/s
        double a = Math.Pow(mu / Math.Pow(n, 2), 1.0 / 3.0);

        Matrix<double> positionMatrix = DenseMatrix.OfArray(new double[,]
        {
        {Math.Cos(ano_rad) - Line2.ecc},
        {Math.Sqrt(1 - Line2.ecc * Line2.ecc) * Math.Sin(ano_rad)},
        {0}
        });

        // Compute ECI coordinates
        Matrix<double> ECI = R3 * R2 * R1 * a * positionMatrix;
        return (ECI[0, 0], ECI[1, 0], ECI[2, 0]);
    }

    public static double JulianDate(string UTC_string)
    {
        // Convert to ISO 8601
        string[] parts = UTC_string.Split('T');
        string UTC_ISO = parts[0] + "T" + parts[1].Replace('-', ':');

        UTC = DateTime.Parse(UTC_ISO, null, System.Globalization.DateTimeStyles.RoundtripKind);

        int year = UTC.Year;
        int month = UTC.Month;
        int day = UTC.Day;
        double hour = UTC.Hour;
        double minute = UTC.Minute;
        double second = UTC.Second;

        if (month <= 2)
        {
            year--;
            month += 12;
        }

        int A = year / 100;
        int B = 2 - A + A / 4;

        double jd = Math.Floor(365.25 * (year + 4716)) + Math.Floor(30.6001 * (month + 1)) + day + B - 1524.5 + (hour / 24.0) + (minute / 1440.0) + (second / 86400.0);

        return jd;
    }

    public static (double, double, double) ECI2ECEF((double, double, double) ECI, string UTC_string)
    {

        // Convert UTC to Julian Date
        double jd = JulianDate(UTC_string);

        // Calculate Julian Century
        double T = (jd - 2451545.0) / 36525.0;

        // Compute Greenwich Mean Sidereal Time (GMST) in degrees
        double GMST = (280.46061837 + 360.98564736629 * (jd - 2451545.0) +
                       T * T * (0.000387933 - T / 38710000.0)) % 360.0;

        // Convert GMST to radians
        double GMSTRad = GMST * Math.PI / 180.0;

        // Compute the rotation matrix
        Matrix<double> R = DenseMatrix.OfArray(new double[,]
        {
            {Math.Cos(GMSTRad), Math.Sin(GMSTRad), 0},
            {-Math.Sin(GMSTRad), Math.Cos(GMSTRad), 0},
            {0, 0, 1}
        });

        Matrix<double> ECI_pos = DenseMatrix.OfArray(new double[,]
        {
        {ECI.Item1},
        {ECI.Item2},
        {ECI.Item3}
        });

        // Transform ECI to ECEF
        Matrix<double> ECEF = R * ECI_pos;

        return (ECEF[0, 0], ECEF[1, 0], ECEF[2, 0]);
    }

    public static (double, double, double) ECEF2LLA((double, double, double) ECEF)
    {
        double a = 6378137.0; // Semi-major axis in meters
        double e = 8.1819190842622e-2; // Eccentricity
        double b = a * Math.Sqrt(1 - e * e);
        double ep = Math.Sqrt((a * a - b * b) / (b * b));
        double p = Math.Sqrt(ECEF.Item1 * ECEF.Item1 + ECEF.Item2 * ECEF.Item2);
        double th = Math.Atan2(a * ECEF.Item3, b * p);
        double lon = Math.Atan2(ECEF.Item2, ECEF.Item1);
        double lat = Math.Atan2(ECEF.Item3 + ep * ep * b * Math.Pow(Math.Sin(th), 3), p - e * e * a * Math.Pow(Math.Cos(th), 3));
        double N = a / Math.Sqrt(1 - e * e * Math.Pow(Math.Sin(lat), 2));
        double alt = p / Math.Cos(lat) - N;

        // Convert radians to degrees
        lon = lon * 180.0 / Math.PI;
        lat = lat * 180.0 / Math.PI;

        return (lat, lon, alt);
    }

    public static (double, double, double) LLA2ECEF((double, double, double) LLA)
    {
        // Constants
        double a = 6378137.0; // Semi-major axis in meters
        double e = 8.1819190842622e-2; // Eccentricity

        // Convert lat/lon from degrees to radians
        double lat = LLA.Item1 * Math.PI / 180.0;
        double lon = LLA.Item2 * Math.PI / 180.0;

        // Compute N (prime vertical radius of curvature)
        double N = a / Math.Sqrt(1 - e * e * Math.Pow(Math.Sin(lat), 2));

        // Compute ECEF coordinates
        double x = (N + LLA.Item3) * Math.Cos(lat) * Math.Cos(lon);
        double y = (N + LLA.Item3) * Math.Cos(lat) * Math.Sin(lon);
        double z = ((1 - e * e) * N + LLA.Item3) * Math.Sin(lat);

        return (x, y, z);
    }

    public bool Above((double, double, double) sat_ECEF)
    {
        // Observer Location is predetermined, Location of site instead of location of receiver is fine due to scale (minimal errors)
        // using earth as sphere is also fine

        var loc_ECEF = LLA2ECEF((LocLLA.Item1, LocLLA.Item2, LocLLA.Item3));

        Matrix<double> loc = DenseMatrix.OfArray(new double[,]{
            {loc_ECEF.Item1},
            {loc_ECEF.Item2},
            {loc_ECEF.Item3}
            });

        Matrix<double> sat = DenseMatrix.OfArray(new double[,]{
            {sat_ECEF.Item1},
            {sat_ECEF.Item2},
            {sat_ECEF.Item3}
            });

        // Compute distance satellite to Earth Center
        double r_earth = 6371e3; //in m
        double r_sat = Math.Sqrt(Math.Pow(sat_ECEF.Item1, 2) + Math.Pow(sat_ECEF.Item2, 2) + Math.Pow(sat_ECEF.Item3, 2));

        double angle_crit = Math.Acos(r_earth / r_sat) * 180.0 / Math.PI;
        double angle_loc = Math.Acos((loc.Transpose() * loc + sat.Transpose() * sat - (sat - loc).Transpose() * (sat - loc)).At(0, 0) / (2 * Math.Sqrt((loc.Transpose() * loc * sat.Transpose() * sat).At(0, 0)))) * 180.0 / Math.PI;

        if (angle_loc <= angle_crit)
        {
            // Satellite is visible aka. Signal is not blocked by earth
            return true;
        }
        else
        {
            // Satellite isn't visible
            return false;
        }
    }

    public static (double, double, double) ECEF2ENU((double, double, double) ECEF)
    {
        var LocECEF = LLA2ECEF(LocLLA);

        // ENU origin in ECEF
        Matrix<double> ENU_origin = DenseMatrix.OfArray(new double[,]
        {
            {LocECEF.Item1},
            {LocECEF.Item2},
            {LocECEF.Item3}
        });

        Matrix<double> ECEF_point = DenseMatrix.OfArray(new double[,]
        {
            {ECEF.Item1},
            {ECEF.Item2},
            {ECEF.Item3}
        });

        // Subtract ECEF coordinates of the reference point
        var local = ECEF_point - ENU_origin;

        // Define rotation matrix
        double phi = LocLLA.Item1 * Math.PI / 180.0;       // Latitude
        double lambda = LocLLA.Item2 * Math.PI / 180.0;    // Longitude

        // Compute the rotation matrix
        Matrix<double> rotation = DenseMatrix.OfArray(new double[,]
        {
            {-Math.Sin(lambda), Math.Cos(lambda), 0},
            {-Math.Sin(phi) * Math.Cos(lambda), -Math.Sin(phi) * Math.Sin(lambda), Math.Cos(phi)},
            {Math.Cos(phi) * Math.Cos(lambda), Math.Cos(phi) * Math.Sin(lambda), Math.Sin(phi)}
        });

        // Apply rotation and add translation Vector
        var ENU = rotation * local;

        // Extract ENU coordinates from the resulting matrix
        double ENU_east = ENU[0, 0];            // Original X
        double ENU_north = ENU[1, 0];           // Original Y
        double ENU_up = ENU[2, 0];              // Original Z

        // adjusted for Unity, since Y and Z are flipped
        return (ENU_east, ENU_up, ENU_north);
    }

    public static (double, double, double) ENU2ECEF((double, double, double) ENU_Unity)
    {
        // Calculate the ENU origin in ECEF
        var LocECEF = LLA2ECEF(LocLLA);

        // Define the inverse rotation matrix
        double phi = LocLLA.Item1 * Math.PI / 180.0;       // Latitude
        double lambda = LocLLA.Item2 * Math.PI / 180.0;    // Longitude

        // Compute the inverse rotation matrix
        Matrix<double> rotationInverse = DenseMatrix.OfArray(new double[,]
        {
        {-Math.Sin(lambda), -Math.Sin(phi) * Math.Cos(lambda), Math.Cos(phi) * Math.Cos(lambda)},
        {Math.Cos(lambda), -Math.Sin(phi) * Math.Sin(lambda), Math.Cos(phi) * Math.Sin(lambda)},
        {0, Math.Cos(phi), Math.Sin(phi)}
        });

        // Convert ENU coordinates to local ECEF coordinates
        Matrix<double> ENU_matrix = DenseMatrix.OfArray(new double[,]
        {
        {ENU_Unity.Item1},
        {ENU_Unity.Item3},    // Adjust for the Y-Axis and Z-Axis flip in Unity
        {ENU_Unity.Item2}
        });

        var localECEF = rotationInverse * ENU_matrix;

        // Add the ENU origin in ECEF to get the final ECEF coordinates
        Matrix<double> ENU_origin = DenseMatrix.OfArray(new double[,]
        {
        {LocECEF.Item1},
        {LocECEF.Item2},
        {LocECEF.Item3}
        });

        var ECEF = localECEF + ENU_origin;

        return (ECEF[0, 0], ECEF[1, 0], ECEF[2, 0]);
    }
}
