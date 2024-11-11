using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Networking;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.IO;
using System.IO.Compression;
using System.Linq;

public class GPS : MonoBehaviour
{
    public GameObject Multicopter;
    public GameObject estMulticopter;
    public List<((double, double, double), (double, double, double))> SatVisibleECEF_ENU = new List<((double, double, double), (double, double, double))>();  // in ECEF, not blocked by Earth
    static int layermask = 1 << 6; // assuming obstacle (incl. its colliders) is on layer 6
    (double, double, double) estMultiPosECEF;
    (double, double, double) estMultiPosENU;
    private static string validGIMpath;

    public List<(double, double, double, double)> plot_list = new List<(double, double, double, double)>();

    public class TECMap
    {
        double[] Config = new double[5];
        double[] TECmap = new double[73];
    }

    #region EventHandler
    private void OnEnable()
    {
        SatelliteSetup.SetupDone += RunGPS;
    }

    private void OnDisable()
    {
        SatelliteSetup.SetupDone -= RunGPS;
    }
    #endregion

    private void RunGPS()
    {
        // Code to execute after setup is complete
        // Acts as main function
        Debug.Log("[GPS] Running GPS script logic.");

        StartCoroutine(GetGIM());                                   // Essentiel
        //double TEC = GetTEC((45, 30, 450));                       // Enter Lat/Long/Alt from IPP

        // Initial estimated Multicopter position in ECEF, Multicopter at ENU Origin / Location
        estMultiPosECEF = SatelliteSetup.LLA2ECEF(SatelliteSetup.LocLLA);
        //estMultiPosECEF = SatelliteSetup.ENU2ECEF((10, 20, 15));

        // Run update loop
        StartCoroutine(Update_Loop());
    }

    IEnumerator Update_Loop()
    {
        while (true)
        {
            var trueMultiPosENU = Multicopter.transform.position;   // Vector3
            var trueMultiPosECEF = SatelliteSetup.ENU2ECEF(((double)trueMultiPosENU.x, (double)trueMultiPosENU.y, (double)trueMultiPosENU.z)); // (double,double,double)

            SatVisibleECEF_ENU.Clear();

            // Add visible satellite to list
            foreach (var sat in SatelliteSetup.SatAboveECEF_ENU)
            {
                if (Visible(trueMultiPosENU, sat.Item2)) // second tuple contains ENU
                {
                    SatVisibleECEF_ENU.Add(sat);
                }
            }

            // At least 4 satellite have to be visible for a calculation
            if (SatVisibleECEF_ENU.Count >= 4)
            {
                Debug.Log($"Anzahl der Sichtbaren Satelliten {SatVisibleECEF_ENU.Count}");
                // Output position and error
                var error = posError(SatVisibleECEF_ENU, estMultiPosECEF, trueMultiPosECEF); //xyz_U and xyz_0
                //Debug.Log($"User Position Error: {error}");

                // Compute the estimated user position, xyz_u = xyz_0 + delta_xyz
                estMultiPosECEF = (estMultiPosECEF.Item1 + error.Item1, estMultiPosECEF.Item2 + error.Item2, estMultiPosECEF.Item3 + error.Item3);

                // Print
                estMultiPosENU = SatelliteSetup.ECEF2ENU(estMultiPosECEF);
                //Debug.Log($"Sim Multi ENU: {SatelliteSetup.ECEF2ENU(trueMultiPosECEF)}");
                //Debug.Log($"Est Multi ENU: {estMultiPosENU}");

                // Visualize
                estMulticopter.transform.position = new Vector3((float)estMultiPosENU.Item1, (float)estMultiPosENU.Item2, (float)estMultiPosENU.Item3);

                // Calculate the distance between vectorA and vectorB
                float distance = Vector3.Distance(trueMultiPosENU, new Vector3((float)estMultiPosENU.Item1, (float)estMultiPosENU.Item2, (float)estMultiPosENU.Item3));
                plot_list.Add(((double)trueMultiPosENU.x, (double)trueMultiPosENU.z, (double)trueMultiPosENU.y, (double)distance));     // convert X-North, Y-Up, Z-East to X-North, Y-East, Z-Up

                // Print the distance to the console
                Debug.Log("Distance between true and est Positation: " + distance);
            }

            yield return new WaitForSeconds(0.25f);
        }
    }

    private bool Visible(Vector3 multi, (double, double, double) satENU)
    {
        Vector3 sat = new Vector3((float)satENU.Item1, (float)satENU.Item2, (float)satENU.Item3);

        // Check if satellite signal is not blocked by any obstacle
        if (Physics.Linecast(multi, sat, layermask))
        {
            return false;
        }

        return true;
    }

    private (double, double, double, double) posError(List<((double, double, double), (double, double, double))> SatVisible, (double, double, double) estMulti, (double, double, double) trueMulti)
    {
        // Need to make sure poserror is only running when at least 4 satellites are visible

        //Debug.Log($"Est. ECEF: {estMulti} | True ECEF: {trueMulti}");

        // Current aka. last estimated position of multicopter
        double x_u = trueMulti.Item1;
        double y_u = trueMulti.Item2;
        double z_u = trueMulti.Item3;

        // Current true position of multicopter
        double x_0 = estMulti.Item1;
        double y_0 = estMulti.Item2;
        double z_0 = estMulti.Item3;

        // Number of visible satellites & counter
        int N = SatVisible.Count;
        int i = 0;

        // Create holder matrices for delta_p (nx1), H (nx4)
        var delta_p = DenseMatrix.Create(N, 1, 0.0);
        var H = DenseMatrix.Create(N, 4, 0.0);

        foreach (var sat in SatVisible)
        {
            var sat_ECEF = sat.Item1;       // ECEF tuple
            var sat_ENU = sat.Item2;

            var x_sat = sat_ECEF.Item1;
            var y_sat = sat_ECEF.Item2;
            var z_sat = sat_ECEF.Item3;

            // Simulate measured pseudoranges
            double p_i = Math.Sqrt(Math.Pow(x_sat - x_u, 2) + Math.Pow(y_sat - y_u, 2) + Math.Pow(z_sat - z_u, 2)) + Error(trueMulti, sat);      // 18.8 error in values, adjust later
            // Compute f_xyz_0
            double f_xyz_0 = Math.Sqrt(Math.Pow(x_sat - x_0, 2) + Math.Pow(y_sat - y_0, 2) + Math.Pow(z_sat - z_0, 2));
            // Compute Pseudo Range Error
            double p_error_i = f_xyz_0 - p_i;

            // Compute a_x, a_y, a_z
            double a_xi = (x_sat - x_0) / f_xyz_0;
            double a_yi = (y_sat - y_0) / f_xyz_0;
            double a_zi = (z_sat - z_0) / f_xyz_0;

            // Add to matrices
            delta_p[i, 0] = p_error_i;

            H[i, 0] = a_xi;
            H[i, 1] = a_yi;
            H[i, 2] = a_zi;
            H[i, 3] = 1.0;

            // increase counter
            i += 1;
        }

        if (N == 4)
        {
            var delta_x = H.Inverse() * delta_p;
            return (delta_x[0, 0], delta_x[1, 0], delta_x[2, 0], delta_x[3, 0]); // Only outputting x,y,z | -ctau_u (clock_offset error?) not needed
        }
        else
        {
            var delta_x = (H.Transpose() * H).Inverse() * H.Transpose() * delta_p;
            return (delta_x[0, 0], delta_x[1, 0], delta_x[2, 0], delta_x[3, 0]);
        }
    }

    public double Error((double, double, double) multi_ECEF, ((double, double, double), (double, double, double)) sat)
    {
        // Compute errors dynamically ECEF

        // Clock
        // Ephemeris

        // Ionospheric
        double I = IonoDelay(multi_ECEF, sat);


        // Tropospheric
        double T = 0;//TropoDelay(multi_ECEF, sat);

        // Multipath
        // Relativity

        // Noise
        float noise = Functions.Random(0, 1); //in m

        double E = (double)noise;

        return I + T;
    }

    IEnumerator GetGIM()
    {
        // Checks if the right GIM exists, if not try do download, if it doesn't exist, download the previous ones until one is valid

        DateTime utc = SatelliteSetup.UTC;  //DateTime.UtcNow; 
        int year = utc.Year;
        int dayofyear = utc.DayOfYear;
        string year_string = year.ToString();
        string dayofyear_string = dayofyear.ToString("D3");

        string username = "username";
        string password = "password";

        // Create the folder GIM under assets if it doesn't exist
        string directoryPath = Path.Combine(Application.dataPath, "GIM");
        if (!Directory.Exists(directoryPath))
        {
            Directory.CreateDirectory(directoryPath);
        }

        // Estimated Filepath of the correct GIM
        string filePath = Path.Combine(directoryPath, $"IGS0OPSRAP_{year_string}{dayofyear_string}0000_01D_02H_GIM.INX");

        // Check if the file exists
        if (File.Exists(filePath))
        {
            Debug.Log($"Found existing IGS0OPSRAP_{year_string}{dayofyear_string}0000_01D_02H_GIM.INX");
            validGIMpath = filePath;                        // confirmed Filepath of the valid GIM file
        }
        else
        {
            Debug.Log($"Requesting GIM for {utc}");

            // Run until a valid GIM is received
            while (true)
            {
                // Create a UnityWebRequest with Basic Authentication
                string resource = $"https://cddis.nasa.gov/archive/gnss/products/ionex/{year_string}/{dayofyear_string}/IGS0OPSRAP_{year_string}{dayofyear_string}0000_01D_02H_GIM.INX.gz";
                UnityWebRequest request = UnityWebRequest.Get(resource);
                string encodedAuth = System.Convert.ToBase64String(System.Text.Encoding.ASCII.GetBytes(username + ":" + password));
                request.SetRequestHeader("Authorization", "Basic " + encodedAuth);

                yield return request.SendWebRequest();

                if (request.result != UnityWebRequest.Result.Success)
                {
                    Debug.Log($"Request failed with error: {request.error} - {resource}");

                    // reduce dayofyear by 1, if it's 0 already, subtract on from the year and set day of year to 365 or 366 depending on the year
                    dayofyear--;
                    dayofyear_string = dayofyear.ToString("D3");

                    if (dayofyear < 1)
                    {
                        year--;
                        year_string = year.ToString();

                        if (DateTime.IsLeapYear(year))                  // Check if the previous year was a leap year
                        {
                            dayofyear = 366;
                            dayofyear_string = dayofyear.ToString("D3");
                        }
                        else
                        {
                            dayofyear = 365;
                            dayofyear_string = dayofyear.ToString("D3");
                        }
                    }
                }
                else
                {
                    // Set the path to save the file, potentially adjusted for previous GIMs, recall for updated path
                    filePath = Path.Combine(directoryPath, $"IGS0OPSRAP_{year_string}{dayofyear_string}0000_01D_02H_GIM.INX.gz");

                    // Save the file
                    File.WriteAllBytes(filePath, request.downloadHandler.data);

                    // Optional: If the file is compressed, decompress it
                    if (filePath.EndsWith(".gz"))
                    {
                        validGIMpath = filePath.Replace(".gz", "");                 // confirmed Filepath of the valid GIM file
                        DecompressGZip(filePath, validGIMpath);

                        // Print
                        DateTime validGIMEpoch = new DateTime(year, 1, 1);
                        validGIMEpoch = validGIMEpoch.AddDays(dayofyear - 1);      // Add the dayOfYear - 1 (because day 1 is already included)
                        Debug.Log($"Next valid GIM found for {validGIMEpoch} and saved as IGS0OPSRAP_{year_string}{dayofyear_string}0000_01D_02H_GIM.INX");

                        // Delete the compressed .gz
                        File.Delete(filePath);
                    }
                    break;
                }
            }
        }
    }

    private void DecompressGZip(string sourceFile, string destinationFile)
    {
        using (FileStream originalFileStream = new FileStream(sourceFile, FileMode.Open, FileAccess.Read))
        using (FileStream decompressedFileStream = new FileStream(destinationFile, FileMode.Create, FileAccess.Write))
        using (GZipStream decompressionStream = new GZipStream(originalFileStream, CompressionMode.Decompress))
        {
            decompressionStream.CopyTo(decompressedFileStream);
        }
    }

    public static double GetTEC((double, double, double) IPP)
    {
        // Ionicsphere Piercing Point in Lat/Long/Alt (here Alt is 450km)
        DateTime utc = SatelliteSetup.UTC;

        int hour = -1;                                            // Holder
        int mapnumber = -1;                                       // Holder
        int maphour;                                              // Target Map
        bool storeTEC = false;
        string TECstring = "";                                    // Holder

        // compute resulting map time (hour)
        if (utc.Hour % 2 == 0)
        {
            maphour = utc.Hour;
        }
        else
        {
            maphour = utc.Hour + 1;                 // if maphour is bigger than 23 aka next day, use different trigger 13th tec map in file
        }

        // Read file
        string[] lines = File.ReadAllLines(validGIMpath);

        foreach (string line in lines)
        {
            if (line.Contains("START OF TEC MAP"))
            {
                string[] lineParts = line.Trim().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                mapnumber = int.Parse(lineParts[0]);
            }
            else if (line.Contains("EPOCH OF CURRENT MAP"))
            {
                // Split the line into an array of strings, separating by whitespace
                string[] dateParts = line.Trim().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                // Get the hour of current map
                hour = int.Parse(dateParts[3]);
            }
            else if (line.Contains("LAT/LON1/LON2/DLON/H"))
            {
                // Parse the latitude and longitude range, e.g., 87.5-180.0 180.0 5.0 450.0
                string newline = line.Replace("-", " -");
                string[] info = newline.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                double lat = double.Parse(info[0]);

                // Debug.Log($"IPP Latitude: {IPP.Item1}, Map Latitude: {lat}, Difference: {Math.Abs(IPP.Item1 - lat)}");
                // Check if hour and latitude are correct
                if (maphour == hour && Math.Abs(IPP.Item1 - lat) < 2.5)                          // Trigger 1
                {
                    storeTEC = true;
                }

                if (maphour > 23 && mapnumber == 13 && Math.Abs(IPP.Item1 - lat) < 2.5)
                {
                    storeTEC = true;                                                             // Trigger 2, next day fits better
                }

                // Reading stops if next "LAT/LON1/LON2/DLON/H" is being read which are invalid, if TEC_string is not empty >> previous map is stored
                if (TECstring != "")
                {
                    storeTEC = false;
                    break;                                                                       // Break 1
                }
            }
            else if (line.Contains("END OF TEC MAP") && TECstring != "")                         // Different trigger for last map in epoch
            {
                storeTEC = false;
                break;                                                                           // Break 2
            }
            else if (storeTEC == true)
            {
                // Store TEC, add all lines to a long string
                TECstring += " " + line;
            }
        }
        // extract TEC from string
        string[] TECs = TECstring.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

        // Mapping Longitude from -180 to 180 into 0 to 72 (Index)
        int TECIndex = (int)Math.Round(Functions.MapValue(IPP.Item2, -180, 180, 0, 72));
        double TEC = double.Parse(TECs[TECIndex]) * 0.1;        // TEC values in 0.1 TECU in used GIM
        return TEC;
    }

    public double IonoDelay((double, double, double) multi_ECEF, ((double, double, double), (double, double, double)) sat)
    {
        var sat_ECEF = sat.Item1;
        var sat_ENU = sat.Item2;
        var multi_ENU = SatelliteSetup.ECEF2ENU(multi_ECEF);

        double r_earth = 6371e3;               //      in m
        double H_max = 450e3;                  //      Height of max Electron density in m (source used GIM)

        // Calculate Satellite Elevation
        // Vector from receiver to satellite
        double deltaX = sat_ENU.Item1 - multi_ENU.Item1;
        double deltaY = sat_ENU.Item2 - multi_ENU.Item2;
        double deltaZ = sat_ENU.Item3 - multi_ENU.Item3;

        // Elevation angle in radians
        double E = Math.Atan2(deltaY, Math.Sqrt(deltaX * deltaX + deltaZ * deltaZ));            // Vertical divided by horitzontal

        //Azimuth
        double A = Math.Atan2(deltaZ, deltaX) * 180.0 / Math.PI;     // in degree

        // Ensure azimuth is within the range 0 to 360 degrees
        if (A < 0)
        {
            A += 360;
        }

        // Convert back to radians
        A = A * Math.PI / 180.0;

        // IPP Calculation

        var multi_LLA = SatelliteSetup.ECEF2LLA(sat_ECEF);   // Geographic Coord of Multicopter
        double lat = multi_LLA.Item1 * Math.PI / 180.0;                 // radians
        double lon = multi_LLA.Item2 * Math.PI / 180.0;

        double psi = Math.PI / 2 - E - Math.Asin(r_earth / (r_earth + H_max) * Math.Cos(E));
        double latIPP = Math.Asin(Math.Sin(lat) * Math.Cos(psi) + Math.Cos(lat) * Math.Sin(psi) * Math.Cos(A));

        double lonIPP = lon + Math.Asin(Math.Sin(psi) * Math.Sin(A) / Math.Cos(latIPP));

        // Convert IPP latitude and longitude from radians to degrees if needed
        latIPP = latIPP * 180.0 / Math.PI;
        lonIPP = lonIPP * 180.0 / Math.PI;

        // Ensure longitude is within the -180 to 180 range
        if (lonIPP > 180.0)
        {
            lonIPP -= 360.0;
        }
        else if (lonIPP < -180.0)
        {
            lonIPP += 360.0;
        }

        // Debug.Log($"Latitude of IPP {latIPP}");
        // Debug.Log($"Longitude of IPP {lonIPP}");

        double TEC = GetTEC((latIPP, lonIPP, H_max));        // in TECU 10e16 Electrons / m^2
        //Debug.Log($"TEC: {TEC}");

        // Mapping Function
        double z = Math.PI / 2 - E;                // in Radian
        double zdash = Math.Asin(r_earth / (r_earth + H_max) * Math.Sin(z));
        double M = 1 / Math.Cos(zdash);

        // Computing STEC
        double STEC = TEC / M;
        // Debug.Log($"STEC: {STEC}");

        // Ionospheric Delay
        double k = 40.3;                 // constant 40.3 m^3/s^2
        double fc = 1575.42e6;            // Carrier Frequency in Hz, this GPS module uses L1 1575.42 MHz
        double I = (k * STEC * 1e16) / (fc * fc);

        //Debug.Log($"Ionospheric Delay: {I}");

        return I;
    }

    public double TropoDelay((double, double, double) multi_ECEF, ((double, double, double), (double, double, double)) sat)
    {
        double P0 = 101325; // Sea level standard atmospheric pressure in Pa
        double alpha = 6.5e-3; // Temperature lapse rate in K/m = °C/m
        double T0 = 288.15; // MSL standard temperature in K
        double g = 9.80665; // Earth's gravitational acceleration in m/s^2
        double M = 0.02896; // Molar mass of Earth's air in kg/mol
        double R = 8.314f; // Universal gas constant in J/(mol·K)

        var sat_ENU = sat.Item2;
        var multi_ENU = SatelliteSetup.ECEF2ENU(multi_ECEF);

        double h = multi_ENU.Item2;

        double p = P0 * Math.Pow((1 - (alpha * h) / T0), (g * M) / (R * alpha));               // Pressure in Pa

        double deltaX = sat_ENU.Item1 - multi_ENU.Item1;
        double deltaY = sat_ENU.Item2 - multi_ENU.Item2;
        double deltaZ = sat_ENU.Item3 - multi_ENU.Item3;

        double E = Math.Atan2(deltaY, Math.Sqrt(deltaX * deltaX + deltaZ * deltaZ));            // Elevation angle in radians

        double T = (T0 - 273.15) - alpha * h;                                                   // Temperature in °C
        double T_k = T + 273.15;                                                                // Temperature in K

        double e = Math.Pow(0.61094, (17.625 * T) / (243.04 + T));                              // Vapor Pressure in kPa

        // Saastamoinen Model, Enter T in K, p and e in mbar
        double E0 = E + 16 / T_k * (p / 100 + (4810 / T_k) * e * 10) * (1 / Math.Tan(E));
        double beta = 1.16 - 0.15e-3 * h + 0.716e-8 * h * h;
        double Td = (0.002277 / Math.Sin(E0)) * (p / 100 + ((1255 / T_k) + 0.05) * e * 10 - (beta / Math.Pow(Math.Tan(E0), 2)));


        Debug.Log($"h = {h}m");
        Debug.Log($"p = {p}Pa");
        Debug.Log($"E = {E * 180 / Math.PI}°");
        Debug.Log($"T = {T}°C");
        Debug.Log($"e = {e}kPa");
        Debug.Log($"Tropospheric Delay: {Td}");
        return Td;
    }

    void OnApplicationQuit()    // Export when scenario is stopped
    {
        Debug.Log("Exporting");
        using (StreamWriter writer = new StreamWriter("Assets/gps_points.csv"))
        {
            foreach (var pair in plot_list)
            {
                string combinedline = $"{pair.Item1},{pair.Item2},{pair.Item3},{pair.Item4}";
                writer.WriteLine(combinedline);
            }
        }
    }
}