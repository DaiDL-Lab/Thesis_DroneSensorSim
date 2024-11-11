using UnityEngine;
using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Globalization;
using Unity.Mathematics;
//using MathNet.Numerics;

public class WMM : MonoBehaviour
{
    static float lat = 52.514726f;  //°
    static float lng = 13.323701f;  //°
    static float alt = 42f;      //m
    /*
        static float lat = 50f;  //°
        static float lng = 13f;  //°
        static float alt = 42f;      //m
    */
    public GameObject multicopter;

    // Setting WGS-84 parameters
    public static float a = 6378.137f; //semi-major axis of the ellipsoid in km
    public static float b = 6356.7523142f; // semi-minor axis of the ellipsoid in km
    public static float f = 1f / 298.257223563f; // flattening
    //static float eps = Mathf.Sqrt(1 - (b * b) / (a * a)); // first eccentricity
    static float eps_sq = f * (2 - f);  // eccentricity squared
    //static float r_e = 6371f; //Earth's radius

    // Sets EGM-96 model file parameters
    static int NumbGeoidCols = 1441; // 360 degrees of longitude at 15 minute spacing
    static int NumbGeoidRows = 721; // 180 degrees of latitude  at 15 minute spacing
    static int ScaleFactor = 4; // 4 grid cells per degree at 15 minute spacing
    static int NumbGeoidElevs = NumbGeoidCols * NumbGeoidRows;
    static string cofFileName = "WMM.cof";

    static float[] geoidHeights;
    static List<(float, float, float, float)[]> Coeff;
    static float h;
    private Vector3 adj_vector;
    private Vector3 bias;
    private float maxNoiseRange;
    private float3x3 scaleMatrix;

    void Start()
    {
        Coeff = read_Coeff();

        // Access the GeoidHeightBuffer variable from the EGM9615 class
        geoidHeights = EGM9615.GeoidHeightBuffer;

        var mag_elements = Pos_Date2MagElements(lat, lng, alt, "1.1.2024");    //Test

        // Convert magnetic vector orientation to unity
        Vector3 wmm_vector_unity = new Vector3(
            mag_elements[0].Item5,          // East as X
            -mag_elements[0].Item6,         // -Down as Y
            mag_elements[0].Item4          // North as Z
        );

        // Rotate horizontal plane so magnetic north aligns correctly and not with geographic/local ENU North
        Quaternion adj_declination = Quaternion.Euler(0, -mag_elements[0].Item1, 0);
        adj_vector = adj_declination * wmm_vector_unity;

        Debug.Log($"Declination: {mag_elements[0].Item1}°");
        Debug.Log($"N: {wmm_vector_unity.z}nT, E: {wmm_vector_unity.x}nT, D: {-wmm_vector_unity.y}nT");


        // Error Model
        float initialBias = 0.3f; // µT = 1000nT
        float tempCoeff = 0.024f; // Zero-B Offset Temperature drift in µT/°K for -40°C - 85°C
        float T_ref = 25f; // Reference temperature in °C

        float T = 30f; // operating Temperature in °C

        // Calculate temperature-induced bias change
        float T_delta = (T - T_ref);    // delta °C = delta K
        float tempInducedBias = tempCoeff * T_delta;

        // Total bias range including temperature effect per axis and convert µT to nT
        float minBiasRange = (-initialBias - tempInducedBias) / Mathf.Sqrt(3) * 1000;
        float maxBiasRange = (initialBias + tempInducedBias) / Mathf.Sqrt(3) * 1000;

        bias = new Vector3(
            Functions.Random(minBiasRange, maxBiasRange),
            Functions.Random(minBiasRange, maxBiasRange),
            Functions.Random(minBiasRange, maxBiasRange)
        );

        Debug.Log($"Bias: {bias} in nT");

        // Define scale matrix
        scaleMatrix = new float3x3(
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f
        );

        maxNoiseRange = 0.03f * 1000 / Mathf.Sqrt(3); // 0.3µT in 300nT for now, not clear from data sheet

        StartCoroutine(print_magnetometer());
    }


    IEnumerator print_magnetometer()
    {
        while (true)
        {
            // Rotation angles
            Vector3 rotation_angles = multicopter.transform.eulerAngles;

            // Combine the rotation quaternions
            Quaternion rotation = Quaternion.Euler(rotation_angles);

            // Noise
            Vector3 noise = new Vector3(
                Functions.Random(-maxNoiseRange, maxNoiseRange),
                Functions.Random(-maxNoiseRange, maxNoiseRange),
                Functions.Random(-maxNoiseRange, maxNoiseRange)
            );

            // Rotate 
            Vector3 rot_adj_vector = rotation * adj_vector;

            // Include scale error
            float3 scaled_adj_vector = math.mul(scaleMatrix, new float3(rot_adj_vector.x, rot_adj_vector.y, rot_adj_vector.z));

            // Add bias and noise
            Vector3 local_mag_vector = new Vector3(scaled_adj_vector.x, scaled_adj_vector.y, scaled_adj_vector.z) + bias + noise; //noise

            float total_field = (float)Math.Sqrt(local_mag_vector.x * local_mag_vector.x + local_mag_vector.y * local_mag_vector.y + local_mag_vector.z * local_mag_vector.z);
            Debug.Log($"Magnetic Vector in Local Coordinates: {local_mag_vector}; Total field strength: {total_field}");
            yield return new WaitForSeconds(0.25f);
        }
    }


    public static List<(float, float, float, float)[]> read_Coeff()
    {

        int n, m;
        float gnm, hnm, dgnm, dhnm;
        // First line n=0 is not defined, 9999 as filler
        List<(float, float, float, float)[]> c = new List<(float, float, float, float)[]>
        {
            new (float, float, float, float)[] { (9999f, 9999f, 9999f, 9999f) }
        };
        List<(float, float, float, float)> n_grouped = new List<(float, float, float, float)>();

        string cofFilePath = Path.Combine(Application.dataPath, cofFileName);
        try
        {
            // Read all lines from the .cof file
            string[] lines = File.ReadAllLines(cofFilePath);

            // Skip first line (header) 
            string[] skippedLines = lines.Skip(1).ToArray();

            // Extract values of each line and store in coef[n][m]
            foreach (string line in skippedLines)
            {
                // End of File
                if (line == "999999999999999999999999999999999999999999999999")
                {
                    break;
                }

                string[] values = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                int.TryParse(values[0], out n);
                int.TryParse(values[1], out m);

                float.TryParse(values[2], NumberStyles.Any, CultureInfo.InvariantCulture, out gnm);
                float.TryParse(values[3], NumberStyles.Any, CultureInfo.InvariantCulture, out hnm);
                float.TryParse(values[4], NumberStyles.Any, CultureInfo.InvariantCulture, out dgnm);
                float.TryParse(values[5], NumberStyles.Any, CultureInfo.InvariantCulture, out dhnm);

                n_grouped.Add((gnm, hnm, dgnm, dhnm));
                if (n == m)
                {
                    c.Add(n_grouped.ToArray());
                    n_grouped.Clear();
                }
            }
        }
        catch (Exception ex)
        {
            Debug.LogError($"An error occurred: {ex.Message}");
        }
        return c;
    }

    public static float GetGeoidHeight(float lat, float lng)
    {
        int Error_Code = 0;
        int Index;
        float OffsetX, OffsetY;
        int PostX, PostY;
        float ElevationSE, ElevationSW, ElevationNE, ElevationNW;
        float DeltaX, DeltaY, UpperY, LowerY;
        float DeltaHeight;

        if ((lat < -90f) || (lat > 90f))
        {
            Error_Code = 1;
            Debug.Log("Invalid Latitude!");
        }

        if ((lng < -180f) || (lng > 360f))
        {
            Error_Code = 1;
            Debug.Log("Invalid Longitude!");
        }

        if (Error_Code != 1)
        {
            // Compute X and Y Offsets into Geoid Height Array:
            if (lng < 0f)
            {
                OffsetX = (lng + 360f) * ScaleFactor;
            }
            else
            {
                OffsetX = lng * ScaleFactor;
            }
            OffsetY = (90f - lat) * ScaleFactor;

            //  Find Four Nearest Geoid Height Cells for specified Latitude, Longitude;
            //  Assumes that (0,0) of Geoid Height Array is at Northwest corner:

            PostX = (int)Mathf.Floor(OffsetX);
            if ((PostX + 1) == NumbGeoidCols)
            {
                PostX--;
            }

            PostY = (int)Mathf.Floor(OffsetY);
            if ((PostY + 1) == NumbGeoidRows)
            {
                PostY--;
            }

            Index = PostY * NumbGeoidCols + PostX;
            ElevationNW = geoidHeights[Index];
            ElevationNE = geoidHeights[Index + 1];

            Index = (PostY + 1) * NumbGeoidCols + PostX;
            ElevationSW = geoidHeights[Index];
            ElevationSE = geoidHeights[Index + 1];

            //  Perform Bi-Linear Interpolation to compute Height above Ellipsoid:
            DeltaX = OffsetX - PostX;
            DeltaY = OffsetY - PostY;

            UpperY = ElevationNW + DeltaX * (ElevationNE - ElevationNW);
            LowerY = ElevationSW + DeltaX * (ElevationSE - ElevationSW);

            DeltaHeight = UpperY + DeltaY * (LowerY - UpperY);
            return DeltaHeight;

        }
        else
        {
            return -999999;
        }
    }

    public static Vector3 Geodetic2Spherical(float lat, float lng, float h)
    {
        // Convert geodetic coordinates, (defined by the WGS-84 reference ellipsoid), to Earth Centered Earth Fixed Cartesian coordinates, and then to spherical coordinates.
        float cos_phi = Mathf.Cos(Mathf.Deg2Rad * lat);
        float sin_phi = Mathf.Sin(Mathf.Deg2Rad * lat);

        // Compute the local radius of curvature on the WGS-84 reference ellipsoid
        float r_c = a / (Mathf.Sqrt(1 - eps_sq * sin_phi * sin_phi));

        // Compute ECEF Cartesian coordinates of specified point (for longitude=0)
        float p = (r_c + h) * cos_phi;
        float z = (r_c * (1 - eps_sq) + h) * sin_phi;

        // Compute spherical radius and angle phi and lambda of specified point
        float r = Mathf.Sqrt(p * p + z * z);
        float phi_g = Mathf.Rad2Deg * Mathf.Asin(z / r); //geocentric latitude
        // longitude (lambda)
        return new Vector3(lng, phi_g, r);
    }

    public static Vector2 Model_Coeff(int n, int m, float date)
    {
        // reference year is 2020.0, main field coeff in nT, secular variation coeff nT/yr
        float g_nm_t = Coeff[n][m].Item1 + (date - 2020.0f) * Coeff[n][m].Item3;
        float h_nm_t = Coeff[n][m].Item2 + (date - 2020.0f) * Coeff[n][m].Item4;

        return new Vector2(g_nm_t, h_nm_t);
    }

    public static float associatedLegendre(int l, int m, float x)
    {
        //Source: (code slightly changed) "Numerical Recipes: The Art of Scientific Computing" p.254
        //Computes the associated Legendre polynomial P m
        //l(x).Here m and l are integers satisfying
        //0 ≤ m ≤ l, while x lies in the range −1 ≤ x ≤ 1.
        {
            float fact, pmm, pmmp1, somx2;
            int i, ll;
            float pll = 0;

            pmm = 1;
            if (m > 0)
            {
                somx2 = Mathf.Sqrt((1 - x) * (1 + x));
                fact = 1;
                for (i = 1; i <= m; i++)
                {
                    pmm *= -fact * somx2;
                    fact += 2;
                }
            }
            if (l == m)
            {
                return pmm;
            }
            else
            {
                pmmp1 = x * (2 * m + 1) * pmm;
                if (l == (m + 1))
                {
                    return pmmp1;
                }
                else
                {
                    for (ll = m + 2; ll <= l; ll++)
                    {
                        pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                        pmm = pmmp1;
                        pmmp1 = pll;
                    }
                    return pll;
                }
            }
        }
    }

    public static float Factorial(int n)
    {
        if (n == 0 || n == 1)
        {
            return 1;
        }
        float result = 1;
        for (int i = 1; i <= n; i++)
        {
            result *= i;
        }
        return result;
    }

    public static (float, float, float, float, float, float) FieldVec(float lng, float phi_g, float r, float date)
    {
        // Computes the field vector components X', Y', Z' and its secular variations X'_dot, Y_'dot, Z'_dot 
        float sum_X_fv_m_loop, sum_Y_fv_m_loop, sum_Z_fv_m_loop, sum_X_fv_dot_m_loop, sum_Y_fv_dot_m_loop, sum_Z_fv_dot_m_loop;
        float sum_X_fv_n_loop, sum_Y_fv_n_loop, sum_Z_fv_n_loop, sum_X_fv_dot_n_loop, sum_Y_fv_dot_n_loop, sum_Z_fv_dot_n_loop;

        float sin_phi_g = Mathf.Sin(Mathf.Deg2Rad * phi_g);
        float cos_phi_g = Mathf.Cos(Mathf.Deg2Rad * phi_g);
        float tan_phi_g = Mathf.Tan(Mathf.Deg2Rad * phi_g);

        float P_n_m_schmidt, P_nplus_m_schmidt;

        sum_X_fv_n_loop = 0;
        sum_Y_fv_n_loop = 0;
        sum_Z_fv_n_loop = 0;

        sum_X_fv_dot_n_loop = 0;
        sum_Y_fv_dot_n_loop = 0;
        sum_Z_fv_dot_n_loop = 0;

        // for n loop
        for (int n = 1; n <= 12; n++)
        {
            sum_X_fv_m_loop = 0;
            sum_Y_fv_m_loop = 0;
            sum_Z_fv_m_loop = 0;

            sum_X_fv_dot_m_loop = 0;
            sum_Y_fv_dot_m_loop = 0;
            sum_Z_fv_dot_m_loop = 0;

            for (int m = 0; m <= n; m++)
            {
                // dynamic start (n, m)
                float g_nm_t = Model_Coeff(n, m, date).x;
                float h_nm_t = Model_Coeff(n, m, date).y;
                float g_nm_dot = Coeff[n][m].Item3;
                float h_nm_dot = Coeff[n][m].Item4;
                float cos_m_lam = Mathf.Cos(m * Mathf.Deg2Rad * lng);
                float sin_m_lam = Mathf.Sin(m * Mathf.Deg2Rad * lng);

                // Calculating P_n,m(sin(phi_g)) using the relation P_n,m(µ) = (-1)^(m)*P_(n)^(m)(µ) >> µ will be substituted with sin(phi_g) not cos()! Calling associated Legendre
                float P_n_m = (float)Math.Pow(-1, m) * associatedLegendre(n, m, sin_phi_g);
                float P_nplus_m = (float)Math.Pow(-1, m) * associatedLegendre(n + 1, m, sin_phi_g);


                // Schmidt semi-normalized associated Legendre functions
                if (m == 0)
                {
                    P_n_m_schmidt = P_n_m;
                    P_nplus_m_schmidt = P_nplus_m;
                }
                else
                {
                    P_n_m_schmidt = Mathf.Sqrt((float)(2 * Factorial(n - m) / Factorial(n + m))) * P_n_m;
                    P_nplus_m_schmidt = Mathf.Sqrt((float)(2 * Factorial(n + 1 - m) / Factorial(n + 1 + m))) * P_nplus_m;
                }

                // Derivative
                float dP_n_m_schmidt = (n + 1) * tan_phi_g * P_n_m_schmidt - Mathf.Sqrt((n + 1) * (n + 1) - m * m) * (1 / cos_phi_g) * P_nplus_m_schmidt;

                // Field Vector components in geocentric coordinates
                sum_X_fv_m_loop += (g_nm_t * cos_m_lam + h_nm_t * sin_m_lam) * dP_n_m_schmidt;
                sum_Y_fv_m_loop += m * (g_nm_t * sin_m_lam - h_nm_t * cos_m_lam) * P_n_m_schmidt;
                sum_Z_fv_m_loop += (g_nm_t * cos_m_lam + h_nm_t * sin_m_lam) * P_n_m_schmidt;

                // Secular variation
                sum_X_fv_dot_m_loop += (g_nm_dot * cos_m_lam + h_nm_dot * sin_m_lam) * dP_n_m_schmidt;
                sum_Y_fv_dot_m_loop += m * (g_nm_dot * sin_m_lam - h_nm_dot * cos_m_lam) * P_n_m_schmidt;
                sum_Z_fv_dot_m_loop += (g_nm_dot * cos_m_lam + h_nm_dot * sin_m_lam) * P_n_m_schmidt;
            }
            sum_X_fv_n_loop += Mathf.Pow(a / r, n + 2) * sum_X_fv_m_loop;
            sum_Y_fv_n_loop += Mathf.Pow(a / r, n + 2) * sum_Y_fv_m_loop;
            sum_Z_fv_n_loop += (n + 1) * Mathf.Pow(a / r, n + 2) * sum_Z_fv_m_loop;

            sum_X_fv_dot_n_loop += Mathf.Pow(a / r, n + 2) * sum_X_fv_dot_m_loop;
            sum_Y_fv_dot_n_loop += Mathf.Pow(a / r, n + 2) * sum_Y_fv_dot_m_loop;
            sum_Z_fv_dot_n_loop += (n + 1) * Mathf.Pow(a / r, n + 2) * sum_Z_fv_dot_m_loop;
        }
        // Returning X_fv, Y_fv, Z_fv, X_fv_dot, Y_fv_dot, Z_fv_dot
        return (-sum_X_fv_n_loop, (1 / cos_phi_g) * sum_Y_fv_n_loop, -sum_Z_fv_n_loop, -sum_X_fv_dot_n_loop, (1 / cos_phi_g) * sum_Y_fv_dot_n_loop, -sum_Z_fv_dot_n_loop);
    }

    public static (float, float, float, float, float, float) geocen2ellip((float, float, float, float, float, float) FieldVec, float lat, float phi_g)
    {
        float X_fv = FieldVec.Item1;
        float Y_fv = FieldVec.Item2;
        float Z_fv = FieldVec.Item3;
        float X_fv_dot = FieldVec.Item4;
        float Y_fv_dot = FieldVec.Item5;
        float Z_fv_dot = FieldVec.Item6;

        float phi_diff = phi_g - lat; // Difference of latitude in geocentric and geodetic

        float sin_phi_diff = Mathf.Sin(Mathf.Deg2Rad * phi_diff);
        float cos_phi_diff = Mathf.Cos(Mathf.Deg2Rad * phi_diff);

        // Main Vector
        float X = X_fv * cos_phi_diff - Z_fv * sin_phi_diff;
        float Y = Y_fv;
        float Z = X_fv * sin_phi_diff + Z_fv * cos_phi_diff;

        // Derivatives
        float X_dot = X_fv_dot * cos_phi_diff - Z_fv_dot * sin_phi_diff;
        float Y_dot = Y_fv_dot;
        float Z_dot = X_fv_dot * sin_phi_diff + Z_fv_dot * cos_phi_diff;

        return (X, Y, Z, X_dot, Y_dot, Z_dot);
    }

    public static (float, float, float, float, float, float, float, float) MagElements((float, float, float, float, float, float) comp)
    {
        float X = comp.Item1;
        float Y = comp.Item2;
        float Z = comp.Item3;
        float X_dot = comp.Item4;
        float Y_dot = comp.Item5;
        float Z_dot = comp.Item6;

        // Magnetic Elements
        float H = Mathf.Sqrt(X * X + Y * Y);
        float F = Mathf.Sqrt(H * H + Z * Z);
        float I = Mathf.Atan2(Z, H);    // Angles in rad, usual output in deg
        float D = Mathf.Atan2(Y, X);    //

        // Secular Variation
        float H_dot = (X * X_dot + Y * Y_dot) / H;
        float F_dot = (X * X_dot + Y * Y_dot + Z * Z_dot) / F;
        float I_dot = (H * Z_dot - Z * H_dot) / (F * F);
        float D_dot = (X * Y_dot - Y * X_dot) / (H * H);
        // Angles in rad, need to be converted to deg

        return (H, F, I * Mathf.Rad2Deg, D * Mathf.Rad2Deg, H_dot, F_dot, I_dot * Mathf.Rad2Deg, D_dot * Mathf.Rad2Deg);
    }

    public static float Date2decYear(string date)
    {
        //input i.e. 15.1.2024
        int year, month, day, extraday;
        float total_days;

        string[] values = date.Split('.');

        int.TryParse(values[0], out day);
        int.TryParse(values[1], out month);
        int.TryParse(values[2], out year);

        if ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0)
        {
            extraday = 1;
            total_days = 366;
        }
        else
        {
            extraday = 0;
            total_days = 365;
        }

        List<int> DaysinMonth = new List<int> { 31, 28 + extraday, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

        int sum_days = day; //add days from previous months to this

        for (int i = 0; i < (month - 1); i++)
        {
            sum_days += DaysinMonth[i];
        }

        float dec_year = year + sum_days / total_days;

        return dec_year;
    }

    public static List<(float, float, float, float, float, float, float)> Pos_Date2MagElements(float lat, float lng, float alt, string date)
    {

        float dec_year = Date2decYear(date);

        h = GetGeoidHeight(lat, lng); // Height above the WGS 84 Ellipsoid in m

        Vector3 spherical = Geodetic2Spherical(lat, lng, h);

        var f_v = FieldVec(spherical.x, spherical.y, spherical.z, dec_year);

        var rotated_comp = geocen2ellip(f_v, lat, spherical.y);     //X, Y, Z, X_dot, Y_dot, Z_dot

        var magElements = MagElements(rotated_comp);    // H, F, I, D, H_dot, F_dot, I_dot, D_dot

        Debug.Log($"Latitude: {lat} | Longitude: {lng} | Altitude: {alt} (Height above Mean Sea Level)");
        Debug.Log("Declination D  |  Inclination I  | Horizontal Intensity H | North Comp X   | East Comp Y | Down Comp Z  | Total Field F");
        Debug.Log($"{Math.Round(magElements.Item4, 4)}°             | {Math.Round(magElements.Item3, 4)}°          | {Math.Round(magElements.Item1, 4)}nT            | {Math.Round(rotated_comp.Item1, 4)}nT   | {Math.Round(rotated_comp.Item2, 4)}nT  | {Math.Round(rotated_comp.Item3, 4)}nT  | {Math.Round(magElements.Item2, 4)}nT");
        Debug.Log($"{Math.Round(magElements.Item8, 4)}°/yr        | {Math.Round(magElements.Item7, 4)}°/yr    | {Math.Round(magElements.Item5, 4)}nT/yr                  | {Math.Round(rotated_comp.Item4, 4)}nT/yr     | {Math.Round(rotated_comp.Item5, 4)}nT/yr | {Math.Round(rotated_comp.Item6, 4)}nT/yr       | {Math.Round(magElements.Item6, 4)}nT/yr");


        // returned in the format of WMM Calculator from NOAA
        // D I H X Y Z F (G grid actually not)
        // changes/year
        List<(float, float, float, float, float, float, float)> mag_elements = new List<(float, float, float, float, float, float, float)>();

        mag_elements.Add(((float)Math.Round(magElements.Item4, 4), (float)Math.Round(magElements.Item3, 4), (float)Math.Round(magElements.Item1, 4), (float)Math.Round(rotated_comp.Item1, 4), (float)Math.Round(rotated_comp.Item2, 4), (float)Math.Round(rotated_comp.Item3, 4), (float)Math.Round(magElements.Item2, 4)));
        mag_elements.Add(((float)Math.Round(magElements.Item8, 4), (float)Math.Round(magElements.Item7, 4), (float)Math.Round(magElements.Item5, 4), (float)Math.Round(rotated_comp.Item4, 4), (float)Math.Round(rotated_comp.Item5, 4), (float)Math.Round(rotated_comp.Item6, 4), (float)Math.Round(magElements.Item6, 4)));

        return mag_elements;
    }

    // measure the angle between main axis X and Z
    // rotate by declination
    //
}

// Transform from global to local ENU
// Include Error Model / uncertainty?


