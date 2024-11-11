using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class Barometer : MonoBehaviour
{

    private const float P0 = 101325f; // Sea level standard atmospheric pressure in Pa
    private const float L = 0.0065f; // Temperature lapse rate in K/m
    private const float T0 = 288.15f; // MSL standard temperature in K
    private const float g = 9.80665f; // Earth's gravitational acceleration in m/s^2
    private const float M = 0.02896f; // Molar mass of Earth's air in kg/mol
    private const float R = 8.314f; // Universal gas constant in J/(mol·K)
    //private readonly float bias; // Bias in Pa
    private float maxNoiseRange = 2f; // Maximum noise range in Pa, derived from ±2.5 mbar errorband
    private float maxBiasRange = 0.5f;  // Maximum bias range in Pa, derived from accuracy ±1.5 mbar
    private float bias;

    public GameObject Multicopter;

    void Start()
    {
        bias = Functions.Random(-maxBiasRange, maxBiasRange) * 100;   // converted mbar to Pa
        //float drift = 0f; // Initialize drift
        StartCoroutine(print_barometer());
    }

    IEnumerator print_barometer()
    {
        while (true)
        {
            float alt = 42.0f + Multicopter.transform.position.y; //Height to Local Origin + Local Multicopter Height

            float noise = Functions.Random(-maxNoiseRange, maxNoiseRange) * 100; //range converted from mbar in Pa

            // Calculate pressure for Berlin
            float p = pressure(alt, T0) + bias + noise; // calculation + errors
            Debug.Log($"Berlin: Altitude: {alt} meters, Atmospheric Pressure: {p} Pa / {p / 100} mbar");

            yield return new WaitForSeconds(0.25f);
        }
    }

    // Function to calculate atmospheric pressure at a given altitude above sea level
    private static float pressure(float h_asl, float T0)
    {
        float pressure = P0 * Mathf.Pow((1 - (L * h_asl) / T0), (g * M) / (R * L)); // works for tropospheric height <11km
        return pressure;
    }
}
