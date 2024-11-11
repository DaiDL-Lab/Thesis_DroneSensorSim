using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class Gyroscope : MonoBehaviour
{
    public GameObject multicopter;
    private Vector3 lastEulerAngles;
    private Vector3 rotationRate;
    private float maxNoiseRange;
    private Vector3 bias;
    private float3x3 scaleMatrix;


    // Start is called before the first frame update
    void Start()
    {
        lastEulerAngles = multicopter.transform.eulerAngles;

        float bandwidth = 100; //in HZ, certain bandwidth depending on sensor model
        float T = 30f; // operating Temperature in °C

        // Depending on which sensor model is being used
        #region ICM-42688-P configuration

        float RNSD = 0.0028f; // Rate Noise Spectral Density in (°/s)/(√Hz)
        float initialBias = 0.5f; // Initial ZRO Tolerance in º/s
        float tempCoeff = 0.005f; // ZRO Variation vs. Temperature in º/s/°C for 0°C 70°C
        float T_ref = 25f; // Reference temperature in °C
        #endregion

        #region BMI055 configuration
        /*
        float RNSD = 0.014f; // Rate Noise Spectral Density in (°/s)/(√Hz)
        float initialBias = 1f; // Initial ZRO Tolerance in º/s
        float tempCoeff = 0.015f; // ZRO Variation vs. Temperature in º/s/K for -40°C 85°C
        float T_ref = 25f; // Reference temperature in °C
        */
        #endregion

        // Calculate temperature-induced bias change
        float T_delta = (T - T_ref);    // delta °C = delta K
        float tempInducedBias = tempCoeff * T_delta;

        // Total bias range including temperature effect per axis
        float minBiasRange = (-initialBias - tempInducedBias) / Mathf.Sqrt(3);
        float maxBiasRange = (initialBias + tempInducedBias) / Mathf.Sqrt(3);

        bias = new Vector3(
            Functions.Random(minBiasRange, maxBiasRange),
            Functions.Random(minBiasRange, maxBiasRange),
            Functions.Random(minBiasRange, maxBiasRange)
        );

        maxNoiseRange = RNSD * Mathf.Sqrt(bandwidth) / Mathf.Sqrt(3); // Calculate RMS noise for each axis

        // Define scale matrix
        scaleMatrix = new float3x3(
            new float3(1.0f, 0.0f, 0.0f),
            new float3(0.0f, 1.0f, 0.0f),
            new float3(0.0f, 0.0f, 1.0f)
        );

        StartCoroutine(print_gyroscope());
    }

    IEnumerator print_gyroscope()
    {
        while (true)
        {
            // Calculate the raw rotation rate
            Vector3 currentEulerAngles = multicopter.transform.eulerAngles;
            Vector3 raw_rotationRate = (currentEulerAngles - lastEulerAngles) / Time.deltaTime;

            // Noise
            Vector3 noise = new Vector3(
                Functions.Random(-maxNoiseRange, maxNoiseRange),
                Functions.Random(-maxNoiseRange, maxNoiseRange),
                Functions.Random(-maxNoiseRange, maxNoiseRange)
            );

            // scaled_rotationRate, converted to float3 for matrix multiplication with float3x3
            float3 s_rR = math.mul(scaleMatrix, new float3(raw_rotationRate.x, raw_rotationRate.y, raw_rotationRate.z));

            // Calculate final rotation rate incl. scale error, bias and noise
            rotationRate = new Vector3(s_rR.x, s_rR.y, s_rR.z) + bias + noise;

            // Display the tilting angles and rotation rate
            Vector3 tiltAngles = multicopter.transform.eulerAngles;
            Debug.Log($"Tilt Angles: {tiltAngles}, Rotation Rate: X-axis: {rotationRate.x} °/s, Y-axis: {rotationRate.y} °/s, Z-axis: {rotationRate.z} °/s");

            lastEulerAngles = currentEulerAngles;   // updating last Angles
            yield return new WaitForSeconds(0.25f);
        }
    }
}
