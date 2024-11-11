using System.Collections;
using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;

public class Accelerometer : MonoBehaviour
{
    public GameObject multicopter;
    private Rigidbody rb;
    private Vector3 lastVelocity;
    private Vector3 acceleration;
    private float maxNoiseRange_XY;
    private float maxNoiseRange_Z;
    private Vector3 bias;
    private float3x3 scaleMatrix;


    private void Start()
    {
        lastVelocity = multicopter.GetComponent<Rigidbody>().velocity;
        rb = multicopter.GetComponent<Rigidbody>();

        float g = 9.81f; // earth acceleration in m/s^2
        float bandwidth = 100f; //in HZ, certain bandwidth depending on sensor model
        float T = 30f; // operating Temperature in °C

        // Depending on which sensor model is being used
        #region ICM-42688-P configuration
        
        float PSD_XY = 65f; // Power Spectral Density for X and Y-Axis in µg/√Hz, g is earth acceleration, not gramm
        float PSD_Z = 70f;  // for Z-Axis
        float initialBias = 20f; // Initial Zero-G Tolerance in mg (milli g-acceleration)
        float tempCoeff = 0.15f; // Zero-G Level Change vs. Temperature in mg/°C
        float T_ref = 25f; // Reference temperature in °C
        
        #endregion

        #region BMI055 configuration
        /*
        float PSD_XY = 150f; // Power Spectral Density/Output Noise Density in µg/√Hz, X=Y=Z no seperation
        float PSD_Z = PSD_XY;  // for Z-Axis
        float initialBias = 70f; // Initial Zero-G Tolerance in mg (milli g-acceleration) Off_x,z and Off_y are equal
        float tempCoeff = 1f; // Zero-G Level Change vs. Temperature in mg/K
        float T_ref = 25f; // Reference temperature in °C
        */
        #endregion

        // Calculate temperature-induced bias change
        float T_delta = (T - T_ref); // delta °C = delta K
        float tempInducedBias = tempCoeff * T_delta;

        // Total bias range including temperature effect convert mg in m/s^2 per Axis
        float minBiasRange = ((-initialBias - tempInducedBias) / 1000 * g) / Mathf.Sqrt(3);
        float maxBiasRange = ((initialBias + tempInducedBias) / 1000 * g) / Mathf.Sqrt(3);

        bias = new Vector3(
            Functions.Random(minBiasRange, maxBiasRange),
            Functions.Random(minBiasRange, maxBiasRange),
            Functions.Random(minBiasRange, maxBiasRange)
        );

        maxNoiseRange_XY = PSD_XY * Mathf.Sqrt(bandwidth) / 1000000 * g; // convert µg to g, then insert g for m/s^2
        maxNoiseRange_Z = PSD_Z * Mathf.Sqrt(bandwidth) / 1000000 * g;

        // Define scale matrix
        scaleMatrix = new float3x3(
            1.0f, 0.0f, 0.0f,
            0.0f, 1.0f, 0.0f,
            0.0f, 0.0f, 1.0f
        );

        StartCoroutine(print_accelerometer());
    }

    IEnumerator print_accelerometer()
    {
        while (true) // Loop indefinitely
        {
            // Calculate the raw acceleration
            Vector3 currentVelocity = rb.velocity;  // aka speed in m/s
            Vector3 raw_Acceleration = (currentVelocity - lastVelocity) / Time.deltaTime; // in m/s^2

            // Add noise
            Vector3 noise = new Vector3(
                Functions.Random(-maxNoiseRange_XY, maxNoiseRange_XY),
                Functions.Random(-maxNoiseRange_XY, maxNoiseRange_XY),
                Functions.Random(-maxNoiseRange_Z, maxNoiseRange_Z)
            );

            // scaledAcceleration, converted to float3 for matrix multiplication with float3x3
            float3 s_acc = math.mul(scaleMatrix, new float3(raw_Acceleration.x, raw_Acceleration.y, raw_Acceleration.z));

            // Calculate final acceleration incl. scale error, bias and noise
            acceleration = new Vector3(s_acc.x, s_acc.y, s_acc.z) + bias + noise; //bias, gravity?

            // Display the acceleration
            Debug.Log($"Acceleration: X: {acceleration.x} m/s^2, Y: {acceleration.y} m/s^2, Z: {acceleration.z} m/s^2");

            lastVelocity = currentVelocity;
            yield return new WaitForSeconds(0.25f);
        }
    }
}
