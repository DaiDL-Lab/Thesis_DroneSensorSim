using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Functions : MonoBehaviour
{
    public static float Random(float minValue, float maxValue)
    {
        // White noise as normal distribution
        //UnityEngine.Random.InitState(seed); // Initialize random state with the provided seed

        // Define mean as the midpoint of the range
        float mean = (minValue + maxValue) / 2.0f;

        // Define stdDev as a fraction of the range width (for example, 1/6th of the range)
        float stdDev = (maxValue - minValue) / 6.0f;

        // Generate two uniform random numbers in the range (0, 1]
        float u1 = 1.0f - UnityEngine.Random.value;
        float u2 = 1.0f - UnityEngine.Random.value;

        // Apply Box-Muller transform
        float randStdNormal = Mathf.Sqrt(-2.0f * Mathf.Log(u1)) * Mathf.Sin(2.0f * Mathf.PI * u2);

        // Return the random number with specified mean and standard deviation
        float randomValue = mean + stdDev * randStdNormal;

        // Clamp the value to ensure it stays within the specified range
        return Mathf.Clamp(randomValue, minValue, maxValue);
    }

    public static double MapValue(double variable, double lowerLimit1, double upperLimit1, double lowerLimit2, double upperLimit2)
    {
        // Ensure the variable is within the source range
        if (variable < lowerLimit1 || variable > upperLimit1)
        {
            throw new ArgumentOutOfRangeException(nameof(variable), "Variable is out of the source range.");
        }

        // Perform the mapping
        double mappedValue = ((variable - lowerLimit1) / (upperLimit1 - lowerLimit1)) * (upperLimit2 - lowerLimit2) + lowerLimit2;

        return mappedValue;
    }

}