using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using MathNet.Numerics;
using MathNet.Numerics.Statistics; // Add this line

public class Noise_Testing : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        float min = -0.06867f; //m/^2 equal to +/- 70mg (milli-earthacceleration)
        float max = 0.06867f; // m/^2
        List<float> p = new List<float>();
        for (int i = 0; i < 13034; i++)
        {
            p.Add(Functions.Random(min, max));
        }

        double mean = Statistics.Mean(p);
        double stdDev = Statistics.StandardDeviation(p);

        Debug.Log($"Average noise: {mean} | Standard Devitatio of noise: {stdDev}");

    }
}
