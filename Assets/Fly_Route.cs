using System;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

public class Fly_Route : MonoBehaviour
{
    public string filePath; // Path to the CSV file
    public GameObject targetObject; // The GameObject to move
    public float moveSpeed = 1.0f; // Speed of movement

    private List<Vector3> positions = new List<Vector3>();
    private int currentIndex = 0;

    #region EventHandler
    private void OnEnable()
    {
        SatelliteSetup.SetupDone += Fly;
    }

    private void OnDisable()
    {
        SatelliteSetup.SetupDone -= Fly;
    }
    #endregion

    void Fly()
    {
        // Read the CSV file
        ReadCSV(filePath);

        // Start the coroutine to move the object along the path
        StartCoroutine(MoveAlongPath());
    }

    void ReadCSV(string filePath)
    {
        using (StreamReader reader = new StreamReader(filePath))
        {
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                string[] values = line.Split(',');

                if (values.Length == 3)
                {
                    float x = float.Parse(values[0]);
                    float y = float.Parse(values[1]);
                    float z = float.Parse(values[2]);

                    Vector3 position = new Vector3(x, y, z);
                    positions.Add(position);
                }
            }
        }
    }

    System.Collections.IEnumerator MoveAlongPath()
    {
        while (currentIndex < positions.Count)
        {
            Vector3 targetPosition = positions[currentIndex];

            // Move the object towards the target position
            while (targetObject.transform.position != targetPosition)
            {
                targetObject.transform.position = Vector3.MoveTowards(targetObject.transform.position, targetPosition, moveSpeed * Time.deltaTime);
                yield return null; // Wait for the next frame
            }

            currentIndex++; // Move to the next position
        }
    }
}
