#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <time.h>
using namespace std;

struct constraint
{
    int a; // An entity
    int b; // Another entity
    int r; // Constraint between both entities
};

// Constraints list computing
vector<constraint> calculate_const_list(vector<vector<int>> *cl_set_const)
{

    vector<constraint> const_list; // Contraints list
    constraint r;                  // A constraint

    for (int i = 0; i < (*cl_set_const).size(); ++i)
    {
        for (int j = i; j < (*cl_set_const)[i].size(); ++j)
        {
            if ((*cl_set_const)[i][j] == 1 || (*cl_set_const)[i][j] == -1)
            {
                r.a = i;
                r.b = j;
                r.r = (*cl_set_const)[i][j];
                const_list.push_back(r);
            }
        }
    }

    return const_list;
}

// Infeasibility computing
int calculate_infeasibility(vector<int> sol, vector<constraint> const_list)
{

    int infeasibility = 0; // Infeasibility value

    for (int i = 0; i < const_list.size(); ++i)
    {
        if (const_list[i].r == 1)
        {
            if (sol[const_list[i].a] != sol[const_list[i].b])
                ++infeasibility;
        }
        else
        {
            if (sol[const_list[i].a] == sol[const_list[i].b])
                ++infeasibility;
        }
    }

    return infeasibility;
}

// General deviation computing
float calculate_desv(vector<vector<float>> *cl_set, vector<int> sol, int k)
{

    vector<vector<float>> centroids; // Vector of centroids
    vector<float> fl_v_aux;          // Float vector assistant
    float fl_aux;                    // Float assistant
    vector<float> cl_counter;        // Clusters entities counter
    vector<float> intracl_dist;      // Clusters mean intra-cluster distance values
    float general_deviation = 0;     // General deviation value

    // Centroids computing

    for (int i = 0; i < (*cl_set)[0].size(); i++)
        fl_v_aux.push_back(0);
    for (int i = 0; i < k; i++)
        centroids.push_back(fl_v_aux);

    for (int i = 0; i < centroids.size(); i++)
        cl_counter.push_back(0);
    for (int i = 0; i < sol.size(); i++)
        cl_counter[sol[i]]++;

    for (int i = 0; i < sol.size(); i++)
        for (int j = 0; j < centroids[0].size(); j++)
            centroids[sol[i]][j] += (*cl_set)[i][j];

    for (int i = 0; i < centroids.size(); i++)
        for (int j = 0; j < centroids[0].size(); j++)
            centroids[i][j] /= cl_counter[i];

    // Mean intra-cluster distance computing

    for (int i = 0; i < k; i++)
        intracl_dist.push_back(0);

    for (int i = 0; i < sol.size(); i++)
    {
        fl_aux = 0;
        for (int j = 0; j < centroids[0].size(); j++)
            fl_aux += pow(centroids[sol[i]][j] - (*cl_set)[i][j], 2);
        intracl_dist[sol[i]] += pow(fl_aux, 0.5);
    }
    for (int i = 0; i < k; i++)
        intracl_dist[i] /= cl_counter[i];

    //General desviation computing

    for (int i = 0; i < k; i++)
        general_deviation += intracl_dist[i];

    return general_deviation /= k;
}

// Lambda value computing
float calculate_lambda(vector<vector<float>> *cl_set, int const_num)
{
    float max_distance = 0; // Maximum distance between entities
    float fl_aux;           // Float assistant

    for (int i = 0; i < (*cl_set).size(); ++i)
    {
        for (int j = i + 1; j < (*cl_set).size(); ++j)
        {
            fl_aux = 0;
            for (int k = 0; j < (*cl_set)[i].size(); ++j)
            {
                fl_aux += pow((*cl_set)[i][k] - (*cl_set)[j][k], 2);
            }
            fl_aux = pow(fl_aux, 0.5);
            if (fl_aux > max_distance)
                max_distance = ceil(fl_aux);
        }
    }
    max_distance /= const_num;

    return max_distance;
}

// Cost function f computing
float calculate_f(vector<vector<float>> *cl_set, vector<int> sol, vector<constraint> const_list, int k, float lambda, float *deviation, int *infeas)
{
    bool empty_condition = true; // True if there are no empty clusters, false otherwise
    int infeasibility;           // Infeasibility value
    float general_deviation;     // General deviation value

    //Empty clusters control

    for (int i = 0; i < k; i++)
    {
        if (count(sol.begin(), sol.end(), i) < 1)
        {
            empty_condition = false;
            break;
        }
    }

    //Infeasibility computing

    if (empty_condition)
        infeasibility = calculate_infeasibility(sol, const_list);
    else
        infeasibility = 9999;

    (*infeas) = infeasibility;

    general_deviation = calculate_desv(cl_set, sol, k);

    (*deviation) = general_deviation;

    return (general_deviation + infeasibility * lambda);
}

void silly_season(vector<vector<float>> *cl_set, vector<vector<int>> *cl_set_const, int k, int seed)
{
    time_t time1, time2;    // Time measure
    float best_dev;         // Best general deviation
    int best_inf;           // Best infeasibility
    int best_f;             // Best cost value

    time1 = clock();

    srand(seed);

    time2 = clock();

    cout << best_dev << " ";
    cout << best_inf << " ";
    cout << best_f << " ";
    cout << (float)time2 / CLOCKS_PER_SEC - (float)time1 / CLOCKS_PER_SEC << endl;
}

int main(int argc, char *argv[])
{

    //Read data

    string data_path = "data/";  // Database file path
    string const_path = "data/"; // Constraints file path
    string string_aux;           // String assistant
    int k;                       // Number of clusters
    int seed;                    // Random seed

    if (argc == 5)
    {
        data_path += argv[1];
        const_path += argv[2];
        k = atoi(argv[3]);
        seed = atoi(argv[4]);
    }
    else
    {
        cout << "\nWrong number of arguments, must be entered below" << endl;
        cout << "Enter the database file name: ";
        cin >> string_aux;
        data_path += string_aux;
        cout << "Enter the constraints file name: ";
        cin >> string_aux;
        const_path += string_aux;
        cout << "Enter the clusters number: ";
        cin >> k;
        cout << "Enter the random seed: ";
        cin >> seed;
    }

    if (k < 1)
    {
        cout << "Clusters number must be positive" << endl;
        return 0;
    }

    // Dimension computing

    ifstream f(data_path); // Data stream
    if (!f)
    {
        cout << "Error opening database file" << endl;
        return 0;
    }

    int dimension = 0; // Dimensions number
    float reading;     // Read value
    char trash = ',';  // File separator

    while (trash == ',')
    {
        dimension++;
        f >> reading;
        f >> trash;
    }

    f.seekg(0);

    // Database reading

    vector<vector<float>> cl_set; // Data set
    vector<float> tuple_cl;       // Read element

    while (!f.eof())
    {
        tuple_cl.clear();
        for (int i = 0; i < dimension; i++)
        {
            f >> reading;
            f.ignore();
            tuple_cl.push_back(reading);
        }
        cl_set.push_back(tuple_cl);
    }

    f.close();

    // Constraints reading

    vector<vector<int>> cl_set_const; // Constraints set
    vector<int> tuple_const;          // Read constraint
    int const_reading;                // Read value

    f.open(const_path);

    if (!f)
    {
        cout << "Error opening constraints file" << endl;
        return 0;
    }

    for (int i = 0; i < cl_set.size(); i++)
    {
        tuple_const.clear();
        for (int j = 0; j < cl_set.size(); j++)
        {
            f >> const_reading;
            f.ignore();
            tuple_const.push_back(const_reading);
        }
        cl_set_const.push_back(tuple_const);
    }

    f.close();

    // Silly season algorithm aplication

    silly_season(&cl_set,&cl_set_const,k,seed);

    return 0;
}