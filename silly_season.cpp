//////////////////////////////////////////////////////////////////////////////////
//  SILLY SEASON ALGORITHM implementation to constrained clustering problem     //
//                                                                              //
//  Antonio José Blánquez Pérez - 3ºCSI                                         //
//  FINAL PRACTICE - METAHEURÍSTICAS - University of Granada                    //
//////////////////////////////////////////////////////////////////////////////////

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

// Standard deviation computing
float calculate_desv(vector<vector<float>> *cl_set, vector<int> sol, int k)
{

    vector<vector<float>> centroids; // Vector of centroids
    vector<float> fl_v_aux;          // Auxiliary float vector
    float fl_aux;                    // Auxiliary float
    vector<float> cl_counter;        // Clusters entities counter
    vector<float> intracl_dist;      // Clusters mean intra-cluster distance values
    float standard_deviation = 0;    // Standard deviation value

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

    //Standard desviation computing

    for (int i = 0; i < k; i++)
        standard_deviation += intracl_dist[i];

    return standard_deviation /= k;
}

// Lambda value computing
float calculate_lambda(vector<vector<float>> *cl_set, int const_num)
{
    float max_distance = 0; // Maximum distance between entities
    float fl_aux;           // Auxiliary float

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
    float standard_deviation;    // Standard deviation value

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

    standard_deviation = calculate_desv(cl_set, sol, k);

    (*deviation) = standard_deviation;

    return (standard_deviation + infeasibility * lambda);
}

// Silly season algorithm
void silly_season(vector<vector<float>> *cl_set, vector<vector<int>> *cl_set_const, int k, int seed)
{
    ////// PARAMETERS ////////////////////////////////////////////////////////////////////////////
    int n_teams = 5;                     // Number of teams
    int n_pilot_per_team = 2;            // Number of pilots per team
    int season_lenght = 10000;           // Season lenght in evaluations
    int max_iters = 100000;              // Maximum number of evaluations
    int max_iters_bl = 1000;             // Maximum number of evaluations for BL
    int prob_mutation = 1050;            // Mutation probability
    int n_team_changes_silly_season = 4; // Number of team changes between seasons
    int n_team_changes_mid_season = 1;   // Number of team changes during a season
    //////////////////////////////////////////////////////////////////////////////////////////////

    time_t time1, time2;            // Time measure
    float best_dev;                 // Best standard deviation
    int best_inf;                   // Best infeasibility
    float best_f;                   // Best cost value
    int n_iters = 0;                // Number of evaluations
    int current_season = 1;         // Current season
    int last_season = 1;            // Last season
    float current_f;                // Current f value
    int current_inf;                // Current infeasibility
    float current_dev;              // Current standard deviation
    vector<vector<float>> f_matrix; // Each team f value
    vector<float> team_hierarchy;   // Team hierarchy representation
    float mean_f;                   // Cut-off point

    float f_act;                    // Actual cost function value for BL
    float f_ini;                    // Best cost function value for BL
    int bl_counter;                 // BL evaluations counter
    vector<vector<int>> neighbours; // Neighbours vector for BL

    int n_mutations = 0; // Number of mutations
    int ind, val;        // Auxiliary variables for mutations

    int n_team_changes;         // Number of team changes
    int ind1, val1, ind2, val2; // Auxiliary variables for team changes

    vector<int> c;                     // Auxiliary vector of int
    vector<float> cf;                  // Auxiliary vector of float
    vector<vector<int>> t;             // Auxiliary vector of vector of int
    vector<vector<vector<int>>> teams; // Teams set
    float i_aux;                       // Auxiliary float
    bool condition;                    // Empty clusters condition
    time1 = clock();

    srand(seed);

    // Initialize teams randomly

    for (int m = 0; m < n_teams; m++)
    {
        t.clear();
        for (int j = 0; j < n_pilot_per_team; j++)
        {
            do
            {
                c.clear();
                for (int i = 0; i < (*cl_set).size(); i++)
                    c.push_back(rand() % k);
                condition = count(c.begin(), c.end(), 0) == 0;
                for (int i = 1; i < k && !condition; i++)
                    condition || count(c.begin(), c.end(), i) == 0;
            } while (condition);
            t.push_back(c);
        }
        teams.push_back(t);
    }

    // Show teams

    /*for(int i=0;i<teams.size();i++){
        for(int j=0;j<n_pilot_per_team;j++){
            for(int m=0;m<c.size();m++){
                cout << teams[i][j][m] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }*/

    // Common values computing

    vector<constraint> const_list = calculate_const_list(cl_set_const);
    float lambda = calculate_lambda(cl_set, const_list.size());

    // Initialize best values

    best_f = calculate_f(cl_set, teams[0][0], const_list, k, lambda, &best_dev, &best_inf);
    n_iters++;

    f_matrix.clear();
    for (int i = 0; i < n_teams; i++)
    {
        cf.clear();
        for (int j = 0; j < n_pilot_per_team; j++)
        {
            current_f = calculate_f(cl_set, teams[i][j], const_list, k, lambda, &current_dev, &current_inf);
            n_iters++;
            if (current_f < best_f)
            {
                best_f = current_f;
                best_inf = current_inf;
                best_dev = current_dev;
            }
            cf.push_back(current_f);
        }
        f_matrix.push_back(cf);
    }

    // Initialize teams hierarchy

    mean_f = 0;
    team_hierarchy.clear();
    for (int i = 0; i < n_teams; i++)
    {
        i_aux = 0;
        for (int j = 0; j < n_pilot_per_team; j++)
        {
            i_aux += f_matrix[i][j];
            mean_f += f_matrix[i][j];
        }
        team_hierarchy.push_back(i_aux / n_pilot_per_team);
    }
    mean_f = mean_f / (n_teams * n_pilot_per_team);

    // Main loop

    do
    {
        // Each team train its pilots according to its hierarchy

        for (int i = 0; i < n_teams; i++)
        {
            if (team_hierarchy[i] > mean_f) // BL
            {
                for (int j = 0; j < n_pilot_per_team; j++)
                {
                    bl_counter = 1;
                    f_act = calculate_f(cl_set, teams[i][j], const_list, k, lambda, &current_dev, &current_inf);
                    n_iters++;
                    do
                    {
                        f_ini = f_act;

                        // Neighborhood creation
                        for (int m = 0; m < teams[i][j].size(); m++)
                        {
                            for (int o = 1; o < k; o++)
                            {
                                c.clear();
                                c.push_back(m);
                                c.push_back((c[m] + o) % k);
                                neighbours.push_back(c);
                            }
                        }

                        random_shuffle(neighbours.begin(), neighbours.end());

                        // Neighborhood exploration
                        for (int m = 0; m < neighbours.size() && f_act >= f_ini; m++)
                        {
                            c = teams[i][j];
                            c[neighbours[m][0]] = neighbours[m][1];
                            f_act = calculate_f(cl_set, c, const_list, k, lambda, &current_dev, &current_inf);
                            bl_counter++;
                            n_iters++;
                        }

                        if (f_act < f_ini)
                            teams[i][j] = c;

                    } while (f_act < f_ini && bl_counter < max_iters_bl);
                }
            }
            else // Mutation
            {
                for (int j = 0; j < n_pilot_per_team; j++)
                {
                    n_mutations += teams[i].size() * teams[i][j].size();
                    if (n_mutations >= prob_mutation)
                    {
                        n_mutations = n_mutations % prob_mutation;
                        ind = rand() % teams[i].size();
                        val = rand() % teams[i][j].size();
                        teams[i][ind][val] = (teams[i][ind][val] + rand() % (k - 1) + 1) % k;
                    }
                }
            }
        }

        // If start new season, pilots can change their team

        if (n_iters > (current_season * season_lenght))
            current_season++;

        if (current_season > last_season)
            n_team_changes = n_team_changes_silly_season;
        else
            n_team_changes = n_team_changes_mid_season;

        for (int i = 0; i < n_team_changes; i++)
        {
            ind1 = rand() % teams.size();
            val1 = rand() % teams[0].size();
            ind2 = rand() % teams.size();
            val2 = rand() % teams[0].size();
            c = teams[ind1][val1];
            teams[ind1][val1] = teams[ind2][val2];
            teams[ind2][val2] = c;
        }

        last_season = current_season;

        // Fix best values and f-matrix
        // Se puede intentar introducir en train para optimizar

        f_matrix.clear();
        for (int i = 0; i < n_teams; i++)
        {
            cf.clear();
            for (int j = 0; j < n_pilot_per_team; j++)
            {
                current_f = calculate_f(cl_set, teams[i][j], const_list, k, lambda, &current_dev, &current_inf);
                n_iters++;
                if (current_f < best_f)
                {
                    best_f = current_f;
                    best_inf = current_inf;
                    best_dev = current_dev;
                }
                cf.push_back(current_f);
            }
            f_matrix.push_back(cf);
        }

        // Fix teams hierarchy

        mean_f = 0;
        team_hierarchy.clear();
        for (int i = 0; i < n_teams; i++)
        {
            i_aux = 0;
            for (int j = 0; j < n_pilot_per_team; j++)
            {
                i_aux += f_matrix[i][j];
                mean_f += f_matrix[i][j];
            }
            team_hierarchy.push_back(i_aux / n_pilot_per_team);
        }
        mean_f = mean_f / (n_teams * n_pilot_per_team);

    } while (n_iters <= max_iters);

    // Show f matrix

    /*for(int i=0;i<n_teams;i++){
        for(int j=0;j<n_pilot_per_team;j++){
            cout << f_matrix[i][j] << " ";
        }
        cout << endl;
    }*/

    // Show hierarchy

    /*for(int i=0;i<team_hierarchy.size();i++) cout << team_hierarchy[i] << " ";
    cout << endl;

    cout << mean_f << endl;*/

    time2 = clock();

    // Output

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
    string string_aux;           // Auxiliary string
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

    silly_season(&cl_set, &cl_set_const, k, seed);

    return 0;
}