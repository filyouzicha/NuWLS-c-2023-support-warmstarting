#ifndef Alg_NUWEIGHTING_h
#define Alg_NUWEIGHTING_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h> //these two h files are for timing in linux
#include <unistd.h>

using namespace std;

#define mypop(stack) stack[--stack##_fill_pointer]
#define mypush(item, stack) stack[stack##_fill_pointer++] = item

const float MY_RAND_MAX_FLOAT = 10000000.0;
const int MY_RAND_MAX_INT = 10000000;
const float BASIC_SCALE = 0.0000001; // 1.0f/MY_RAND_MAX_FLOAT;
const int USING_NEIGHBOR_MODE = 3;   // 1. using 2. don't use 3. depends on ins

// Define a data structure for a literal.
struct lit
{
    int clause_num; // clause num, begin with 0
    int var_num;    // variable num, begin with 1
    bool sense;     // is 1 for true literals, 0 for false literals.
};

struct clauselit
{
    int var_num; // variable num, begin with 1
    bool sense;  // is 1 for true literals, 0 for false literals.
};

struct varlit
{
    int clause_num; // clause num, begin with 0
    bool sense;     // is 1 for true literals, 0 for false literals.
};

static struct tms start_time;
static double get_runtime()
{
    struct tms stop;
    times(&stop);
    return (double)(stop.tms_utime - start_time.tms_utime + stop.tms_stime - start_time.tms_stime) / sysconf(_SC_CLK_TCK);
}
static void start_timing()
{
    times(&start_time);
}

class NUWEIGHTING
{
public:
    /***********non-algorithmic information ****************/
    int problem_weighted;
    int partial; // 1 if the instance has hard clauses, and 0 otherwise.
    int pure_sat;

    int max_clause_length;
    int min_clause_length;

    // size of the instance
    int num_vars;    // var index from 1 to num_vars
    int num_clauses; // clause index from 0 to num_clauses-1
    int num_hclauses;
    int num_sclauses;

    // steps and time
    int tries;
    int max_flips;
    int max_non_improve_flip;
    int step;

    int print_time;
    int cutoff_time;
    int prioup_time;
    double opt_time;

    /**********end non-algorithmic information*****************/
    /* literal arrays */
    varlit **var_lit;       // var_lit[i][j] means the j'th literal of var i.
    int *var_lit_count;     // amount of literals of each var
    clauselit **clause_lit; // clause_lit[i][j] means the j'th literal of clause i.
    int *clause_lit_count;  // amount of literals in each clause

    /* Information about the variables. */
    double *score;
    int *time_stamp;
    int **var_neighbor;
    int *var_neighbor_count;
    int *neighbor_flag = NULL;
    int *temp_neighbor = NULL;
    bool if_using_neighbor = false;

    /* Information about the clauses */
    unsigned long long top_clause_weight;
    long long *org_clause_weight;
    long long total_soft_weight;
    double *clause_weight;
    int *sat_count;
    int *sat_var;
    int *best_soft_clause = NULL;

    // original unit clause stack
    // lit *unit_clause;
    int unit_clause_count;

    // unsat clauses stack
    int *hardunsat_stack;          // store the unsat clause number
    int *index_in_hardunsat_stack; // which position is a clause in the unsat_stack
    int hardunsat_stack_fill_pointer;

    int *softunsat_stack;          // store the unsat clause number
    int *index_in_softunsat_stack; // which position is a clause in the unsat_stack
    int softunsat_stack_fill_pointer;

    // variables in unsat clauses
    int *unsatvar_stack;
    int unsatvar_stack_fill_pointer;
    int *index_in_unsatvar_stack;
    int *unsat_app_count = NULL; // a varible appears in how many unsat clauses

    // good decreasing variables (dscore>0 and confchange=1)
    int *goodvar_stack;
    int goodvar_stack_fill_pointer;
    int *already_in_goodvar_stack;
    int *score_change_stack = NULL;
    int score_change_stack_fill_pointer;
    bool *if_score_change = NULL;

    /* Information about solution */
    vector<int> binary_array;
    int *cur_soln; // the current solution, with 1's for True variables, and 0's for False variables
    int *best_soln;
    int *local_opt_soln = NULL;
    int best_soln_feasible; // when find a feasible solution, this is marked as 1.
    int local_soln_feasible;
    int hard_unsat_nb;
    unsigned long long soft_unsat_weight;
    unsigned long long opt_unsat_weight;
    unsigned long long local_opt_unsat_weight;

    // clause weighting
    int *large_weight_clauses;
    int large_weight_clauses_count;
    int large_clause_count_threshold;

    int *soft_large_weight_clauses;
    int *already_in_soft_large_weight_stack;
    int soft_large_weight_clauses_count;
    int soft_large_clause_count_threshold;

    // tem data structure used in algorithm
    int *best_array = NULL;
    int best_count;
    int *temp_lit = NULL;

    // parameters used in algorithm
    string init_soln_file;
    float rwprob;
    float rdprob;
    float smooth_probability;
    int hd_count_threshold;
    double h_inc;
    double softclause_weight_threshold;

    // clause weight tuned
    int coe_soft_clause_weight;
    double *tuned_org_clause_weight;
    double coe_tuned_weight;
    float soft_smooth_probability;
    double s_inc;
    long long total_soft_length;
    int NUWEIGHTING_TIME_LIMIT;

    // function used in algorithm
    void build_neighbor_relation();
    void allocate_memory();
    bool verify_sol();
    bool verify_goodvarstack(int flipvar);
    void smooth_weights();
    void hard_smooth_weights();
    void soft_smooth_weights();
    void hard_increase_weights();
    void soft_increase_weights();
    void update_clause_weights();
    void unsat(int clause);
    void sat(int clause);
    void init(vector<int> &init_solution);
    void flip(int flipvar);
    void flip2(int flipvar);
    void update_goodvarstack1(int flipvar);
    void update_goodvarstack2(int flipvar);
    int pick_var();
    long long floorToPowerOfTen(double x);
    long long closestPowerOfTen(double num);

    void soft_increase_weights_partial();
    void soft_increase_weights_not_partial();
    int *soft_clause_num_index;

    NUWEIGHTING();
    void settings();
    void build_instance(char *filename);
    void build_instance(int numVars, int numClauses, unsigned long long topClauseweight, clauselit **nuweighting_clause, int *nuweighting_clause_lit_count, long long *nuweighting_clause_weight);
    // void local_search(vector<int> &init_solution);
    // void local_search_with_decimation(char *inputfile);
    void simple_print();
    void print_best_solution();
    void free_memory();
    bool parse_parameters(int argc, char **argv);
    void get_init_solution(string& filename, vector<int>& binary_array);
};

inline NUWEIGHTING::NUWEIGHTING() {}

inline bool NUWEIGHTING::parse_parameters(int argc, char **argv)
{
    int i = 0;
    for (i = 1; i < argc; i++)
    {
        if (0 == strcmp(argv[i], "-init_soln_file"))
        {
            i++;
            if (i >= argc)
                return false;
            init_soln_file = argv[i];
            cout << "c The file containing the initial assignment for NuWLS: " << init_soln_file << endl;
            get_init_solution(init_soln_file, binary_array);
        }
    }
    
    return true;
}

inline void NUWEIGHTING::get_init_solution(string& filename, vector<int>& binary_array)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "c Unable to open file: " << filename << endl;
        exit(1);
    }

    string line;
    if (getline(file, line)) {
        size_t pos = line.find('v');
        if (pos != string::npos && pos + 1 < line.length()) {
            string binary_string = line.substr(pos + 1);
            for (char bit : binary_string) {
                if (bit == '0')
                    binary_array.push_back(0);
                else if (bit == '1')
                    binary_array.push_back(1);
                else if(!std::isspace(bit)) {
                    cerr << "c Invalid binary digit encountered: " << bit << endl;
                    exit(1);
                }
            }
        } else {
            cerr << "c Invalid format in the file." << endl;
            exit(1);
        }
    } else {
        cerr << "c Empty file: " << filename << endl;
        exit(1);
    }
    if (binary_array.size() != num_vars)
    {
        cerr << "c The number of variables in the file " << filename << " is incorrect." << endl;
    }
    file.close();
    return;
}

inline long long NUWEIGHTING::closestPowerOfTen(double num)
{
    if (num <= 1)
        return 1;

    int n = ceil(log10(num));
    int x = round(num / pow(10, n - 1));

    if (x == 10)
    {
        x = 1;
        n += 1;
    }
    return pow(10, n - 1) * x;
}

inline long long NUWEIGHTING::floorToPowerOfTen(double x)
{
    if (x <= 0.0) // if x <= 0, then return 0.
    {
        return 0;
    }
    int exponent = (int)log10(x);
    double powerOfTen = pow(10, exponent);
    long long result = (long long)powerOfTen;
    if (x < result)
    {
        result /= 10;
    }
    return result;
}

inline void NUWEIGHTING::settings()
{
    local_soln_feasible = 1;
    NUWEIGHTING_TIME_LIMIT = 15;
    if (1 == problem_weighted)
    {
        cout << "c problem weighted = 1" << endl;
        max_flips = 10000000;
        max_non_improve_flip = 10000000;
        large_clause_count_threshold = 0;
        soft_large_clause_count_threshold = 0;

        if (num_hclauses > 0) // weighted partial
        {
            NUWEIGHTING_TIME_LIMIT = 20;
            h_inc = 5;
            s_inc = 6;
            // coe_soft_clause_weight = 1000;
            hd_count_threshold = 50;
            rdprob = 0.036;
            rwprob = 0.48;
            soft_smooth_probability = 2E-6;
            smooth_probability = 2E-5;
            softclause_weight_threshold = 50;
            coe_tuned_weight = 1.0 / (double(top_clause_weight - 1) / (double)(num_sclauses));
            for (int c = 0; c < num_clauses; c++)
            {
                if (org_clause_weight[c] != top_clause_weight)
                {
                    tuned_org_clause_weight[c] = (double)org_clause_weight[c] * coe_tuned_weight;
                }
            }
        }
        else // weighted not partial
        {
            softclause_weight_threshold = 10;
            s_inc = 3.0;
            NUWEIGHTING_TIME_LIMIT = 25;
            soft_smooth_probability = 1E-3;
            hd_count_threshold = 22;
            rdprob = 0.036;
            rwprob = 0.48;

            coe_soft_clause_weight = 1000;
            coe_tuned_weight = ((double)coe_soft_clause_weight) / ((double(top_clause_weight - 1) / (double)(num_sclauses)));
            // cout << "c coe_tuned_weight: " << coe_tuned_weight << endl;
            for (int c = 0; c < num_clauses; c++)
            {
                tuned_org_clause_weight[c] = org_clause_weight[c] * coe_tuned_weight;
            }
        }
    }
    else
    {
        cout << "c problem weighted = 0" << endl;
        max_flips = 10000000;
        max_non_improve_flip = 10000000;

        large_clause_count_threshold = 0;
        soft_large_clause_count_threshold = 0;

        h_inc = 1;
        s_inc = 1;

        if (num_hclauses > 0) // unweighted partial
        {
            hd_count_threshold = 50;
            coe_soft_clause_weight = 1;
            rdprob = 0.079;
            rwprob = 0.087;
            soft_smooth_probability = 1E-5;
            softclause_weight_threshold = 500;
            smooth_probability = 1E-4;
        }
        else // unweighted not partial
        {
            s_inc = 1;
            NUWEIGHTING_TIME_LIMIT = 15;
            hd_count_threshold = 94;
            coe_soft_clause_weight = 397;
            rdprob = 0.007;
            rwprob = 0.047;
            soft_smooth_probability = 0.002;
            softclause_weight_threshold = 550;
        }
    }
}

inline void NUWEIGHTING::allocate_memory()
{
    int malloc_var_length = num_vars + 5;
    int malloc_clause_length = num_clauses + 5;
    int i = 0;

    var_lit = new varlit *[malloc_var_length];
    var_lit_count = new int[malloc_var_length];
    // clause_lit = new lit *[malloc_clause_length];
    // clause_lit_count = new int[malloc_clause_length];
    score = new double[malloc_var_length];
    var_neighbor = new int *[malloc_var_length];
    for (i = 0; i < malloc_var_length; i++)
    {
        var_lit[i] = NULL;
        var_neighbor[i] = NULL;
        var_lit_count[i] = 0;
    }
    var_neighbor_count = new int[malloc_var_length];
    time_stamp = new int[malloc_var_length];
    neighbor_flag = new int[malloc_var_length];
    temp_neighbor = new int[malloc_var_length];

    // org_clause_weight = new long long[malloc_clause_length];
    clause_weight = new double[malloc_clause_length];
    sat_count = new int[malloc_clause_length];
    sat_var = new int[malloc_clause_length];
    // best_soft_clause = new int[malloc_clause_length];

    hardunsat_stack = new int[malloc_clause_length];
    index_in_hardunsat_stack = new int[malloc_clause_length];
    softunsat_stack = new int[malloc_clause_length];
    index_in_softunsat_stack = new int[malloc_clause_length];

    unsatvar_stack = new int[malloc_var_length];
    index_in_unsatvar_stack = new int[malloc_var_length];
    // unsat_app_count = new int[malloc_var_length];

    goodvar_stack = new int[malloc_var_length];
    already_in_goodvar_stack = new int[malloc_var_length];
    // score_change_stack = new int[malloc_var_length];
    // if_score_change = new bool[malloc_var_length];

    cur_soln = new int[malloc_var_length];
    best_soln = new int[malloc_var_length];
    // local_opt_soln = new int[malloc_var_length];

    large_weight_clauses = new int[malloc_clause_length];
    soft_large_weight_clauses = new int[malloc_clause_length];
    already_in_soft_large_weight_stack = new int[malloc_clause_length];
    soft_clause_num_index = new int[malloc_clause_length];

    // best_array = new int[malloc_var_length];
    // temp_lit = new int[malloc_var_length];

    tuned_org_clause_weight = new double[malloc_clause_length];
}

inline void NUWEIGHTING::free_memory()
{
    int i;
    for (i = 0; i < num_clauses; i++)
        delete[] clause_lit[i];

    for (i = 1; i <= num_vars; ++i)
    {
        delete[] var_lit[i];
        delete[] var_neighbor[i];
    }

    delete[] var_lit;
    delete[] var_lit_count;
    delete[] clause_lit;
    delete[] clause_lit_count;

    delete[] score;
    delete[] var_neighbor;
    delete[] var_neighbor_count;
    delete[] time_stamp;
    delete[] neighbor_flag;
    delete[] temp_neighbor;

    delete[] org_clause_weight;
    delete[] clause_weight;
    delete[] sat_count;
    delete[] sat_var;
    // delete[] best_soft_clause;

    delete[] hardunsat_stack;
    delete[] index_in_hardunsat_stack;
    delete[] softunsat_stack;
    delete[] index_in_softunsat_stack;

    delete[] unsatvar_stack;
    delete[] index_in_unsatvar_stack;
    // delete[] unsat_app_count;

    delete[] goodvar_stack;
    delete[] already_in_goodvar_stack;
    // delete[] if_score_change;
    // delete[] score_change_stack;

    // delete [] fix;
    delete[] cur_soln;
    delete[] best_soln;
    // delete[] local_opt_soln;

    delete[] large_weight_clauses;
    delete[] soft_large_weight_clauses;
    delete[] already_in_soft_large_weight_stack;

    // delete[] best_array;
    // delete[] temp_lit;
    delete[] soft_clause_num_index;
    delete[] tuned_org_clause_weight;
}

inline void NUWEIGHTING::build_neighbor_relation()
{
    int i, j, count;
    int v, c, n;
    int temp_neighbor_count;

    for (v = 1; v <= num_vars; ++v)
    {
        neighbor_flag[v] = 1;
        temp_neighbor_count = 0;

        for (i = 0; i < var_lit_count[v]; ++i)
        {
            c = var_lit[v][i].clause_num;
            for (j = 0; j < clause_lit_count[c]; ++j)
            {
                n = clause_lit[c][j].var_num;
                if (neighbor_flag[n] != 1)
                {
                    neighbor_flag[n] = 1;
                    temp_neighbor[temp_neighbor_count++] = n;
                }
            }
        }

        neighbor_flag[v] = 0;

        var_neighbor[v] = new int[temp_neighbor_count];
        var_neighbor_count[v] = temp_neighbor_count;

        count = 0;
        for (i = 0; i < temp_neighbor_count; i++)
        {
            var_neighbor[v][count++] = temp_neighbor[i];
            neighbor_flag[temp_neighbor[i]] = 0;
        }
    }
}

inline void NUWEIGHTING::build_instance(char *filename)
{
    cout << "org build instance function" << endl;
}

inline void NUWEIGHTING::build_instance(int numVars, int numClauses, unsigned long long topClauseweight, clauselit **nuweighting_clause, int *nuweighting_clause_lit_count, long long *nuweighting_clause_weight)
{
    istringstream iss;
    string line;
    char tempstr1[10];
    char tempstr2[10];

    /*** build problem data structures of the instance ***/
    start_timing();
    total_soft_length = 0;
    num_vars = numVars;
    num_clauses = numClauses;
    top_clause_weight = topClauseweight;
    clause_lit = nuweighting_clause;
    clause_lit_count = nuweighting_clause_lit_count;
    org_clause_weight = nuweighting_clause_weight;

    allocate_memory();
    int v, c;

    int cur_lit;
    partial = 0;
    num_hclauses = num_sclauses = 0;
    max_clause_length = 0;
    min_clause_length = 100000000;
    unit_clause_count = 0;
    // int *redunt_test = new int[num_vars + 1];
    // memset(redunt_test, 0, sizeof(int) * num_vars + 1);
    // Now, read the clauses, one at a time.
    bool clause_reduent = false;
    long long total_lit_count = 0;
    for (int i = 0; i < num_clauses; ++i)
    {
        for (int j = 0; j < clause_lit_count[i]; ++j)
        {
            // int temv = clause_lit[i][j].var_num;
            // int temsense = clause_lit[i][j].sense;
            var_lit_count[clause_lit[i][j].var_num]++;
        }

        if (org_clause_weight[i] != top_clause_weight)
        {
            total_soft_weight += org_clause_weight[i];
            // total_soft_length += clause_lit_count[i];
            soft_clause_num_index[num_sclauses++] = i;
        }
        else
        {
            num_hclauses++;
            partial = 1;
        }
        // if (clause_lit_count[i] == 1)
        //     unit_clause[unit_clause_count++] = clause_lit[i][0];
        total_lit_count += clause_lit_count[i];
    }

    double total_memory = 0;
    // creat var literal arrays
    total_lit_count = 0;
    for (v = 1; v <= num_vars; ++v)
    {
        var_lit[v] = new varlit[var_lit_count[v] + 1];
        total_lit_count += (var_lit_count[v] + 1);
        var_lit_count[v] = 0; // reset to 0, for build up the array
    }
    //cout << "c total_lit_count " << total_lit_count;
    //cout << endl;
    // scan all clauses to build up var literal arrays
    for (int i = 0; i < num_clauses; ++i)
    {
        for (int j = 0; j < clause_lit_count[i]; ++j)
        {
            v = clause_lit[i][j].var_num;
            var_lit[v][var_lit_count[v]].clause_num = i;
            // var_lit[v][var_lit_count[v]].var_num = v;
            var_lit[v][var_lit_count[v]].sense = clause_lit[i][j].sense;
            ++var_lit_count[v];
        }
    }

    cout << "c before build neighbor" << endl;

    for (v = 1; v <= num_vars; ++v)
        var_lit[v][var_lit_count[v]].clause_num = -1;

    cout << "c build instime is " << get_runtime() << endl;

    if (USING_NEIGHBOR_MODE == 2 || 1 == problem_weighted || (USING_NEIGHBOR_MODE == 3 && ((get_runtime() > 1.0 || num_clauses > 10000000) && 0 == problem_weighted)))
    {
        if_using_neighbor = false;
    }
    else
    {
        cout << "c using neighbor " << endl;
        if_using_neighbor = true;
        build_neighbor_relation();
    }

    best_soln_feasible = 0;
    opt_unsat_weight = total_soft_weight + 1;
}

inline void NUWEIGHTING::init(vector<int> &init_solution)
{
    soft_large_weight_clauses_count = 0;
    // Initialize clause information
    if (1 == problem_weighted)
    {
        if (num_hclauses > 0) // weighted partial
        {
            for (int c = 0; c < num_clauses; c++)
            {
                already_in_soft_large_weight_stack[c] = 0;

                if (org_clause_weight[c] == top_clause_weight)
                    clause_weight[c] = 1;
                else
                {
                    clause_weight[c] = 0;
                }
            }
        }
        else // weighted not partial
        {
            for (int c = 0; c < num_clauses; c++)
            {
                already_in_soft_large_weight_stack[c] = 0;
                clause_weight[c] = tuned_org_clause_weight[c];
                if (clause_weight[c] > s_inc && already_in_soft_large_weight_stack[c] == 0)
                {
                    already_in_soft_large_weight_stack[c] = 1;
                    soft_large_weight_clauses[soft_large_weight_clauses_count++] = c;
                }
            }
        }
    }
    else
    {
        if (num_hclauses > 0) // unweighted partial
        {
            for (int c = 0; c < num_clauses; c++)
            {
                clause_weight[c] = 1;
            }
        }
        else // unweighted not partial
        {
            for (int c = 0; c < num_clauses; c++)
            {
                already_in_soft_large_weight_stack[c] = 0;

                clause_weight[c] = coe_soft_clause_weight;
                if (clause_weight[c] > 1 && already_in_soft_large_weight_stack[c] == 0)
                {
                    already_in_soft_large_weight_stack[c] = 1;
                    soft_large_weight_clauses[soft_large_weight_clauses_count++] = c;
                }
            }
        }
    }

    // init solution
    if (init_solution.size() == 0)
    {
        for (int v = 1; v <= num_vars; v++)
        {
            cur_soln[v] = rand() % 2;
            time_stamp[v] = 0;
        }
    }
    else
    {
        for (int v = 1; v <= num_vars; v++)
        {
            cur_soln[v] = init_solution[v];
            time_stamp[v] = 0;
        }
    }
    local_soln_feasible = 0;
    // init stacks
    hard_unsat_nb = 0;
    soft_unsat_weight = 0;
    hardunsat_stack_fill_pointer = 0;
    softunsat_stack_fill_pointer = 0;
    unsatvar_stack_fill_pointer = 0;
    large_weight_clauses_count = 0;

    /* figure out sat_count, sat_var and init unsat_stack */
    for (int c = 0; c < num_clauses; ++c)
    {
        sat_count[c] = 0;
        for (int j = 0; j < clause_lit_count[c]; ++j)
        {
            if (cur_soln[clause_lit[c][j].var_num] == clause_lit[c][j].sense)
            {
                sat_count[c]++;
                sat_var[c] = clause_lit[c][j].var_num;
            }
        }
        if (sat_count[c] == 0)
        {
            unsat(c);
        }
    }

    /*figure out score*/
    for (int v = 1; v <= num_vars; v++)
    {
        score[v] = 0.0;
        for (int i = 0; i < var_lit_count[v]; ++i)
        {
            int c = var_lit[v][i].clause_num;
            if (sat_count[c] == 0)
                score[v] += clause_weight[c];
            else if (sat_count[c] == 1 && var_lit[v][i].sense == cur_soln[v])
                score[v] -= clause_weight[c];
        }
    }

    // init goodvars stack
    goodvar_stack_fill_pointer = 0;
    score_change_stack_fill_pointer = 0;
    for (int v = 1; v <= num_vars; v++)
    {
        if (score[v] > 0)
        {
            already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
            mypush(v, goodvar_stack);
        }
        else
            already_in_goodvar_stack[v] = -1;
    }
}

inline void NUWEIGHTING::smooth_weights()
{
    int i, clause, v;

    for (i = 0; i < large_weight_clauses_count; i++)
    {
        clause = large_weight_clauses[i];
        if (sat_count[clause] > 0)
        {
            clause_weight[clause] -= h_inc;

            if (clause_weight[clause] == 1)
            {
                large_weight_clauses[i] = large_weight_clauses[--large_weight_clauses_count];
                i--;
            }
            if (sat_count[clause] == 1)
            {
                v = sat_var[clause];
                score[v] += h_inc;
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }

    for (i = 0; i < soft_large_weight_clauses_count; i++)
    {
        clause = soft_large_weight_clauses[i];
        if (sat_count[clause] > 0)
        {
            clause_weight[clause]--;
            if (clause_weight[clause] == 1 && already_in_soft_large_weight_stack[clause] == 1)
            {
                already_in_soft_large_weight_stack[clause] = 0;
                soft_large_weight_clauses[i] = soft_large_weight_clauses[--soft_large_weight_clauses_count];
                i--;
            }
            if (sat_count[clause] == 1)
            {
                v = sat_var[clause];
                score[v]++;
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }
}
inline void NUWEIGHTING::hard_increase_weights()
{
    int i, c, v;
    for (i = 0; i < hardunsat_stack_fill_pointer; ++i)
    {
        c = hardunsat_stack[i];

        clause_weight[c] += h_inc;

        if (clause_weight[c] == (h_inc + 1))
            large_weight_clauses[large_weight_clauses_count++] = c;

        for (clauselit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
        {
            score[v] += h_inc;
            if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
            {
                already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                mypush(v, goodvar_stack);
            }
        }
    }
    return;
}

inline void NUWEIGHTING::soft_increase_weights()
{
    int i, c, v;

    if (1 == problem_weighted)
    {
        for (i = 0; i < softunsat_stack_fill_pointer; ++i)
        {
            c = softunsat_stack[i];
            if (clause_weight[c] >= tuned_org_clause_weight[c] + softclause_weight_threshold)
                continue;
            else
                clause_weight[c] += s_inc;

            if (clause_weight[c] > s_inc && already_in_soft_large_weight_stack[c] == 0)
            {
                already_in_soft_large_weight_stack[c] = 1;
                soft_large_weight_clauses[soft_large_weight_clauses_count++] = c;
            }
            for (clauselit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
            {
                score[v] += s_inc;
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }
    else
    {
        for (i = 0; i < softunsat_stack_fill_pointer; ++i)
        {
            c = softunsat_stack[i];
            if (clause_weight[c] >= coe_soft_clause_weight + softclause_weight_threshold)
                continue;
            else
                clause_weight[c] += s_inc;

            if (clause_weight[c] > s_inc && already_in_soft_large_weight_stack[c] == 0)
            {
                already_in_soft_large_weight_stack[c] = 1;
                soft_large_weight_clauses[soft_large_weight_clauses_count++] = c;
            }
            for (clauselit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
            {
                score[v] += s_inc;
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }

    return;
}

inline void NUWEIGHTING::hard_smooth_weights()
{
    int i, clause, v;
    for (i = 0; i < large_weight_clauses_count; i++)
    {
        clause = large_weight_clauses[i];
        if (sat_count[clause] > 0)
        {
            clause_weight[clause] -= h_inc;

            if (clause_weight[clause] == 1)
            {
                large_weight_clauses[i] = large_weight_clauses[--large_weight_clauses_count];
                i--;
            }
            if (sat_count[clause] == 1)
            {
                v = sat_var[clause];
                score[v] += h_inc;
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }
    return;
}

inline void NUWEIGHTING::soft_smooth_weights()
{
    int i, clause, v;

    for (i = 0; i < soft_large_weight_clauses_count; i++)
    {
        clause = soft_large_weight_clauses[i];
        if (sat_count[clause] > 0)
        {
            clause_weight[clause] -= s_inc;
            if (clause_weight[clause] <= s_inc && already_in_soft_large_weight_stack[clause] == 1)
            {
                already_in_soft_large_weight_stack[clause] = 0;
                soft_large_weight_clauses[i] = soft_large_weight_clauses[--soft_large_weight_clauses_count];
                i--;
            }
            if (sat_count[clause] == 1)
            {
                v = sat_var[clause];
                score[v] += s_inc;
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }
    return;
}

inline void NUWEIGHTING::soft_increase_weights_partial()
{
    int i, c, v;

    if (1 == problem_weighted)
    {
        for (i = 0; i < num_sclauses; ++i)
        {
            c = soft_clause_num_index[i];
            clause_weight[c] += tuned_org_clause_weight[c];
            if (sat_count[c] <= 0) // unsat
            {
                for (clauselit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
                {
                    score[v] += tuned_org_clause_weight[c];
                    if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                    {
                        already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                        mypush(v, goodvar_stack);
                    }
                }
            }
            else if (sat_count[c] < 2) // sat
            {
                for (clauselit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
                {
                    if (p->sense == cur_soln[v])
                    {
                        score[v] -= tuned_org_clause_weight[c];
                        if (score[v] <= 0 && -1 != already_in_goodvar_stack[v])
                        {
                            int index = already_in_goodvar_stack[v];
                            int last_v = mypop(goodvar_stack);
                            goodvar_stack[index] = last_v;
                            already_in_goodvar_stack[last_v] = index;
                            already_in_goodvar_stack[v] = -1;
                        }
                    }
                }
            }
        }
    }
    else
    {
        for (i = 0; i < num_sclauses; ++i)
        {
            c = soft_clause_num_index[i];
            clause_weight[c] += s_inc;

            if (sat_count[c] <= 0) // unsat
            {
                for (clauselit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
                {
                    score[v] += s_inc;
                    if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                    {
                        already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                        mypush(v, goodvar_stack);
                    }
                }
            }
            else if (sat_count[c] < 2) // sat
            {
                for (clauselit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
                {
                    if (p->sense == cur_soln[v])
                    {
                        score[v] -= s_inc;
                        if (score[v] <= 0 && -1 != already_in_goodvar_stack[v])
                        {
                            int index = already_in_goodvar_stack[v];
                            int last_v = mypop(goodvar_stack);
                            goodvar_stack[index] = last_v;
                            already_in_goodvar_stack[last_v] = index;
                            already_in_goodvar_stack[v] = -1;
                        }
                    }
                }
            }
        }
    }
    return;
}

inline void NUWEIGHTING::soft_increase_weights_not_partial()
{
    int i, c, v;

    if (1 == problem_weighted)
    {
        for (i = 0; i < softunsat_stack_fill_pointer; ++i)
        {
            c = softunsat_stack[i];
            if (clause_weight[c] >= tuned_org_clause_weight[c] + softclause_weight_threshold)
                continue;
            else
                clause_weight[c] += s_inc;

            if (clause_weight[c] > s_inc && already_in_soft_large_weight_stack[c] == 0)
            {
                already_in_soft_large_weight_stack[c] = 1;
                soft_large_weight_clauses[soft_large_weight_clauses_count++] = c;
            }
            for (clauselit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
            {
                score[v] += s_inc;
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }
    else
    {
        for (i = 0; i < softunsat_stack_fill_pointer; ++i)
        {
            c = softunsat_stack[i];
            if (clause_weight[c] >= coe_soft_clause_weight + softclause_weight_threshold)
                continue;
            else
                clause_weight[c] += s_inc;

            if (clause_weight[c] > s_inc && already_in_soft_large_weight_stack[c] == 0)
            {
                already_in_soft_large_weight_stack[c] = 1;
                soft_large_weight_clauses[soft_large_weight_clauses_count++] = c;
            }
            for (clauselit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
            {
                score[v] += s_inc;
                if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
        }
    }
    return;
}

inline void NUWEIGHTING::update_clause_weights()
{
    if (num_hclauses > 0)
    {
        // update hard clause weight
        hard_increase_weights();
        if (0 == hard_unsat_nb)
        {
            soft_increase_weights_partial();
        }
    }
    else
    {
        if (((rand() % MY_RAND_MAX_INT) * BASIC_SCALE) < soft_smooth_probability && soft_large_weight_clauses_count > soft_large_clause_count_threshold)
        {
            soft_smooth_weights();
        }
        else
        {
            soft_increase_weights_not_partial();
        }
    }
}

inline int NUWEIGHTING::pick_var()
{
    int i, v;
    int best_var;

    if (goodvar_stack_fill_pointer > 0)
    {
        if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < rdprob)
            return goodvar_stack[rand() % goodvar_stack_fill_pointer];

        if (goodvar_stack_fill_pointer < hd_count_threshold)
        {
            best_var = goodvar_stack[0];
            for (i = 1; i < goodvar_stack_fill_pointer; ++i)
            {
                v = goodvar_stack[i];
                if (score[v] > score[best_var])
                    best_var = v;
                else if (score[v] == score[best_var])
                {
                    if (time_stamp[v] < time_stamp[best_var])
                        best_var = v;
                }
            }
            return best_var;
        }
        else
        {
            best_var = goodvar_stack[rand() % goodvar_stack_fill_pointer];
            for (i = 1; i < hd_count_threshold; ++i)
            {
                v = goodvar_stack[rand() % goodvar_stack_fill_pointer];
                if (score[v] > score[best_var])
                    best_var = v;
                else if (score[v] == score[best_var])
                {
                    if (time_stamp[v] < time_stamp[best_var])
                        best_var = v;
                }
            }
            return best_var;
        }
    }

    update_clause_weights();

    int sel_c;
    clauselit *p;

    if (hardunsat_stack_fill_pointer > 0)
    {
        sel_c = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];
    }
    else
    {
        sel_c = softunsat_stack[rand() % softunsat_stack_fill_pointer];
    }
    if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < rwprob)
        return clause_lit[sel_c][rand() % clause_lit_count[sel_c]].var_num;

    best_var = clause_lit[sel_c][0].var_num;
    p = clause_lit[sel_c];
    for (p++; (v = p->var_num) != 0; p++)
    {
        if (score[v] > score[best_var])
            best_var = v;
        else if (score[v] == score[best_var])
        {
            if (time_stamp[v] < time_stamp[best_var])
                best_var = v;
        }
    }

    return best_var;
}

inline void NUWEIGHTING::update_goodvarstack1(int flipvar)
{
    int v;
    // remove the vars no longer goodvar in goodvar stack
    for (int index = goodvar_stack_fill_pointer - 1; index >= 0; index--)
    {
        v = goodvar_stack[index];
        if (score[v] <= 0)
        {
            goodvar_stack[index] = mypop(goodvar_stack);
            already_in_goodvar_stack[goodvar_stack[index]] = index;
            already_in_goodvar_stack[v] = -1;
        }
    }

    // add goodvar
    for (int i = 0; i < var_neighbor_count[flipvar]; ++i)
    {
        v = var_neighbor[flipvar][i];
        if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
        {
            already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
            mypush(v, goodvar_stack);
        }
    }
}
/*
inline void NUWEIGHTING::update_goodvarstack2(int flipvar)
{
    int v;
    // remove the vars no longer goodvar in goodvar stack
    for (int index = goodvar_stack_fill_pointer - 1; index >= 0; index--)
    {
        v = goodvar_stack[index];
        if (score[v] <= 0)
        {
            goodvar_stack[index] = mypop(goodvar_stack);
            already_in_goodvar_stack[goodvar_stack[index]] = index;
            already_in_goodvar_stack[v] = -1;
        }
    }

    lit *clause_c;
    int c;

    for (int i = 0; i < score_change_stack_fill_pointer; i++)
    {
        int v = score_change_stack[i];
        // cout << "v is " << v << " score change time is " << score_change_time[v] << endl;
        if (!if_score_change[v])
            continue;

        if_score_change[v] = false;

        if (score[v] > 0 && already_in_goodvar_stack[v] == -1)
        {
            already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
            mypush(v, goodvar_stack);
        }
    }

    if (already_in_goodvar_stack[flipvar] != 0 && score[flipvar] > 0)
    {
        if (goodvar_stack[already_in_goodvar_stack[flipvar]] == flipvar)
        {
            int tem_v = mypop(goodvar_stack);
            goodvar_stack[already_in_goodvar_stack[flipvar]] = tem_v;
            already_in_goodvar_stack[tem_v] = already_in_goodvar_stack[flipvar];
            already_in_goodvar_stack[flipvar] = -1;
        }
    }
}
*/
inline void NUWEIGHTING::flip(int flipvar)
{
    int i, v, c, clen = 0;
    int index;
    clauselit *clause_c;
    score_change_stack_fill_pointer = 0;
    double org_flipvar_score = score[flipvar];
    // cout << "c org_flipvar_score: " <<org_flipvar_score << endl;
    cur_soln[flipvar] = 1 - cur_soln[flipvar];

    for (i = 0; i < var_lit_count[flipvar]; ++i)
    {
        c = var_lit[flipvar][i].clause_num;
        clause_c = clause_lit[c];

        if (cur_soln[flipvar] == var_lit[flipvar][i].sense)
        {
            ++sat_count[c];
            if (sat_count[c] == 2) // sat_count from 1 to 2
            {
                v = sat_var[c];
                score[v] += clause_weight[c];
                if (score[v] > 0 && -1 == already_in_goodvar_stack[v])
                {
                    already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                    mypush(v, goodvar_stack);
                }
            }
            else if (sat_count[c] == 1) // sat_count from 0 to 1
            {
                sat_var[c] = flipvar; // record the only true lit's var
                for (clen = 0; clen < clause_lit_count[c]; clen++)
                {
                    v = clause_lit[c][clen].var_num;
                    score[v] -= clause_weight[c];
                    if (score[v] <= 0 && -1 != already_in_goodvar_stack[v])
                    {
                        int index = already_in_goodvar_stack[v];
                        int last_v = mypop(goodvar_stack);
                        goodvar_stack[index] = last_v;
                        already_in_goodvar_stack[last_v] = index;
                        already_in_goodvar_stack[v] = -1;
                    }
                }
                sat(c);
            }
        }
        else // cur_soln[flipvar] != cur_lit.sense
        {
            --sat_count[c];
            if (sat_count[c] == 1) // sat_count from 2 to 1
            {
                for (clen = 0; clen < clause_lit_count[c]; clen++)
                {
                    v = clause_lit[c][clen].var_num;
                    if (clause_lit[c][clen].sense == cur_soln[v])
                    {
                        score[v] -= clause_weight[c];
                        if (score[v] <= 0 && -1 != already_in_goodvar_stack[v])
                        {
                            int index = already_in_goodvar_stack[v];
                            int last_v = mypop(goodvar_stack);
                            goodvar_stack[index] = last_v;
                            already_in_goodvar_stack[last_v] = index;
                            already_in_goodvar_stack[v] = -1;
                        }
                        sat_var[c] = v;
                        break;
                    }
                }
            }
            else if (sat_count[c] == 0) // sat_count from 1 to 0
            {
                for (clen = 0; clen < clause_lit_count[c]; clen++)
                {
                    v = clause_lit[c][clen].var_num;
                    score[v] += clause_weight[c];
                    if (score[v] > 0 && -1 == already_in_goodvar_stack[v])
                    {
                        already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
                        mypush(v, goodvar_stack);
                    }
                }
                unsat(c);
            } // end else if
        }     // end else
    }

    // update information of flipvar
    score[flipvar] = -org_flipvar_score;
    if (score[flipvar] > 0 && already_in_goodvar_stack[flipvar] == -1)
    {
        already_in_goodvar_stack[flipvar] = goodvar_stack_fill_pointer;
        mypush(flipvar, goodvar_stack);
    }
    else if (score[flipvar] <= 0 && already_in_goodvar_stack[flipvar] != -1)
    {
        int index = already_in_goodvar_stack[flipvar];
        int last_v = mypop(goodvar_stack);
        goodvar_stack[index] = last_v;
        already_in_goodvar_stack[last_v] = index;
        already_in_goodvar_stack[flipvar] = -1;
    }
}

inline void NUWEIGHTING::flip2(int flipvar)
{
    int i, v, c;
    int index;
    clauselit *clause_c;

    score_change_stack_fill_pointer = 0;
    double org_flipvar_score = score[flipvar];
    cur_soln[flipvar] = 1 - cur_soln[flipvar];

    for (i = 0; i < var_lit_count[flipvar]; ++i)
    {
        c = var_lit[flipvar][i].clause_num;
        clause_c = clause_lit[c];

        if (cur_soln[flipvar] == var_lit[flipvar][i].sense)
        {
            ++sat_count[c];
            if (sat_count[c] == 2) // sat_count from 1 to 2
            {
                score[sat_var[c]] += clause_weight[c];
            }
            else if (sat_count[c] == 1) // sat_count from 0 to 1
            {
                sat_var[c] = flipvar; // record the only true lit's var
                for (clauselit *p = clause_c; (v = p->var_num) != 0; p++)
                {
                    score[v] -= clause_weight[c];
                }
                sat(c);
            }
        }
        else // cur_soln[flipvar] != cur_lit.sense
        {
            --sat_count[c];
            if (sat_count[c] == 1) // sat_count from 2 to 1
            {
                for (clauselit *p = clause_c; (v = p->var_num) != 0; p++)
                {
                    if (p->sense == cur_soln[v])
                    {
                        score[v] -= clause_weight[c];
                        sat_var[c] = v;
                        break;
                    }
                }
            }
            else if (sat_count[c] == 0) // sat_count from 1 to 0
            {
                for (clauselit *p = clause_c; (v = p->var_num) != 0; p++)
                {
                    score[v] += clause_weight[c];
                }
                unsat(c);
            } // end else if
        }     // end else
    }

    // update information of flipvar
    score[flipvar] = -org_flipvar_score;
    update_goodvarstack1(flipvar);
}

inline void NUWEIGHTING::print_best_solution()
{
    if (best_soln_feasible == 0)
        return;

    printf("v");
    for (int i = 1; i <= num_vars; i++)
    {
        printf(" ");
        if (best_soln[i] == 0)
            printf("-");
        printf("%d", i);
    }
    printf("\n");
}

inline bool NUWEIGHTING::verify_sol()
{
    int c, j, flag;
    long long verify_unsat_weight = 0;

    for (c = 0; c < num_clauses; ++c)
    {
        flag = 0;
        for (j = 0; j < clause_lit_count[c]; ++j)
        {
            if (cur_soln[clause_lit[c][j].var_num] == clause_lit[c][j].sense)
            {
                flag = 1;
                break;
            }
        }
        if (flag == 0)
        {
            if (org_clause_weight[c] == top_clause_weight) // verify hard clauses
            {
                // output the clause unsatisfied by the solution
                cout << "c Error: hard clause " << c << " is not satisfied" << endl;

                cout << "c ";
                for (j = 0; j < clause_lit_count[c]; ++j)
                {
                    if (clause_lit[c][j].sense == 0)
                        cout << "-";
                    cout << clause_lit[c][j].var_num << " ";
                }
                cout << endl;
                cout << "c ";
                for (j = 0; j < clause_lit_count[c]; ++j)
                    cout << cur_soln[clause_lit[c][j].var_num] << " ";
                cout << endl;
                return 0;
            }
            else
            {
                verify_unsat_weight += org_clause_weight[c];
            }
        }
    }

    if (verify_unsat_weight == opt_unsat_weight)
    {
        cout << "c yes " << verify_unsat_weight << endl;
    }
    else
    {
        cout << "c Error: find opt=" << opt_unsat_weight << ", but verified opt=" << verify_unsat_weight << endl;
    }
    return 0;
}

inline bool NUWEIGHTING::verify_goodvarstack(int flipvar)
{
    for (int i = 1; i <= num_vars; ++i)
    {
        if (i == flipvar)
            continue;
        if (score[i] > 0 && already_in_goodvar_stack[i] == -1)
        {
            cout << "wrong 1 :" << endl;
            cout << "var is " << i << endl;
        }
        else if (score[i] <= 0 && already_in_goodvar_stack[i] != -1)
        {
            cout << "wrong 2 :" << endl;
            cout << "var is " << i << endl;
        }
        /*if (if_score_change[i] != 0)
        {
            cout << "wrong 3 :" << endl;
            cout << "var is " << i << endl;
        }*/
    }
    if (score[flipvar] > 0 && already_in_goodvar_stack[flipvar] != -1)
    {
        cout << "wrong flipvar in good var " << flipvar << endl;
        cout << score[flipvar] << endl;
    }
    return 1;
}

inline void NUWEIGHTING::simple_print()
{
    if (best_soln_feasible == 1)
    {
        if (verify_sol() == 1)
            cout << opt_unsat_weight << '\t' << opt_time << endl;
        else
            cout << "solution is wrong " << endl;
    }
    else
        cout << -1 << '\t' << -1 << endl;
}

inline void NUWEIGHTING::unsat(int clause)
{
    if (org_clause_weight[clause] == top_clause_weight)
    {
        index_in_hardunsat_stack[clause] = hardunsat_stack_fill_pointer;
        mypush(clause, hardunsat_stack);
        hard_unsat_nb++;
    }
    else
    {
        index_in_softunsat_stack[clause] = softunsat_stack_fill_pointer;
        mypush(clause, softunsat_stack);
        soft_unsat_weight += org_clause_weight[clause];
    }
}

inline void NUWEIGHTING::sat(int clause)
{
    int index, last_unsat_clause;

    if (org_clause_weight[clause] == top_clause_weight)
    {
        last_unsat_clause = mypop(hardunsat_stack);
        index = index_in_hardunsat_stack[clause];
        hardunsat_stack[index] = last_unsat_clause;
        index_in_hardunsat_stack[last_unsat_clause] = index;

        hard_unsat_nb--;
    }
    else
    {
        last_unsat_clause = mypop(softunsat_stack);
        index = index_in_softunsat_stack[clause];
        softunsat_stack[index] = last_unsat_clause;
        index_in_softunsat_stack[last_unsat_clause] = index;

        soft_unsat_weight -= org_clause_weight[clause];
    }
}

#endif
