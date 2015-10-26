/**
 * @file  JSD.cpp
 * @brief calculate distance of two k-mer distributions using Jensen-Shannon divergence (JSD)
 *
 * Copyright (c) 2015 Oak Ridge National Lab. All rights reserved.<br />
 *
 * Author: Tae-Hyuk (Ted) Ahn 
 * Created: 01/23/2014
 * Modified: V0.0.2 use unordered_map instead of map for speeding up
 */


#include <stdio.h>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <map>
#include <unordered_map>


// acceptable
using namespace std;


// Help usage
void usage();

// If error, then print error and exit program 
void exitWithError(const std::string & error);

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss 
std::string currentDateTime();

// declare global variables
string program_name;
string version = "0.0.3";


// =============================================================================
// Help usage
// =============================================================================
void usage()
{
    cout << std::endl
         << "  [Usage]" << std::endl
         << "    " << program_name 
         << " -f <first ffp profile> -s <second ffp profile> -o <output file>" << endl
         << endl
         << "  [Inputs]" << endl
         << "    ffp profile file should have \"Sequence Frequency\" format" << endl
         << endl
         << "  [Outputs]" << endl
         << "    report distance" << endl
         << endl
         << "  [Options]" << endl
         << "    -t <int> : number of threads (default: 1)" << endl
         << "    -h       : print help" << endl
         << "    -v       : print version" << endl
         << endl;
}


// ============================================================================
// Initialize arguments
// ============================================================================
void initializeArguments(int argc, char **argv,  // (in) argc, argv
                         int & number_threads,   // (out) number of multi-threads
                         string & ffp1_file,     // (out) ffp1 input file
                         string & ffp2_file,     // (out) ffp2 input file
                         string & output_file)   // (out) output file
{
    // initialize variables
    number_threads = 1;

    // grab command line arguments and push
    std::vector<std::string> arguments_vector;
    while(argc--) {
        arguments_vector.push_back(*argv++);
    }

    // get argements
    for(int i = 1; i <= (int)arguments_vector.size()-1; i++)
    {
        if (arguments_vector[i] == "-f") {      // ffp1_file
            ffp1_file = arguments_vector[++i];
        }
        else if (arguments_vector[i] == "-s") { // ffp2_file
            ffp2_file = arguments_vector[++i];
        }
        else if (arguments_vector[i] == "-o") { // output
            output_file = arguments_vector[++i];
        }
        else if (arguments_vector[i] == "-t") { // multi-threads
            std::stringstream(arguments_vector[++i]) >> number_threads;
            if (number_threads < 1)
                exitWithError("*** Error: check -t option value.");
        }
        else if (arguments_vector[i] == "-v") { // version
            cout << endl << program_name << " V" << version << endl << endl;
            exit(0);
        }
        else if (arguments_vector[i] == "-h") { // help
            usage();
            exit(0);
        }
        else {// unknown option
            cerr << "*** Error: Unknown option " << arguments_vector[i] << endl << endl;
            usage();
            exit(1);
        }
    }
}


//=============================================================================
// Get current date/time, format is YYYY-MM-DD.HH:mm:ss 
//=============================================================================
std::string currentDateTime() 
{   
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "[%Y-%m-%d.%X]", &tstruct);

    // return 
    return buf;
}  


//=============================================================================
// Exit program when error happens 
//=============================================================================
void exitWithError(const string &error) 
{   
    cerr << endl << "  " << error << endl;
    exit(EXIT_FAILURE);
}   


//=============================================================================
// readFFP
//=============================================================================
void readFFP(const string &ffp_file, 
             unordered_map<string, double> &ffp_map,
             double &ffp_sum)
{
    ifstream ffp_input (ffp_file.c_str());    // ffp reads
    if (ffp_input) {   // check open file
        string line;
        while (getline(ffp_input, line))
        {
            std::istringstream ss(line);
            string sequence;
            double count;
            if (ss >> sequence >> count) {
                /*
                // if not exist
                if (ffp_map.find(sequence) == ffp_map.end() ) { 
                    ffp_map.insert ( std::pair<string, double>(sequence, count) );
                    ffp_sum += count;
                }   
                // else (exist)
                else {
                    cerr << endl << "*** Error: " << sequence << " exist more than one time!" << endl;
                    exit(EXIT_FAILURE);
                }   
                */ 
                // direct maching has a better performance
                ffp_map[sequence] = count;
                ffp_sum += count;
            } 
        }
    } 
    else {  // cannot open file
        usage();
        exit(0);
    }
    // close input_file                                                                                                     
    ffp_input.close();           
    
}

//=============================================================================
// runJSD
//=============================================================================
void runJSD(const int &number_threads, 
            const string &ffp1_file, 
            const string &ffp2_file, 
            const string &output_file)
{
    // set number of threads
    omp_set_num_threads(number_threads);

    // type def of unordered_map
    typedef unordered_map<string, double> string_double_map;

    // read first ffp profile
    string_double_map ffp1_map; // map
    double ffp1_sum;
    readFFP(ffp1_file, ffp1_map, ffp1_sum);

    // read second ffp profile
    string_double_map ffp2_map; // map
    double ffp2_sum;
    readFFP(ffp2_file, ffp2_map, ffp2_sum);

    // divided by sum to nomalize
    string_double_map::iterator it;

    for(it = ffp1_map.begin(); it != ffp1_map.end(); it++) {
        it->second = (it->second)/ffp1_sum;
    }

    for(it = ffp2_map.begin(); it != ffp2_map.end(); ++it) {
        it->second = (it->second)/ffp2_sum;
    }

    // Jensen-Shannon divergence
    string_double_map affp_map; // map
    for(it = ffp1_map.begin(); it != ffp1_map.end(); ++it) {
        if (ffp2_map.find(it->first) != ffp2_map.end() ) {  // key exist
            affp_map[it->first] = (it->second + ffp2_map[it->first])/2;
        }
        else {
            affp_map[it->first] = (it->second)/2;
        }
    }

    for(it = ffp2_map.begin(); it != ffp2_map.end(); ++it) {
        if (ffp1_map.find(it->first) == ffp1_map.end() ) {
            affp_map[it->first] = (it->second)/2;
        }
    }

    double JSD_value = 0.0;
    for(it = affp_map.begin(); it != affp_map.end(); ++it) {
        JSD_value = JSD_value - it->second * log(it->second)/log(2);
        if (ffp1_map.find(it->first) != ffp1_map.end() ) { // key exist
            JSD_value = JSD_value + 0.5 * ffp1_map[it->first] * log(ffp1_map[it->first])/log(2);
        }
        if (ffp2_map.find(it->first) != ffp2_map.end() ) { // key exist
            JSD_value = JSD_value + 0.5 * ffp2_map[it->first] * log(ffp2_map[it->first])/log(2);
        }
    }

    // print the JSD value to screen and output file
    ofstream output (output_file.c_str());
    cout << "JSD value=" << JSD_value << endl;
    output << "JSD value=" << JSD_value << endl;
}


//=============================================================================
// Main
//=============================================================================
int main(int argc, char** argv)
{
    // get program name
    string program_path = argv[0];
    program_name = program_path;
    size_t sep_pos= program_path.find_last_of('/');
    if (sep_pos != string::npos)
        program_name.assign(program_path.begin()+sep_pos+1,program_path.end());

    // initialize variables
    int number_threads;
    double start_time, finish_time, elapsed_time;
    string ffp1_file, ffp2_file;
    string output_file = "JSD_out.txt";

    // initialize arguments
    initializeArguments(argc, argv,  // (in)
        number_threads, ffp1_file, ffp2_file, output_file); // (out)

    // record start time (openMP)
    start_time = omp_get_wtime();

    // display work start and time record
    cout << endl
         << "********************************************************************************" << endl
         << currentDateTime() << " Beginning " << program_name << " V" << version <<  endl;

    // run JSD
    runJSD(number_threads, ffp1_file, ffp2_file, output_file);

    // finish time
    finish_time = omp_get_wtime();
    elapsed_time = double(finish_time - start_time);

    // display elapsed time
    cout << currentDateTime() << " Ending " << program_name << endl
        << "Total Elapsed Time =  " << elapsed_time << " [seconds]" << endl
        << "********************************************************************************"
        << endl << endl;
}

