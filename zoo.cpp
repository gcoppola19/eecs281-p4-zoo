// Project identifier: 3E33912F8BAA7542FC4A1585D2DB6FE0312725B9

#include <vector>
#include <iostream>
#include <algorithm> // std::sort
#include <getopt.h>
#include <string> // needed for VS
#include <deque>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <numeric>

using namespace std;
struct PrimData
{
    double d;  // lowest distance to MST
    int32_t p; // index of vertex's parent in MST
    bool k;    // whether or not vertex already in MST
    size_t index;
    PrimData() : d{std::numeric_limits<double>::infinity()}, p{-1}, k{false}, index{0} {}
};
enum class TypeVert
{
    wall,
    wild,
    other
};
struct Cage
{
    int x;
    int y;
    int vertice;
    TypeVert type;
};
enum class Mode
{
    kNone = 0,
    kMST,
    kFASTTSP,
    kOPTTSP
};

struct Options
{
    Mode mode = Mode::kNone;
};

class Zoo
{
public:
    void get_options(int arc, char *argv[], Options &options);
    void read_input(Options &options);
    double runMST();
    vector<size_t> runFSTTSP();
    void runOPTTSP();
    void genPerms(size_t permLength);
    double calculateD(size_t index1, size_t index2);
    double calcD(size_t index1, size_t index2);
    void clear();
    double calcCost(size_t end);
    double getMST(size_t start);
    bool promising(size_t permLength);
    void calculateDistances();
    Zoo()
    {
        numVertices = 0;
        toPrint = false;
        currCost = 0;
        currLen = 0;
        bestLen = 0;
        bestCost = 0;
    };

private:
    size_t numVertices;
    vector<Cage *> cages;
    bool toPrint;
    double currCost;
    size_t currLen;
    size_t bestLen;
    double bestCost;
    vector<size_t> bestPath;
    vector<size_t> path;
    vector<vector<double>> distances;
};

// Read and process command line options.
void Zoo::get_options(int argc, char *argv[], Options &options)
{

    opterr = false;
    int option_index = 0;
    int choice;

    // use getopt to find command line options
    option longOpts[] = {{"mode", required_argument, nullptr, 'm'},
                         {"help", no_argument, nullptr, 'h'},
                         {nullptr, 0, nullptr, '\0'}};

    // if output has no arg
    while ((choice = getopt_long(argc, argv, "m:h", longOpts, &option_index)) != -1)
    {
        switch (choice)
        {
        case 'm':
        {
            string arg{optarg};
            if (arg == "MST")
            {
                options.mode = Mode::kMST;
                // toPrint = true;
                break;
            }
            // figure out what to do if given both
            else if (arg == "FASTTSP")
            {
                options.mode = Mode::kFASTTSP;
                toPrint = true;
                break;
            }
            else if (arg == "OPTTSP")
            {
                options.mode = Mode::kOPTTSP;
                break;
            }
            else
            {

                cerr << "Error: Invalid mode argument";
                exit(1);
            }
            break;
        }
        case 'h':
            std::cout << "This program reads in a number of vertices followed ,\n"
                      << "by vertice pairs and then depending on the mode, \n"
                      << " it outputs a solution. IF MST, it finds the minimum spanning tree\n"
                      << "edges and returns the edges, if FSTTSP finds a fast solution, but not necessarily optimal\n"
                      << "for the travelingsalesperson problem, OPTTSP returns the most optimal solution.\n"
                      << "Usage: \'./zoo\n\t[--mode | -m] <MST or FSTTSP or OPTTSP\n"
                      << "\t[--help | -h]\n"
                      << std::endl;
            exit(0);
        default:
            cerr << "Error: invalid option"
                 << "\n";
            exit(1);
        }
    }

    if (options.mode == Mode::kNone)
    {
        cerr << "Error: invalid mode option";
        exit(1);
    }
}

void Zoo::read_input(Options &options)
{
    string temp;
    getline(cin, temp);
    string x;
    string y;
    int num = stoi(temp);
    numVertices = size_t(num);
    int count = 0;
    while (getline(cin, temp))
    {
        std::istringstream iss(temp);
        iss >> x >> y;
        Cage *cage = new Cage;
        cage->x = stoi(x);
        cage->y = stoi(y);
        cage->vertice = count;
        if (cage->x < 0 && cage->y < 0)
        {
            cage->type = TypeVert::wild;
        }
        else if ((cage->x == 0 && cage->y <= 0) || (cage->y == 0 && cage->x <= 0))
        {
            cage->type = TypeVert::wall;
        }
        else
        {
            cage->type = TypeVert::other;
        }
        cages.push_back(cage);
        count++;
    }

    if (options.mode == Mode::kMST)
    {
        runMST();
    }
    else if (options.mode == Mode::kFASTTSP)
    {
        runFSTTSP();
    }
    else
    {
        runOPTTSP();
    }

    // to do read in vertices and locations
}
double Zoo::calculateD(size_t index1, size_t index2)
{

    // Check if either vertex is located at negative coordinates
    if ((cages[index1]->type == TypeVert::wild && cages[index2]->type == TypeVert::other) || (cages[index2]->type == TypeVert::wild && cages[index1]->type == TypeVert::other))
    {
        return std::numeric_limits<double>::infinity();
    }

    double dx = static_cast<double>(cages[index1]->x - cages[index2]->x);
    double dy = static_cast<double>(cages[index1]->y - cages[index2]->y);

    return dx * dx + dy * dy;
}

double Zoo::runMST()
{
     vector<PrimData> weiners(numVertices);
    weiners[0].d = 0;

    for (size_t i = 0; i < numVertices; ++i) {
        int v_min = -1;
        for (size_t j = 0; j < numVertices; ++j) {
            if (!weiners[j].k && (v_min == -1 || weiners[j].d < weiners[static_cast<size_t>(v_min)].d)) {
                v_min = static_cast<int>(j);
            }
        }

        if (v_min == -1) 
            break; // All vertices processed

        weiners[static_cast<size_t>(v_min)].k = true;

        for (size_t n = 0; n < numVertices; ++n) {
            if (n != static_cast<size_t>(v_min)) {
                double weight = calculateD(static_cast<size_t>(v_min), n);

                if (!weiners[n].k && weight < weiners[n].d) {
                    weiners[n].d = weight;
                    weiners[n].p = v_min;
                }
            }
        }
    }

    double totalWeight = 0;

    if (totalWeight == numeric_limits<double>::infinity())
    {
        cerr << "Cannot construct MST\n";
        exit(1);
    }
    for (size_t i = 0; i < numVertices; ++i)
    {
        totalWeight += sqrt(weiners[i].d);
    }
    cout << totalWeight << "\n";
    for (size_t i = 1; i < numVertices; ++i)
    {
        if (static_cast<int>(i) < weiners[i].p)
        {
            cout << i << " " << weiners[i].p << "\n";
        }
        else
        {
            cout << weiners[i].p << " " << i << "\n";
        }
    }
    return totalWeight;
}

double Zoo::calcD(size_t index1, size_t index2)
{
    int x1 = cages[index1]->x;
    int y1 = cages[index1]->y;
    int x2 = cages[index2]->x;
    int y2 = cages[index2]->y;

    double dx = static_cast<double>(x1 - x2);
    double dy = static_cast<double>(y1 - y2);

    return sqrt(dx * dx + dy * dy);
}

vector<size_t> Zoo::runFSTTSP()
{
    vector<size_t> current_path;
    current_path.reserve(numVertices + 1); 

    current_path.push_back(0); // Start with the first vertex
    
    current_path.push_back(1);
    current_path.push_back(0);

    double totalDist = 0;

    for (size_t i = 1; i < numVertices; ++i)
    {
        double minIncrease = numeric_limits<double>::infinity();
            int insertIndex = -1;
        if (find(current_path.begin(), current_path.end(), i) == current_path.end())
        {
            for (size_t j = 1; j < current_path.size(); ++j)
            {
                double distIncrease = calcD(static_cast<size_t>(current_path[j - 1]), i) +
                                      calcD(i, static_cast<size_t>(current_path[j])) -
                                      calcD(static_cast<size_t>(current_path[j - 1]), static_cast<size_t>(current_path[j]));

                if (distIncrease < minIncrease)
                {
                    minIncrease = distIncrease;
                    insertIndex = static_cast<int>(j);
                }
            }
            if (insertIndex != -1)
            {
                current_path.insert(current_path.begin() + insertIndex, static_cast<vector<size_t>::value_type>(i));
            }
        }
    }
    // Calculate total distance of the path
    for (size_t i = 0; i < current_path.size() - 1; ++i)
    {
        totalDist += calcD(static_cast<size_t>(current_path[i]), static_cast<size_t>(current_path[i + 1]));
    }
   totalDist += calcD(static_cast<size_t>(current_path.back()), 0);

    if (toPrint)
    {
        cout << totalDist << "\n";
        for (size_t i = 0; i < current_path.size() - 1; ++i)
        {
            cout << current_path[i] << " ";
        }
        cout << "\n";
    }
    bestCost = totalDist;
    return current_path;
}


double Zoo::getMST(size_t start)
{
if(path.size()-start <=2){
        return distances[path[path.size()-1]][path[start]];
    }
    vector<PrimData> weiners(numVertices);

    // set d of starting vertex to 0
    for(size_t i = start; i < numVertices; ++i)
    {
        weiners[i].d = numeric_limits<double>::infinity();
        weiners[i].k = false;
        weiners[i].index = static_cast<size_t>(path[i]);
    }

    weiners[start].d = 0;

    
    for (size_t i = start; i < numVertices; ++i) {
        int v_min = -1;
        for (size_t j = start; j < numVertices; ++j) {
            if (!weiners[j].k && (v_min == -1 || weiners[j].d < weiners[static_cast<size_t>(v_min)].d)) {
                v_min = static_cast<int>(j);
            }
        }

        if (v_min == -1) 
            break; // All vertices processed

        weiners[static_cast<size_t>(v_min)].k = true;

        for (size_t n = start; n < numVertices; ++n) {
            if (n != static_cast<size_t>(v_min)) {
                double weight = calcD(weiners[static_cast<size_t>(v_min)].index, weiners[n].index);

                if (!weiners[n].k && weight < weiners[n].d) {
                    weiners[n].d = weight;
                    weiners[n].p = v_min;
                }
            }
        }
    }
    double totalWeight = 0;
    for (size_t i = start; i < numVertices; ++i)
    {
        totalWeight += weiners[i].d;
    }
    return totalWeight;
    
}

 
 void Zoo::genPerms(size_t permLength)
{
    if (permLength == path.size())
    {
        double addedcost1 = distances[path[path.size() - 1]][path[0]];
        currCost += addedcost1;
        if (currCost < bestCost)
        {
            bestCost = currCost;
            bestPath = path;
        }
        currCost -= addedcost1;
        return;
    } // if ..complete path

    if (!promising(permLength))
    {
        return;
    } // if ..not promising

    for (size_t i = permLength; i < path.size(); ++i)
    {
        swap(path[permLength], path[i]);
        double addedcost = distances[path[permLength - 1]][path[permLength]];
        currCost += addedcost;
        genPerms(permLength + 1);
        currCost -= addedcost;
        swap(path[permLength], path[i]);
    } // for ..unpermuted elements
} // genPerms()


bool Zoo::promising(size_t permLength)
{
    if (path.size() - permLength < 5)
    {
        return true;
    }
    // use lowerbound strategy
    //double mst = getMST(permLength);
    double lower = currCost;
    if (lower > bestCost)
    {
        return false;
    }
    // lower += calcCost(path, permLength);

    double arm1 = numeric_limits<double>::infinity();
    double arm2 = numeric_limits<double>::infinity();

    for (size_t i = permLength; i < bestLen; ++i)
    {
        double one = distances[path[permLength - 1]][path[i]];
        double two = distances[path[0]][path[i]];
        if (one < arm1)
        {
            arm1 = one;
        }
        if (two < arm2)
        {
            arm2 = two;
        }
    }
    lower += arm1;
    lower+=arm2;
    if (lower > bestCost)
    {
        return false;
    }
    lower+=getMST(permLength);
    
    //lower += currCost;
    /*bool promise = lower < bestCost;
    for (size_t i = 0; i < path.size(); ++i)
        cerr << setw(2) << path[i] << ' ';
    cerr << setw(4) << permLength << setw(10) << currCost;
    cerr << setw(10) << arm1 << setw(10) << arm2;
    cerr << setw(10) << mst << setw(10) << lower << "  " << promise << '\n';*/
    /* if( lower < bestCost)
     {
         bestCost = lower;
     }*/
    return lower < bestCost;
}

void Zoo::calculateDistances() {
        distances.resize(numVertices, vector<double>(numVertices, 0.0));

        for (size_t i = 0; i < numVertices; ++i) {
            for (size_t j = i + 1; j < numVertices; ++j) {
                distances[i][j] = calcD(i, j);
                distances[j][i] = distances[i][j]; // Store symmetrically
            }
        }
    }

void Zoo::runOPTTSP()
{
    calculateDistances();
    bestPath = runFSTTSP();
    bestPath.pop_back();
    std::reverse(bestPath.begin() + 1, bestPath.end());
    //std::iota(bestPath.begin(), bestPath.end(), 0);
   // bestLen = bestPath.size();
    path = bestPath;

    genPerms(1);

    cout << bestCost << "\n";
    for (size_t i = 0; i < bestPath.size() - 1; ++i)
    {
        cout << bestPath[i] << " ";
    }
    cout << bestPath[bestPath.size() - 1] << "\n";
}


void Zoo::clear()
{
    for (auto &element : cages)
    {
        delete element;
        element = nullptr;
    }
    cages.clear();
}

int main(int argc, char **argv)
{

    ios_base::sync_with_stdio(false);
    cout << std::setprecision(2); // Always show 2 decimal places
    cout << std::fixed;           // Disable scientific notation for large numbers
   // cerr << std::boolalpha;
    Zoo zoo;
    Options options;

    // Read and process the command line options.
    zoo.get_options(argc, argv, options);
    zoo.read_input(options);

    zoo.clear();
    return 0;
}


/*vector<size_t> Zoo::runFSTTSP()
{
    double totalDist = 0;
    vector<size_t> current_path;
    current_path.reserve(numVertices + 1);

    current_path.push_back(0); // Start with the first vertex
    current_path.push_back(1);
    current_path.push_back(0);

   totalDist += 2 * calcD(0, 1);

    for (size_t i = 2; i < numVertices; ++i)
    {
        double minIncrease = numeric_limits<double>::infinity();
        int insertIndex = -1;
        for (size_t j = 0; j < current_path.size() - 1; ++j)
        {
            double distIncrease = calcD(static_cast<size_t>(current_path[j]), i) +
                                  calcD(i, static_cast<size_t>(current_path[j + 1])) -
                                  calcD(static_cast<size_t>(current_path[j]), static_cast<size_t>(current_path[j + 1]));

            if (distIncrease < minIncrease)
            {
                minIncrease = distIncrease;
                insertIndex = static_cast<int>(j + 1);
            }
        }
        if (insertIndex != -1)
        {
            totalDist += minIncrease;
            current_path.insert(current_path.begin() + insertIndex, static_cast<vector<size_t>::value_type>(i));
        }
    }
    // Calculate total distance of the path

        cout << totalDist << "\n";
        for (size_t i = 0; i < current_path.size() - 1; ++i)
        {
            cout << current_path[i] << " ";
        }
        cout << "\n";
    
    bestCost = totalDist;
    return current_path;
}

*/
 
