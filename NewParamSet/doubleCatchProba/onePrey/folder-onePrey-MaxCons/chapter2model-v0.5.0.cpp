#include <iostream>  // input and output from the console
#include <string>    // manipulate char strings
#include <fstream>   // handle files
#include <time.h>    // get the time for random number generator
#include <stdlib.h>  // some random generator tools
#include <algorithm> // random shuffle
#include <vector>    // use vector objects (useful for short lists bc no need for data management)
#include <cmath>     // ceil
#include <random>    // random generation functions
#include <chrono>    // other needed tools for random generetion

using namespace std; // not to have to write std:: in front of every call

/* ---------------------------- Global variables ---------------------------- */

/* give a name to the simulation */
string simulationName;

/* landscape variables */
int worldSize;                // side of the squared lanscape in number of cells [0;+inf[
int resourceTypesNb;          // Number of resources available
vector<string> resourceTypes; // res1..n
vector<int> maxResources;     // maximum resources per cell

/* prey variables */
int preyTypesNb;                   // Number of preys modelled
vector<int> preyInitialDensities;  //
vector<string> preyTypes;          // prey1..n
vector<float> preyMaxMove;         // NOT GOOD should rather be defined in % of landscape size
vector<int> preyMaxConsume;        // in resource units
vector<int> preyMaintenanceCost;   //
vector<int> preyExpectedOffspring; //
vector<int> preyReproCost;         //
// add time of introduction?
// frequency of reproduction
// frequency of survival trial

/* predator variables */
int predatorTypesNb;               // Number of predators modelled
vector<string> predatorTypes;      // res1..n
vector<int> predInitialDensities;  //
vector<float> predMaxMove;         // NOT GOOD should rather be defined in % of landscape size CAREFUL the product has to stay int
vector<int> predMaxConsume;        // in equivalent resources
vector<int> predMaintenanceCost;   //
vector<int> predExpectedOffspring; //
vector<int> predReproCost;         //
vector<int> predIntro;             // time of introduction of the predator
// vector<float> predAsymm;           // ratio of preys' catch conversion into resources
vector<float> predCatchProba;      // pred catching probability
vector<int> predConvRate;         // number of resources a pred gets from 
vector<bool> predOportunistic;     // t/1 or f/0, is the predator oportunistic? (ranked prey types by conversion rate)
vector<bool> predSpecific;         // t/1 or f/0, is the predator specific? (hunts prey 1 in priority regardless of conversion rates)

/* time variables */
int timeMaxi; // simulation time
int freqRepr; // frequency of reproduction trial
int freqSurv; // frequency of survival trial
int freqResu; // frequency of save of the results
int freqSnap; // frequency of snapshot of the landscape
int freqRfll; // frequency of resources replenishment

/* seed for random number generator */
unsigned int randomSeed;

/* structures for the ibm */
int memberMatchingListsSize; // total number of types

/* pointers to members' matching lists (arrays) + size variables */
string *memberTypes;   //
int *typeTags;         //
int *indexInLandscape; // column index of each type in the landscape table

/* pointer to diets' table: will indicate who eats who and who does not */
int **dietsTable;

/* ---------------------------- Global functions ---------------------------- */

int sumColumn(int **table, int maxRow, int columnIndex) // sum of row values in a column
{
    int sum = 0;
    int r = 0;
    while (r < maxRow)
    {
        sum += table[r][columnIndex];
        r++;
    }

    return sum;
}

double randomNumberGenerator0to1(double min, double max) // generates a random number between min and max
{
    double res;
    res = min + (double)(rand()) / (RAND_MAX / (max - min));

    return res;
}

double randomNumberGeneratorAny(double min, double max) // generates a random number between min and max
{
    double res;
    res = min + (double)rand() * (max - min + 1) / (RAND_MAX - 1);

    return res;
}

int randomSampleFromPoissonDist(int Mean)
{
    int res;

    default_random_engine generator(randomNumberGeneratorAny(0, 99999999)); // random sample

    poisson_distribution<int> distribution(Mean); // generate a Poisson distribution

    res = distribution(generator); // random sample from the distribution

    return res;
}

vector<int> shuffleOrder(int populationSize)
{

    vector<int> popVector; // initialise a population vector

    for (int i = 0; i < populationSize; ++i) // create a vector of indexes of size of current population
        popVector.push_back(i);

    random_shuffle(popVector.begin(), popVector.end()); // shuffle indexes

    /* debug : OK
    cout << "popVector contains:";
    for (std::vector<int>::iterator it = popVector.begin(); it != popVector.end(); ++it)
        cout << ' ' << *it;

    cout << '\n';
    */

    return popVector;
}

void assignTagsIndexes() // matches names, tags and column index in landscapeTablePtr table.
{

    /* allocate memory */
    memberTypes = new string[memberMatchingListsSize];
    typeTags = new int[memberMatchingListsSize];
    indexInLandscape = new int[memberMatchingListsSize];

    /* assign values */

    /* choose a tag system: I have chose this one, allowing for 99 types of resources, preys an predators */
    int resourceTagStart = 101; // tag of the first resource on the list
    int preyTagStart = 201;     // tag of the first prey on the list
    int predatorTagStart = 301; // tag of the first predator on the list

    /* assign tags iteratively */
    int r = 0; // initialise row counter
    while (r < memberMatchingListsSize)
    {
        for (int res = 0; res < resourceTypesNb; res++)
        {
            memberTypes[r] = resourceTypes[res];
            typeTags[r] = resourceTagStart + res;
            indexInLandscape[r] = 3 + res; // before: cellCode - x - y

            /* debug : OK
            cout << "memberTypes[" << r << "] is " << memberTypes[r] << endl;
            cout << "typeTags[" << r << "] is " << typeTags[r] << endl;
            cout << "indexInLandscape[" << r << "] is " << indexInLandscape[r] << endl;
            */

            r++;
        }

        for (int prey = 0; prey < preyTypesNb; prey++)
        {
            memberTypes[r] = preyTypes[prey];
            typeTags[r] = preyTagStart + prey;
            indexInLandscape[r] = 3 + resourceTypesNb + prey; // before: x - y - resource densities

            /* debug : OK
            cout << "memberTypes[" << r << "] is " << memberTypes[r] << endl;
            cout << "typeTags[" << r << "] is " << typeTags[r] << endl;
            cout << "indexInLandscape[" << r << "] is " << indexInLandscape[r] << endl;
            */

            r++;
        }

        for (int pred = 0; pred < predatorTypesNb; pred++)
        {
            memberTypes[r] = predatorTypes[pred];
            typeTags[r] = predatorTagStart + pred;
            indexInLandscape[r] = 3 + resourceTypesNb + 2 * preyTypesNb + pred; // before: x - y - resource densities - prey densities - prey catches

            /* debug : OK
            cout << "memberTypes[" << r << "] is " << memberTypes[r] << endl;
            cout << "typeTags[" << r << "] is " << typeTags[r] << endl;
            cout << "indexInLandscape[" << r << "] is " << indexInLandscape[r] << endl;
            */

            r++;
        }
    }
}

int getMemberIndexFromTag(int TypeTag)
{

    int index = 0; // declare integer to be returned

    int row = 0;
    int rowMax = memberMatchingListsSize;
    while (row < rowMax)
    {
        if (TypeTag == typeTags[row]) // if tag corresponds
        {
            index = row;
            break; // exits the nested loops
        }

        row++;
    }

    return index; // returns the row index of the tag. It will be the same in all members matching lists.
}

void makeDietsTable(bool display) // this table allows for each member of the system to eat and be eaten by another
{

    int dietsTableSize = memberMatchingListsSize;
    // no need for headers when the memberMatchingListsIndex is known

    /* allocate memory. Make sure it is freed at the end of the main function! */
    dietsTable = new int *[dietsTableSize];

    for (int row = 0; row < dietsTableSize; row++)
    {
        dietsTable[row] = new int[dietsTableSize];
    }

    /* assign lines and columns headers : not needed here but useful for debug
    for (int row = 1; row < dietsTableSize; row++)
    {
        dietsTable[row][0] = typeTags[row - 1]; // -1 because we start row at 1 while we want typeTags[0]
    }

    for (int col = 1; col < dietsTableSize; col++)
    {
        dietsTable[0][col] = typeTags[col - 1];
    }
    */

    /* intialise all values to 0 */
    for (int row = 0; row < dietsTableSize; row++)
    {
        for (int col = 0; col < dietsTableSize; col++)
            dietsTable[row][col] = 0;
    }

    /* Set to 1 when one feeds on the other */
    dietsTable[0][2] = 1;                                             // prey1 feeds on resource1
    dietsTable[1][3] = 1;                                             // prey2 feeds on resource2
    dietsTable[2][4] = predConvRate[0];     // predator1 feeds on prey1 with a conversion rate of . HARD CODED not good
    dietsTable[3][4] = predConvRate[1];     // predator1 feeds on prey2 with a conversion rate  of

    if (display == true)
    {
        cout << "dietsTable" << endl;
        for (int row = 0; row < dietsTableSize; row++)
        {
            for (int col = 0; col < dietsTableSize; col++)
            {
                cout << dietsTable[row][col];
                if (col == (dietsTableSize - 1)) // last loop
                    cout << endl;
                else
                    cout << " ";
            }
        }
    }
}

vector<int> getDietLandscapeIndexes(int MembersMatchingListsIndex) // get the diet of a particular member of the food chain from its typeTag.
{

    vector<int> dietIndexes; // declare the vector to be returned

    // int typesInDiet = 0;
    // for (int i = 0; i < memberMatchingListsSize; i++)
    // {
    //     if (dietsTable[i][MembersMatchingListsIndex] > 0)
    //         typesInDiet++; // if animal is present in diet increase counter
    // }

    /* loop over lines to fill vector with diet's members column index in landscape table */
    for (int row = 0; row < memberMatchingListsSize; row++)
    {
        if (dietsTable[row][MembersMatchingListsIndex] > 0) // && dietIndexes.size() <= typesInDiet) // control not to add to the array more than its size
        {
            dietIndexes.push_back(indexInLandscape[row]); // the row number of this match is enough to know its index in the landscape table!
        }
    }

    /* debug : OK
    cout << memberTypes[MembersMatchingListsIndex] << "'s diet's member(s) column index(es) in landscape table is(are) : ";
    for (std::vector<int>::iterator it = dietIndexes.begin(); it != dietIndexes.end(); ++it)
        cout << ' ' << *it;

    cout << '\n';
    */

    return dietIndexes;
}

int getCellCode(string *xyCoordinates, int *CellCodes, int LandscapeSize, int x, int y)
{

    int cellCode = 0; // declare the integer to be returned

    /* turn x y coordinates into string */
    string XYcoord = to_string(x) + ";" + to_string(y);

    /* iterate through xyCoordinates to find match */
    int row = 0;
    int rowMax = LandscapeSize * LandscapeSize; // WARNING only works if squared
    while (row < rowMax)
    {
        if (XYcoord == xyCoordinates[row]) // if XY corresponds
        {
            cellCode = CellCodes[row]; // the row index is the same in all landscape matching lists
            break;
        }

        row++;
    }

    return cellCode;
}

int transmissiveBoundaries(int coordinate, int LandscapeSize)
{

    if (coordinate >= LandscapeSize)
    {
        coordinate -= LandscapeSize;
    }
    else if (coordinate < 0)
    {
        coordinate += LandscapeSize;
    }

    return coordinate;
}

void sortWhileKeepingIndex(vector<int> ConversionRateVec, vector<int> IndexesVec)
{
    int row = 0; // initiate row count

    while (row < (ConversionRateVec.size() - 1)) // until we reach the line before last
    {
        if (ConversionRateVec[row] >= ConversionRateVec[row + 1]) // if the focal element is greater than the next
        {
            row++; // leave as is and go to next line
        }
        else // if not switch positions and restart to the first line
        {
            int co1 = ConversionRateVec[row];
            int co2 = ConversionRateVec[row + 1];
            int in1 = IndexesVec[row];
            int in2 = IndexesVec[row + 1];

            ConversionRateVec[row] = co2;
            ConversionRateVec[row + 1] = co1;
            IndexesVec[row] = in2;
            IndexesVec[row + 1] = in1;

            row = 0; // restart row count
        }
    }
}

/* ----------------------------- Object classes ----------------------------- */

/*
create a class:
- all the variables that define it
- all the functions that manipulate these variables
- make sure that the values can be used by other classes
*/

class landscape
{
    /* list of landscape-specific variables */
private:                     // variables that should not be modified directly, nor accessed from the main function
    int Size;                // side of the squared lanscape in number of cells [0;+inf[
    vector<int> MaxResource; // max amount of resources on a cell

    int rowNb;               // row number
    int columnNb;            // column number
    int resColumnStart;      // indexes for table building convinience
    int preyColumnStart;     //
    int predatorColumnStart; //
    fstream resultsTable;    // file to write results in
    fstream snapshotTable;   // file to save all relevant landscape cell info

protected:                          // variables that should not be modified directly, nor accessed from the main function, but accessible to the other classes
    int landscapeMatchingListsSize; //

public:                      // can be used / called in the main function
    int **landscapeTablePtr; // pointer to the landscape table
    string *XYcoordinates;   // landscape matching lists
    int *cellCode;           //

    landscape(int size, vector<int> maxResourceVec) // constructor: function that creates objects of the class by assigning values to or initializing the variables
    {
        Size = size;
        MaxResource = maxResourceVec;

        /* column index where every group of info starts in the table*/
        resColumnStart = 3;                                      // before: cellCode - x - y
        preyColumnStart = resColumnStart + resourceTypesNb;      // before: x - y - resource densities
        predatorColumnStart = preyColumnStart + 2 * preyTypesNb; // before: x - y - resource densities - prey densities - prey catches

        /* table size */
        rowNb = Size * Size;
        columnNb = resColumnStart + resourceTypesNb + 2 * preyTypesNb + predatorTypesNb;

        /* create landscape matching lists */
        landscapeMatchingListsSize = rowNb;
        XYcoordinates = new string[landscapeMatchingListsSize];
        cellCode = new int[landscapeMatchingListsSize];

        /* create a dynamic 2D array
        Super clear video: https://www.youtube.com/watch?v=mGl9LO-je3o&list=PL43pGnjiVwgSSRlwfahAuIqoJ8TfDIlHq&index=6&ab_channel=CodeBeauty */

        landscapeTablePtr = new int *[rowNb]; // define the pointer to an array of int pointers for each row

        for (int row = 0; row < rowNb; row++) // for each row, create a pointer to an array of size columnNb
        {
            landscapeTablePtr[row] = new int[columnNb];
        }

        /* fill the table with info

        x/y cell coordinates */
        int r = 0; // initialise row counter
        int x = 0; // intialise x counter
        int y = 0;
        while (r < rowNb)
        {
            if (y > (Size - 1)) // we reach y = Size-1, reset y to 0 and increment x
            {
                y = 0;
                x++;
            }

            landscapeTablePtr[r][0] = r;
            landscapeTablePtr[r][1] = x;
            landscapeTablePtr[r][2] = y;

            /* match in xy coordinates to cellCode */
            XYcoordinates[r] = to_string(x) + ";" + to_string(y);
            cellCode[r] = r;

            /* debug : OK
            cout << "XYcoordinates[" << r << "] is " << XYcoordinates[r] << " cellCode[" << r << "] is " << cellCode[r] << endl;
            */

            y++;
            r++;
        }

        /* initialise resources to maximum */

        r = 0;
        while (r < rowNb)
        {

            for (int res = 0; res < resourceTypesNb; res++)
                landscapeTablePtr[r][resColumnStart + res] = MaxResource[res];

            r++;
        }

        /* initialise prey densities and catches to 0 */

        r = 0;
        while (r < rowNb)
        {
            for (int prey = 0; prey < (2 * preyTypesNb); prey++)
                landscapeTablePtr[r][preyColumnStart + prey] = 0;

            /* a specific number to check if the info is where expected: OK
            landscapeTablePtr[r][preyColumnStart + prey] = 1 + prey; */

            r++;
        }

        /* initialise predator densities to 0 */

        r = 0;
        while (r < rowNb)
        {
            for (int pred = 0; pred < predatorTypesNb; pred++)
                landscapeTablePtr[r][predatorColumnStart + pred] = 0;

            /* a specific number to check if the info is where expected: OK
            landscapeTablePtr[r][predatorColumnStart + pred] = 11 + pred; */

            r++;
        }
    }

    ~landscape() // free memory allocated to landscape table
    {

        /* free matching list memory */
        delete[] XYcoordinates;
        delete[] cellCode;

        /* free landscape tabel memory */
        for (int row = 0; row < rowNb; row++) // free memory allocated to each row
        {
            delete[] landscapeTablePtr[row];
        }

        delete[] landscapeTablePtr; // free the memory allocated to the array of pointer to each row

        /* debug :OK
        cout << "landscape destructor has been called successfully" << endl << endl;
        */
    }

    void resetCounts() // function to reset resources to maximum and counts to 0
    {

        /* debug
        cout << "reinitialising landscape ..." << endl
             << endl;
        */

        int r = 0;
        while (r < rowNb)
        {
            /* resetting counts to O */
            for (int col = preyColumnStart; col < columnNb; col++)
                landscapeTablePtr[r][col] = 0;

            r++;
        }
    }

    void resetResources() // function to reset resources to maximum and counts to 0
    {

        /* debug
        cout << "reinitialising resources to maximum ..." << endl
             << endl;
        */

        int r = 0;
        while (r < rowNb)
        {
            /* refill resources */
            for (int res = 0; res < resourceTypesNb; res++)
                landscapeTablePtr[r][resColumnStart + res] = MaxResource[0];

            r++;
        }
    }

    void getInfo() // function to cast out the landscape table and check if all good
    {
        /* cast out the column names */
        cout << "cellCode"
             << " "
             << "xCell"
             << " "
             << "yCell"
             << " ";
        for (int i = 0; i < resourceTypesNb; i++)
            cout << resourceTypes[i] << " ";
        for (int i = 0; i < preyTypesNb; i++)
            cout << preyTypes[i]
                 << " ";
        for (int i = 0; i < preyTypesNb; i++)
            cout << preyTypes[i] << "catches"
                 << " ";
        for (int i = 0; i < predatorTypesNb; i++)
            cout << predatorTypes[i] << " ";
        cout << endl;

        /* iterate through lines and column to cast out the values */
        int r = 0;
        while (r < rowNb)
        {
            for (int column = 0; column < columnNb; column++)
                cout << landscapeTablePtr[r][column] << " ";
            cout << endl;

            r++;
        }
        cout << endl;
    }

    void createResultsTable(string name) // creates a resultsTable
    {

        /* debug
        cout << "creating " << name << endl
             << endl;
        */

        /* write headers */
        resultsTable.open(name, ios::out);
        if (resultsTable.is_open())
        {
            resultsTable << "timeStep"
                         << ",";
            for (int i = 0; i < resourceTypesNb; i++)
                resultsTable << resourceTypes[i] << "amount"
                             << ",";
            for (int i = 0; i < preyTypesNb; i++)
                resultsTable << preyTypes[i] << "PopulationSize"
                             << ",";
            for (int i = 0; i < preyTypesNb; i++)
                resultsTable << preyTypes[i] << "catches"
                             << ",";
            for (int i = 0; i < predatorTypesNb; i++)
            {
                resultsTable << predatorTypes[i] << "PopulationSize";
                if (i == (predatorTypesNb - 1))
                    resultsTable << "\n";
                else
                    resultsTable << ",";
            }

            /* close file when finished */
            resultsTable.close();
        }
    }

    void createSnapshotTable(string name) // creates a snapshot file
    {

        /* debug
        cout << "creating " << name << endl
             << endl;
        */

        /* write headers */
        snapshotTable.open(name, ios::out);
        if (snapshotTable.is_open())
        {
            snapshotTable << "timeStep"
                          << ","
                          << "cellCode"
                          << ","
                          << "xCell"
                          << ","
                          << "yCell"
                          << ",";
            for (int i = 0; i < resourceTypesNb; i++)
                snapshotTable << resourceTypes[i] << "amount"
                              << ",";
            for (int i = 0; i < preyTypesNb; i++)
                snapshotTable << preyTypes[i]
                              << ",";
            for (int i = 0; i < preyTypesNb; i++)
                snapshotTable << preyTypes[i] << "catches"
                              << ",";
            for (int i = 0; i < predatorTypesNb; i++)
            {
                snapshotTable << predatorTypes[i];
                if (i == (predatorTypesNb - 1))
                    snapshotTable << "\n";
                else
                    snapshotTable << ",";
            }

            /* close file when finished */
            snapshotTable.close();
        }
    }

    void saveMeasures(string name, int ts) // write sum of results columns in the results file
    {
        /* debug
        cout << "appending " << name << endl
             << endl;
        */

        resultsTable.open(name, ios::app);
        if (resultsTable.is_open())
        {
            resultsTable << ts
                         << ",";

            /* write the sum of the measure columns */

            for (int i = 0; i < resourceTypesNb; i++)
                resultsTable << sumColumn(landscapeTablePtr, rowNb, resColumnStart + i)
                             << ",";
            for (int i = 0; i < preyTypesNb; i++)
                resultsTable << sumColumn(landscapeTablePtr, rowNb, preyColumnStart + i)
                             << ",";
            for (int i = 0; i < preyTypesNb; i++)
                resultsTable << sumColumn(landscapeTablePtr, rowNb, preyColumnStart + preyTypesNb + i)
                             << ",";
            for (int i = 0; i < predatorTypesNb; i++)
            {
                resultsTable << sumColumn(landscapeTablePtr, rowNb, predatorColumnStart + i);
                if (i == (predatorTypesNb - 1))
                    resultsTable << "\n";
                else
                    resultsTable << ",";
            }

            /* close file when finished */
            resultsTable.close();
        }
    }

    void snapshot(string name, int ts) // write all info columns in the snapshot file
    {
        /* debug
        cout << "appending " << name << endl
             << endl;
        */

        snapshotTable.open(name, ios::app);
        if (snapshotTable.is_open())
        {

            /* iterate through lines and column to cast out the values */
            int r = 0;
            while (r < rowNb)
            {
                snapshotTable << ts
                              << ",";
                for (int column = 0; column < columnNb; column++)
                {
                    if (column != (columnNb - 1))
                        snapshotTable << landscapeTablePtr[r][column] << ",";
                    else
                        snapshotTable << landscapeTablePtr[r][column] << "\n";
                }

                r++;
            }

            /* close file when finished */
            snapshotTable.close();
        }
    }
};

class animals
{
private:
    /* useful variables for memory allocation */
    int columnNb; // number of info stored in the table

    /* population level constants that might have different values according to animal types */
    int typeTag;           // a integer tag for each animal type
    float maxMove;         // max number of cells an animal can move by in each direction
    int reproductionCost;  // resources needed to reproduce
    int expectedOffspring; // max number of offspring when passing reproduction trial

protected:
    int landscapeSize;
    int maintenanceCost; // resources needed to survive a time step

    /* matching informations */
    int membersMatchingListsIndex;    //
    vector<int> dietLandscapeIndexes; // array containing their column index in the landscape table

public:
    /* population level variables */
    int currentPopulationSize; // current population size of the animal
    int initialDensity;
    int **populationTablePtr; //

    animals(int InitialDensity, int TypeTag, float MaxMove, int MaintenanceCost, int ReproductionCost, int LandscapeSize, int ExpectedOffspring, string *XYcoordinates, int *CellCodes) // initialise the constants shared by all animal types
    {

        initialDensity = InitialDensity;
        typeTag = TypeTag;
        reproductionCost = ReproductionCost;
        maintenanceCost = MaintenanceCost;
        expectedOffspring = ExpectedOffspring;
        landscapeSize = LandscapeSize;
        maxMove = ceil(MaxMove * (float)landscapeSize);

        if (maxMove >= landscapeSize)
            cout << "Warning animal's moving range is larger than the landscape size, can mess cell coordinates when moving" << endl;

        currentPopulationSize = initialDensity; // initialise current population size variable

        membersMatchingListsIndex = getMemberIndexFromTag(typeTag);
        dietLandscapeIndexes = getDietLandscapeIndexes(membersMatchingListsIndex);

        /* create animals table */
        columnNb = 7; // x - y - cellCode - resourceConsumed - deadOrAlive - offspring - age

        populationTablePtr = new int *[currentPopulationSize]; // define the pointer to an array of int pointers for each row

        for (int row = 0; row < currentPopulationSize; row++) // for each row, create a pointer to an integer array of size columnNb
        {
            populationTablePtr[row] = new int[columnNb];
        }

        /* add the animals with their info */
        int r = 0;
        while (r < initialDensity)
        {
            /* debug
            cout << "float(landscapeSize) is " << float(landscapeSize) << endl;
            cout << "random number between zero and landsape size without imposing data type is " << randomNumberGenerator(0, landscapeSize) << endl;
            cout << "random number between zero and landsape size imposing float(size) is " << randomNumberGenerator(0, float(landscapeSize)) << endl;
            cout << "the value that goes in the table is " << int(randomNumberGenerator(0, float(landscapeSize))) << endl;
            */

            /* assign random coordinates between 0 and Size */
            populationTablePtr[r][0] = int(randomNumberGeneratorAny(0, landscapeSize - 1));
            populationTablePtr[r][1] = int(randomNumberGeneratorAny(0, landscapeSize - 1));

            /* update populationTablePtr with prey info - columns indexes hard coded NOT GOOD */
            populationTablePtr[r][2] = getCellCode(XYcoordinates, CellCodes, landscapeSize, populationTablePtr[r][0], populationTablePtr[r][1]);
            populationTablePtr[r][3] = 0; // initialise resource consumed
            populationTablePtr[r][4] = 1; // the animal is alive
            populationTablePtr[r][5] = 0; // the animal has no offspring yet
            populationTablePtr[r][6] = 0; // initial population is of age 0

            r++;
        }
    }

    ~animals() // free memory allocated to animals table
    {

        for (int row = 0; row < currentPopulationSize; row++) // free memory allocated to each row
        {
            delete[] populationTablePtr[row];
        }

        delete[] populationTablePtr; // free the memory allocated to the array of pointer to each row

        /* debug : OK
        cout << "An animal destructor has been called successfully" << endl << endl;
        */
    }

    void randomMove(string *XYcoordinates, int *CellCodes)
    {
        /* debug
        cout << memberTypes[membersMatchingListsIndex] << " are moving at random" << endl
             << endl;
        */

        /* iterate through individuals */
        int ind = 0;

        while (ind < currentPopulationSize)
        {
            /* update position with random number in moving range */
            int xCell = populationTablePtr[ind][0];
            int yCell = populationTablePtr[ind][1];

            /* debug
            cout << "individual #" << ind << " was on cell [" << xCell << ";" << yCell << "]" << endl;
            */

            double xRand = randomNumberGeneratorAny(-1 * maxMove, maxMove);
            double yRand = randomNumberGeneratorAny(-1 * maxMove, maxMove);

            /* debug
            cout << "individual #" << ind << " moved by [" << xCell << " + (" << xRand << ");" << yCell << " + (" << yRand << ")]" << endl;
            */

            xCell += xRand;
            yCell += yRand;

            /* check boundaries */
            populationTablePtr[ind][0] = transmissiveBoundaries(xCell, landscapeSize); // xCell
            populationTablePtr[ind][1] = transmissiveBoundaries(yCell, landscapeSize); // yCell

            /* debug
            cout << "after transmissive boundaries correction: [" << xCell << ";" << yCell << "]" << endl;
            */

            populationTablePtr[ind][2] = getCellCode(XYcoordinates, CellCodes, landscapeSize, populationTablePtr[ind][0], populationTablePtr[ind][1]); // update cell code

            ind++; // next individual
        }
    }

    void survivalTrial()
    {
        /* debug
        cout << memberTypes[membersMatchingListsIndex] << " survival trials" << endl
             << endl;
        */

        /* iterate through individuals */
        int ind = 0;
        int zz = 0; // deaths counter

        while (ind < currentPopulationSize)
        {
            if (populationTablePtr[ind][4] == 1) // if individual is alive
            {
                /* if resonsume < maintenance cost -> change deadOrAlive to 0 else remove maintenance cost from the individual's resource pool */
                if (populationTablePtr[ind][3] < maintenanceCost)
                {
                    populationTablePtr[ind][4] = 0;
                    zz += 1;
                }
                else
                    populationTablePtr[ind][3] -= maintenanceCost;
            }

            ind++; // next individual
        }

        /* debug
        if (zz > 0)
            cout << zz << " " << memberTypes[membersMatchingListsIndex] << " had not enough resources to survive" << endl
                 << endl;
        */
    }

    void reproductionTrial()
    {
        /* debug
        cout << memberTypes[membersMatchingListsIndex] << " reproduction trials" << endl
             << endl;
        */

        /* iterate through individuals */
        int ind = 0;

        while (ind < currentPopulationSize)
        {
            if (populationTablePtr[ind][4] == 1 && populationTablePtr[ind][3] >= reproductionCost) // if individual is alive and has enough resources to reproduce
            {
                populationTablePtr[ind][5] = randomSampleFromPoissonDist(expectedOffspring); // even if it has the resource it might not reproduce (and thus not paying the cost)

                /* debug 
                cout << "animal of type " << typeTag << " #" << ind << "has enough resources to reproduce. Offspring nb: " << populationTablePtr[ind][5] << endl
                     << endl;
                */

                if (populationTablePtr[ind][5] > 0)
                    populationTablePtr[ind][3] -= reproductionCost; // if at least one offspring subtract reproduction cost from resource pool
            }

            ind++; // next individual
        }
    }

    void updatePopulationTable(int timeStep, bool debug) // updates current pop size, creates a new table with updated information, replaces popTable with the updated version, reset offspring to 0, deletes the older version
    {

        /* update current population size with deaths and offspring */
        int oldPopulationSize = currentPopulationSize;                                                  // store current population size for later purpose                                                        // store previous population size
        int deadInds = currentPopulationSize - sumColumn(populationTablePtr, currentPopulationSize, 4); // number of dead ind is pop size - sum of alive individuals
        int newInds = sumColumn(populationTablePtr, currentPopulationSize, 5);                          // sum of all offspring produces

        currentPopulationSize += newInds - deadInds; // update current population size/* debug */

        // /* make warnings if low pop */
        // if (initialDensity != 0)
        // {
        //     if (float(currentPopulationSize) <= 0.1 * float(initialDensity))
        //         cout << "time step #" << timeStep << ": " << memberTypes[membersMatchingListsIndex] << "'s population is under 10 percent of its initial size" << endl
        //              << endl;
        //     if (currentPopulationSize == 0)
        //         cout << "time step #" << timeStep << ": " << memberTypes[membersMatchingListsIndex] << "'s population has gone extinct" << endl
        //              << endl;
        // }

        if (debug == true)
        {
            cout << memberTypes[membersMatchingListsIndex] << ": " << newInds << " offspring and " << deadInds << " deaths" << endl
                 << endl;
        }

        /* allocate memory to a new table of the right size */
        int **newTablePtr = new int *[currentPopulationSize];
        for (int row = 0; row < currentPopulationSize; row++)
        {
            newTablePtr[row] = new int[columnNb];
        }

        /* replicate all living individuals from the previous list in the new table and add offspring */
        int newTabRow = 0; // initialise a row counter for the new table

        for (int oldPopRow = 0; oldPopRow < oldPopulationSize; oldPopRow++) // iterate through previous pop table
        {
            int indDoAstatus = populationTablePtr[oldPopRow][4]; // get dead or alive status
            int indOffspring = populationTablePtr[oldPopRow][5]; // get offspring number
            int indAge = populationTablePtr[oldPopRow][6];       // get age

            if (indDoAstatus == 1 && newTabRow < currentPopulationSize) // if individual is alive and we don't override the new table's dimensions
            {
                for (int col = 0; col < 5; col++)
                    newTablePtr[newTabRow][col] = populationTablePtr[oldPopRow][col]; // copy x y position cellCode resource pool and DoA status into new table

                newTablePtr[newTabRow][5] = 0;          // reset offspring to zero in the new table
                newTablePtr[newTabRow][6] = indAge + 1; // individual is one time step older in the new table

                newTabRow++; // increment new table row counter
            }

            /* add all the offspring at the cell of their parent */
            if (indOffspring > 0 && newTabRow < currentPopulationSize)
            {
                for (int i = 0; i < indOffspring; i++) // for each offspring
                {
                    for (int col = 0; col < 3; col++) // copy only position variables
                        newTablePtr[newTabRow][col] = populationTablePtr[oldPopRow][col];

                    newTablePtr[newTabRow][3] = 0; // new individual has not consumed any resource yet
                    newTablePtr[newTabRow][4] = 1; // new individual is alive
                    newTablePtr[newTabRow][5] = 0; // new individual has not reproduced yet
                    newTablePtr[newTabRow][6] = 0; // new individual is of age 0

                    newTabRow++; // increment new table row counter
                }
            }
        }

        /* older table not needed anymore, delete it */
        for (int row = 0; row < oldPopulationSize; row++)
        {
            delete[] populationTablePtr[row];
        }
        delete[] populationTablePtr;

        populationTablePtr = newTablePtr;

        newTablePtr = NULL;
    }

    void measureDensity(int **LandscapeTable)
    {

        /* debug
        cout << "measuring " << memberTypes[membersMatchingListsIndex] << " densities" << endl
             << endl;
        */

        /* get animal's column index in landscapeTable */
        int colIndex = indexInLandscape[membersMatchingListsIndex];

        /* resetting densities to 0 */
        for (int landRow = 0; landRow < (landscapeSize * landscapeSize); landRow++)
        {
            LandscapeTable[landRow][colIndex] = 0;
        }

        /* iterate through population table */
        int ind = 0;
        while (ind < currentPopulationSize)
        {
            /* increment landscape table at the right position */
            int rowIndex = populationTablePtr[ind][2];
            LandscapeTable[rowIndex][colIndex] += 1;

            /* debug
            cout << "incrementing cell [" << rowIndex << "] [" << colIndex << "]" << endl
                 << endl;
            */

            ind++; // next individual
        }
    }

    void getInfo() // function to cast out the animalsTable and check if all good
    {

        /* cast out the column names */
        cout << memberTypes[membersMatchingListsIndex] << "'s population table" << endl;
        cout << "xCell"
             << " "
             << "yCell"
             << " "
             << "cellCode"
             << " "
             << "resourcesConsumed"
             << " "
             << "deadOrAlive"
             << " "
             << "offspring"
             << " "
             << "age"
             << endl;

        /* iterate through lines and column to cast out the values */
        int r = 0;
        while (r < currentPopulationSize)
        {
            for (int column = 0; column < columnNb; column++)
                cout << populationTablePtr[r][column] << " ";
            cout << endl;

            r++;
        }
        cout << endl;
    }
};

class prey : public animals // preys are one type of animal object, they share the same charasteristics but also have specificities
{

private:
    int maxConsumption; // max resource units the prey can consume in a time step

public:
    prey(int preyInitialDensity, int preyTypeTag, float preyMaxMovingDistance, int preyMaintenanceCost, int preyReproductionCost, int LandscapeSize, int preyExpectedOffspring, string *XYcoordinates, int *CellCodes) : animals(preyInitialDensity, preyTypeTag, preyMaxMovingDistance, preyMaintenanceCost, preyReproductionCost, LandscapeSize, preyExpectedOffspring, XYcoordinates, CellCodes) //
    {
    }

    void assignPreyVariables(int preyMaxConsumption)
    {
        maxConsumption = preyMaxConsumption;
    }

    void feed(int **landscapeTable, bool debug)
    {

        if (debug == true)
            cout << memberTypes[membersMatchingListsIndex] << " are feeding" << endl
                 << endl;

        /* shuffle the order of the individuals */
        vector<int> shuffledPop = shuffleOrder(currentPopulationSize);

        int ind = 0;                     // initialise individual counter
        while (ind < shuffledPop.size()) // iterate through individuals
        {

            int rowIndex = shuffledPop[ind]; // get shuffled row index

            int indCellCode = populationTablePtr[rowIndex][2]; // get individual's cellCode

            /* iteration through diet and feed */
            int dietSize = dietLandscapeIndexes.size();

            if (debug == true)
                cout << "prey #" << rowIndex << ", located on cell " << indCellCode << " feeds on " << dietSize << " kind of resources" << endl
                     << endl;

            for (int i = 0; i < dietSize; i++)
            {
                /* get resource amount */
                int resourceColIndex = dietLandscapeIndexes[i]; // resource column index in landscape table

                int resourceAvailable = landscapeTable[indCellCode][resourceColIndex];

                if (debug == true)
                    cout << "feeding on resource situated col " << resourceColIndex << " of which " << resourceAvailable << "are left" << endl
                         << endl;

                /* consumption */
                int cons = 0;                            // initialise a consumption variable
                if (resourceAvailable >= maxConsumption) // if there is more than maxConsumption rate, consume max
                {
                    cons = randomNumberGeneratorAny(0, maxConsumption); // generate a random number between 0 and maxConsume
                    // landscapeTable[indCellCode][resourceColIndex] -= maxConsumption; // subtract consumption to landscape cell
                    // populationTablePtr[rowIndex][3] += maxConsumption;
                }
                else // else consume what is left (should also work if resourceAvailable=0)
                {
                    cons = randomNumberGeneratorAny(0, resourceAvailable); // generate a random number between 0 and maxConsume
                    // landscapeTable[indCellCode][resourceColIndex] -= resourceAvailable;
                    // populationTablePtr[rowIndex][3] += resourceAvailable;
                }
                landscapeTable[indCellCode][resourceColIndex] -= cons; // subtract consumption to landscape cell
                populationTablePtr[rowIndex][3] += cons;

                if (debug == true)
                    // cout << "prey consumed" << maxConsumption << "resources, " << landscapeTable[indCellCode][resourceColIndex] << "left on this cell, prey resource pool is now" << populationTablePtr[rowIndex][3] << endl
                    cout << "prey consumed" << cons << "resources, " << landscapeTable[indCellCode][resourceColIndex] << "left on this cell, prey resource pool is now" << populationTablePtr[rowIndex][3] << endl
                         << endl;
            }

            ind++; // next individual
        }
    }

    void countCatches(int **LandscapeTable, bool debug)
    {

        int catchesColumn = indexInLandscape[membersMatchingListsIndex] + preyTypesNb; // get the catches column index in landscape table

        /* debug */
        if (debug == true)
        {
            // cout << "counting " << memberTypes[membersMatchingListsIndex] << " catches" << endl
            //     << "catches counting column is " << catchesColumn << endl
            cout << memberTypes[membersMatchingListsIndex] << ": " << sumColumn(LandscapeTable, landscapeSize * landscapeSize, catchesColumn) << " catches" << endl
                 << endl;
        }

        /* found the cells where prey have been caught by predators in landscape table */
        for (int landRow = 0; landRow < (landscapeSize * landscapeSize); landRow++) // iterate through landscape rows
        {

            if (LandscapeTable[landRow][catchesColumn] > 0) // if there was at least 1 catch
            {

                int catches = LandscapeTable[landRow][catchesColumn]; // store the number of catches

                /* debug */
                if (debug == true)
                    cout << "there was " << catches << " " << memberTypes[membersMatchingListsIndex] << " caught on cell " << landRow << endl
                         << endl;

                vector<int> shuffledPop = shuffleOrder(currentPopulationSize); // shuffle the indexes of prey table

                for (int ind = 0; ind < shuffledPop.size(); ind++) // iterate through population table individuals
                {
                    int popRow = shuffledPop[ind];                   // get shuffled row index
                    int indCellCode = populationTablePtr[popRow][2]; // get the cell it is on
                    int DoAstatus = populationTablePtr[popRow][4];   // DoA?

                    if (indCellCode == landRow && DoAstatus == 1) // if it is on the same cell, and is not already dead for some reason
                    {

                        /* debug */
                        if (debug == true)
                            cout << memberTypes[membersMatchingListsIndex] << " number " << popRow << " is on cell " << landRow << " with DoA status " << DoAstatus << endl
                                 << endl;

                        populationTablePtr[popRow][4] = 0; // update individual's dead or alive status
                        catches--;                         // one less catch to take into account

                        /* debug */
                        if (debug == true)
                        {
                            cout << memberTypes[membersMatchingListsIndex] << " number " << popRow << "'s DoA status is now " << populationTablePtr[popRow][4] << endl;
                            // if (populationTablePtr[popRow][3] < maintenanceCost)
                            //     cout << "but it did not have enough resources to survive anyway" << endl;
                            cout << catches << " " << memberTypes[membersMatchingListsIndex] << " catches left to count on cell " << landRow << endl
                                 << endl;
                        }

                        if (catches == 0) // if all catches have been taken into account break out of the innermost loop
                            if (debug == true)
                                cout << "breaking out of loop over population table" << endl
                                     << endl;
                        break;
                    }
                }
            }
        }
    }
};

class predator : public animals // predators are another type of animal object
{

protected:
    int predMaxConsumption;
    int dailyPredMaxConsumption;
    vector<int> conversionRates; // conversion of each prey catch in resources
    vector<int> maxCatches;      // max number of catches per day and per prey
    vector<float> catchProbas;   // catch probability of each prey in diet

public:
    predator(int initialDensity, int typeTag, float maxMovingDistance, int predatorMaintenanceCost, int predatorReproductionCost, int LandscapeSize, int predatorExpectedOffspring, string *XYcoordinates, int *CellCodes) : animals(initialDensity, typeTag, maxMovingDistance, predatorMaintenanceCost, predatorReproductionCost, LandscapeSize, predatorExpectedOffspring, XYcoordinates, CellCodes) //
    {
    }

    void getDietInfo()
    {
        /* set maxConsume */
        predMaxConsumption = 3 * maintenanceCost;

        /* set daily maxConsume */
        dailyPredMaxConsumption = predMaxConsume[0];

        // /* get preys conversion rates in dietsTable */
        // for (int i = 0; i < memberMatchingListsSize; i++)
        // {
        //     if (dietsTable[i][membersMatchingListsIndex] > 0)
        //         conversionRates.push_back(dietsTable[i][membersMatchingListsIndex]);
        // }

        /* get preys conversion rates from parameters A AMELIORER */
        for (int i = 0; i < predConvRate.size(); i++)
        {
            conversionRates.push_back(predConvRate[i]);
        }

        /* compute max catches */
        for (int i = 0; i < conversionRates.size(); i++)
        {
            int res = ceil(float(3 * dailyPredMaxConsumption) / float(conversionRates[i]));
            // int res = float(maintenanceCost) * 3 / (float(conversionRates[i]) * float(freqSurv));
            maxCatches.push_back(res);
        }

        /* debug */
        cout << "Max catches per day should be ceil " << 3 * dailyPredMaxConsumption << "/" << conversionRates[0] << ";" << 3 * dailyPredMaxConsumption << "/" << conversionRates[1] << " and are " << maxCatches[0] << " " << maxCatches[1] << endl
             << endl;

        catchProbas.push_back(predCatchProba[0]);
        catchProbas.push_back(predCatchProba[1]);
    }

    void hunt(int **LandscapeTable, bool debug)
    {

        /* debug */
        if (debug == true)
        {
            cout << memberTypes[membersMatchingListsIndex] << " are hunting" << endl
                 << endl;
        }

        /* shuffle the order of the individuals */
        vector<int> shuffledPop = shuffleOrder(currentPopulationSize);

        int ind = 0;                        // initialise individual counter
        while (ind < currentPopulationSize) // iterate through individuals
        {
            int rowIndex = shuffledPop[ind]; // get shuffled row index

            int indCellCode = populationTablePtr[rowIndex][2]; // get individual's cellCode

            /* debug */
            if (debug == true)
            {
                cout << "predator number " << rowIndex << " is hunting" << endl
                     << "it is on cell " << indCellCode << endl
                     << endl;
            }

            /* generate vectors preys info */
            vector<int> preysIndexes;
            vector<int> preysConversionRates;
            vector<int> preysLandscapeIndexes;
            vector<int> preysMaxCatches;
            vector<float> preysCatchProbas;

            for (int i = 0; i < conversionRates.size(); i++)
            {
                preysIndexes.push_back(i);
            }

            if (predSpecific[0] == false)
            {
                if (debug == true)
                    cout << "predator is generalist: randomizing prey types in diet." << endl;

                /* generate vector of indexes corresponding to conversion rates */
                vector<int> shuffledIndexes = preysIndexes;
                vector<int> shuffledConversionRates;
                vector<int> shuffledpreysLandscapeIndexes;
                vector<int> shuffledMaxCatches;
                vector<float> shuffledCatchProbas;

                /* shuffle indexes and update conversionRates and maxCatches accordingly */
                random_shuffle(shuffledIndexes.begin(), shuffledIndexes.end());
                for (int i = 0; i < shuffledIndexes.size(); i++)
                {
                    int index = shuffledIndexes[i];
                    shuffledConversionRates.push_back(conversionRates[index]);
                    shuffledpreysLandscapeIndexes.push_back(dietLandscapeIndexes[index]);
                    shuffledMaxCatches.push_back(maxCatches[index]);
                    shuffledCatchProbas.push_back(catchProbas[index]);
                }

                if (predOportunistic[0] == true)
                {
                    if (debug == true)
                        cout << "predator is opportunistic: ranking prey types per resources/catch after randomizing them." << endl;

                    /* sort the diet by conversion rate while keeping index info: if all the same, still random, otherwise most nourrishing prey first */
                    vector<int> sortedIndexes = shuffledIndexes;
                    vector<int> sortedConversionRates = shuffledConversionRates;
                    vector<int> sortedPreyLandscapeIndexes;
                    vector<int> sortedMaxCatches;
                    vector<float> sortedCatchProbas;

                    /* sorting according to conversion rate while keeping the indexes */
                    int row = 0; // initiate row count

                    while (row < (sortedConversionRates.size() - 1)) // until we reach the line before last
                    {
                        if (sortedConversionRates[row] >= sortedConversionRates[row + 1]) // if the focal element is greater than the next
                        {
                            row++; // leave as is and go to next line
                        }
                        else // if not switch positions and restart to the first line
                        {
                            int co1 = sortedConversionRates[row];
                            int co2 = sortedConversionRates[row + 1];
                            int in1 = sortedIndexes[row];
                            int in2 = sortedIndexes[row + 1];

                            sortedConversionRates[row] = co2;
                            sortedConversionRates[row + 1] = co1;
                            sortedIndexes[row] = in2;
                            sortedIndexes[row + 1] = in1;

                            row = 0; // restart row count
                        }
                    }

                    for (int i = 0; i < sortedIndexes.size(); i++)
                    {
                        int index = sortedIndexes[i];
                        sortedPreyLandscapeIndexes.push_back(dietLandscapeIndexes[index]);
                        sortedMaxCatches.push_back(maxCatches[index]);
                        sortedCatchProbas.push_back(catchProbas[index]);
                    }

                    /* update vectors preys info */
                    preysIndexes = sortedIndexes;
                    preysConversionRates = sortedConversionRates;
                    preysLandscapeIndexes = sortedPreyLandscapeIndexes;
                    preysMaxCatches = sortedMaxCatches;
                    preysCatchProbas = sortedCatchProbas;
                } // end if opportunistic 
                else 
                {
                    /* update vectors preys info */
                    preysIndexes = shuffledIndexes;
                    preysConversionRates = shuffledConversionRates;
                    preysLandscapeIndexes = shuffledpreysLandscapeIndexes;
                    preysMaxCatches = shuffledMaxCatches;
                    preysCatchProbas = shuffledCatchProbas;
                }
            }
            else
            {
                if (debug == true)
                    cout << "Predator is specific: always hunts in diets order" << endl
                         << endl;

                /* update vectors preys info */
                preysConversionRates = conversionRates;
                preysLandscapeIndexes = dietLandscapeIndexes;
                preysMaxCatches = maxCatches;
                preysCatchProbas = catchProbas;
            }

            /* iterate through prey columns and while catches < maxCatches and
               that there are prey available, catch them */

            int i = 0;                                      // initiate index for conversion rate
            int dailyCons = 0;                              // initiate a tracker for pred consumption (in resources)
            int predCons = populationTablePtr[rowIndex][3]; // get predator resource pool

            while (i < preysConversionRates.size())
            {
                /* get preys info */
                int densColumn = preysLandscapeIndexes[i];
                int catchColumn = densColumn + preyTypesNb;
                int dens = LandscapeTable[indCellCode][densColumn];
                int convRate = preysConversionRates[i];
                int maxCatch = preysMaxCatches[i];
                float catchProb = preysCatchProbas[i];

                int catches = 0; // initialise a catch counter for this particular prey

                if (debug == true)
                {
                    cout << "searching for prey situated column " << densColumn << " in landscape table" << endl
                         << "its catch counting column in landscape table is number " << catchColumn << endl
                         << "its density on this cell is " << dens << endl
                         << endl;
                }

                while (dens > 0 && catches < maxCatch && dailyCons < dailyPredMaxConsumption && predCons < predMaxConsumption)
                {
                    /* sample a random number between 0 and 1 in an uniform distribution */
                    float randomNb = randomNumberGenerator0to1(0, 1);

                    if (debug == true)
                        cout << "random sample is: " << randomNb << " ; proba is: " << catchProb << endl
                             << endl;

                    /* compare to the catch probability, if < catch, if > fail */
                    if (randomNb < catchProb)
                    {
                        LandscapeTable[indCellCode][catchColumn] += 1;  // increment corresponding catch cell in landscape table
                        LandscapeTable[indCellCode][densColumn] -= 1;   // decrement density on the cell such that another predator cannot catch more individuals than there actually are on the cell
                        populationTablePtr[rowIndex][3] += convRate;    // increment predator resource pool
                                                                        //
                        predCons = populationTablePtr[rowIndex][3];     // update predCons variable
                        dailyCons += convRate;                          // update daily consuption variable
                        catches++;                                      // increment catches counter
                        dens = LandscapeTable[indCellCode][densColumn]; // update prey's density on this cell

                        /* debug */
                        if (debug == true)
                        {
                            cout << "a " << memberTypes[membersMatchingListsIndex] << " caught a prey (situated column " << densColumn << " in landscape table) on cell " << indCellCode << endl
                                 << "its resource pool now got " << predCons << " in it" << endl;

                            if (dens == 0)
                                cout << "No more prey situated column " << densColumn << " in landscape table on cell " << indCellCode << "." << endl
                                     << endl;
                            else if (catches >= maxCatch)
                                cout << memberTypes[membersMatchingListsIndex] << "situated on cell " << indCellCode << " has eaten enough prey situated column " << densColumn << " in landscape table for this time step." << endl
                                     << endl;
                            else if (dailyCons >= dailyPredMaxConsumption)
                                cout << memberTypes[membersMatchingListsIndex] << "situated on cell " << indCellCode << " has eaten enough for this time step." << endl
                                     << endl;
                            else if (predCons >= predMaxConsumption)
                                cout << memberTypes[membersMatchingListsIndex] << "situated on cell " << indCellCode << " has eaten enough for this moving+feeding sequence." << endl
                                     << endl;
                        }
                    }
                    else
                    {
                        catches++; // increase catches to prevent preadators to continue until they catch something

                        if (debug == true)
                            cout << "predator catch attempt failed" << endl
                                 << endl;
                    }

                } // end of while loop over conditions for hunting

                i++; // increment prey index

            } // end of while loop over preys

            ind++; // next individual

        } // end of while loop over predators

    } // end of hunt function
};    // end of predator class

/* ------------------------------ Main program ------------------------------ */

int main(int argc, char **argv)
{

    /* ---- assign values to the variables from the command line ---- */

    int p = 1; // initialize parameter counter

    /* give a name to the simulation */
    simulationName = argv[p];

    p++; // increment parameter counter

    /* landscape variables */
    worldSize = atoi(argv[p]);

    p++; // increment parameter counter

    resourceTypesNb = atoi(argv[p]);

    p++;

    for (int i = 0; i < resourceTypesNb; i++)
        resourceTypes.push_back("resource" + to_string(i + 1)); // res1..n

    maxResources.push_back(atoi(argv[p]));

    p++;

    maxResources.push_back(atoi(argv[p]));

    p++;

    /* prey variables */
    preyTypesNb = atoi(argv[p]);

    p++;

    preyInitialDensities.push_back(atoi(argv[p]));

    p++;

    preyInitialDensities.push_back(atoi(argv[p]));

    p++;

    for (int i = 0; i < preyTypesNb; i++)
        preyTypes.push_back("prey" + to_string(i + 1)); // prey1..n

    preyMaxMove.push_back(atof(argv[p]));

    p++;

    preyMaxMove.push_back(atof(argv[p]));

    p++;

    preyMaxConsume.push_back(atoi(argv[p]));

    p++;

    preyMaxConsume.push_back(atoi(argv[p]));

    p++;

    preyMaintenanceCost.push_back(atoi(argv[p]));

    p++;

    preyMaintenanceCost.push_back(atoi(argv[p]));

    p++;

    preyExpectedOffspring.push_back(atoi(argv[p]));

    p++;

    preyExpectedOffspring.push_back(atoi(argv[p]));

    p++;

    preyReproCost.push_back(atoi(argv[p]));

    p++;

    preyReproCost.push_back(atoi(argv[p]));

    p++;
    

    /* predator variables */
    predatorTypesNb = atoi(argv[p]);

    p++;

    for (int i = 0; i < predatorTypesNb; i++)
        predatorTypes.push_back("predator" + to_string(i + 1)); // res1..n

    predInitialDensities.push_back(atoi(argv[p]));

    p++;

    predMaxMove.push_back(atof(argv[p]));

    p++;

    predMaxConsume.push_back(atoi(argv[p]));

    p++;

    predMaintenanceCost.push_back(atoi(argv[p]));

    p++;

    predExpectedOffspring.push_back(atoi(argv[p]));

    p++;

    predReproCost.push_back(atoi(argv[p]));

    p++;

    predIntro.push_back(atoi(argv[p]));

    p++;

    // predAsymm.push_back(atof(argv[p]));

    // p++;

    predCatchProba.push_back(atof(argv[p]));

    p++;

    predCatchProba.push_back(atof(argv[p]));
    
    p++;

    predConvRate.push_back(atoi(argv[p]));

    p++;

    predConvRate.push_back(atoi(argv[p]));
    
    p++;

    predOportunistic.push_back(atoi(argv[p]));
    
    p++;

    predSpecific.push_back(atoi(argv[p]));

    p++;

    /* time variables */
    timeMaxi = atoi(argv[p]); // simulation time

    p++;

    freqRepr = atoi(argv[p]);

    p++;

    freqSurv = atoi(argv[p]);

    p++;

    freqRfll = atoi(argv[p]); // let the animals feed for a while before "daily" death trial

    p++;

    /* assessment frequency variables */
    freqResu = atoi(argv[p]);

    p++;

    freqSnap = atoi(argv[p]);

    p++;

    /* seed */
    randomSeed = atoi(argv[p]);

    /* ---- construct matching structures ---- */

    memberMatchingListsSize = resourceTypesNb + preyTypesNb + predatorTypesNb;

    assignTagsIndexes(); // creates individuals' matching tables. NEED 3 delete[] : OK

    makeDietsTable(true); // NEED 1 delete[] : OK

    srand(randomSeed); // set random generator seed with instant time

    /* ---- initiate world ---- */

    /* contruct and create landscape */
    landscape world(worldSize, maxResources); // NOT GOOD add a variable for different max for each resource

    /* check if everything is where expected: OK
    world.getInfo();
    */

    /* ---- construct animals populations ---- */

    prey *prey1 = new prey(preyInitialDensities[0], 201, preyMaxMove[0], preyMaintenanceCost[0], preyReproCost[0], worldSize, preyExpectedOffspring[0], world.XYcoordinates, world.cellCode); // construct prey1 population = assigning values to the constants, intitialise some variables, compute others
    prey1->assignPreyVariables(preyMaxConsume[0]);

    prey *prey2 = new prey(preyInitialDensities[1], 202, preyMaxMove[1], preyMaintenanceCost[1], preyReproCost[1], worldSize, preyExpectedOffspring[1], world.XYcoordinates, world.cellCode); // construct prey1 population = assigning values to the constants, intitialise some variables, compute others
    prey2->assignPreyVariables(preyMaxConsume[1]);

    /* create prey pointer group */
    prey *preys[2] = {prey1, prey2};

    predator *pred1 = new predator(predInitialDensities[0], 301, predMaxMove[0], predMaintenanceCost[0], predReproCost[0], worldSize, predExpectedOffspring[0], world.XYcoordinates, world.cellCode);
    pred1->getDietInfo();

    /* if more than one predator
    predator *predators[predatorTypesNb] = {};
    pred1->getInfo();
    */

    /* check create function : OK*/
    // prey1->getInfo();
    // prey2->getInfo();
    // pred1->getInfo();

    /* ---- create results and snapshot csv files ---- */

    /* file names */
    string resultsTableName = simulationName + "-ResultsTable.csv";
    string snapshotTableName = simulationName + "-SnapshotTable.csv";

    world.createResultsTable(resultsTableName);
    world.createSnapshotTable(snapshotTableName);

    /* ---- start simulation ---- */

    /* debug */
    cout << "starting simulation over " << timeMaxi << " time steps" << endl
         << endl;

    int timeStep = 0; // initialise time step variable

    while (timeStep <= timeMaxi)
    {
        /* debug
        cout << "time step " << timeStep << endl
             << endl;
        */

        if (timeStep > 0)
        {
            /* ---- life events ---- */

            /* -- moving -- */

            /* preys */
            for (int i = 0; i < preyTypesNb; i++)
                preys[i]->randomMove(world.XYcoordinates, world.cellCode);

            // prey1->getInfo();
            // prey2->getInfo();

            /* predators */
            if (timeStep >= predIntro[0])
                pred1->randomMove(world.XYcoordinates, world.cellCode);

            // for (int i = 0; i < predatorTypesNb; i++)
            // {
            //     predators[i]->randomMove(world.XYcoordinates, world.cellCode);
            // }

            // pred1->getInfo();

            /* -- measure densities -- */

            /* preys */
            for (int i = 0; i < preyTypesNb; i++)
                preys[i]->measureDensity(world.landscapeTablePtr);

            // prey1.getInfo();
            // prey2.getInfo();

            /* predators */
            if (timeStep >= predIntro[0])
                pred1->measureDensity(world.landscapeTablePtr);

            // for (int i = 0; i < predatorTypesNb; i++)
            // {
            //     predators[i]->measureDensity(world.landscapeTablePtr);
            // }

            // pred1->getInfo();

            /* -- feeding -- */

            /* preys */
            for (int i = 0; i < preyTypesNb; i++)
                preys[i]->feed(world.landscapeTablePtr, false);

            // prey1->getInfo();
            // prey2->getInfo();

            /* predators */
            if (timeStep >= predIntro[0])
                pred1->hunt(world.landscapeTablePtr, false);

            // for (int i = 0; i < predatorTypesNb; i++)
            // {
            //     predators[i]->hunt(world.landscapeTablePtr, false);
            // }

            // pred1->getInfo();

            /* -- counting catches -- */

            for (int i = 0; i < preyTypesNb; i++)
                preys[i]->countCatches(world.landscapeTablePtr, false);

            // prey1.getInfo();
            // prey2.getInfo();

            /* -- surviving -- */

            if (timeStep % freqSurv == 0)
            {
                /* preys */
                for (int i = 0; i < preyTypesNb; i++)
                    preys[i]->survivalTrial();

                // prey1.getInfo();
                // prey2.getInfo();

                /* predators */
                if (timeStep >= predIntro[0])
                    pred1->survivalTrial();

                // for (int i = 0; i < predatorTypesNb; i++)
                // {
                //     predators[i]->survivalTrial();
                // }

                // pred1->getInfo();
            }

            /* -- reproducing -- */
            if (timeStep % freqRepr == 0)
            {
                /* preys */
                for (int i = 0; i < preyTypesNb; i++)
                    preys[i]->reproductionTrial();

                // prey1.getInfo();
                // prey2.getInfo();

                /* predators */
                if (timeStep >= predIntro[0])
                    pred1->reproductionTrial();

                // for (int i = 0; i < predatorTypesNb; i++)
                // {
                //     predators[i]->reproductionTrial();
                // }

                // pred1->getInfo();
            }
        }

        /* -- update animals dynamic arrays with birth catches and deaths -- */

        /* preys */
        for (int i = 0; i < preyTypesNb; i++)
            preys[i]->updatePopulationTable(timeStep, false);

        // prey1->getInfo();
        // prey2->getInfo();

        /* predators */
        if (timeStep >= predIntro[0])
            pred1->updatePopulationTable(timeStep, false);

        // for (int i = 0; i < predatorTypesNb; i++)
        // {
        //     predators[i]->updatePopulationTable(false);
        // }

        // pred1->getInfo();

        /* -- measure densities -- */

        /* preys */
        for (int i = 0; i < preyTypesNb; i++)
            preys[i]->measureDensity(world.landscapeTablePtr);

        // prey1.getInfo();
        // prey2.getInfo();

        /* predators */
        if (timeStep >= predIntro[0])
            pred1->measureDensity(world.landscapeTablePtr);

        // for (int i = 0; i < predatorTypesNb; i++)
        // {
        //     predators[i]->measureDensity(world.landscapeTablePtr);
        // }

        // pred1->getInfo();

        /* ---- save measures and snapshot in files ---- */
        if (timeStep % freqResu == 0)
            world.saveMeasures(resultsTableName, timeStep);

        // if (timeStep % freqSnap == 0)
        //     world.snapshot(snapshotTableName, timeStep);

        /* ---- check extinctions ---- */
        // if (prey1->currentPopulationSize == 0 | prey2->currentPopulationSize == 0 | pred1->currentPopulationSize == 0)
        // if (prey1->currentPopulationSize == 0 && prey1->initialDensity != 0 | prey2->currentPopulationSize == 0 && prey2->initialDensity != 0 | pred1->currentPopulationSize == 0 && pred1->initialDensity != 0)
        // {
        //     cout << "Time step #" << timeStep << ": at least one population got extinct, stop simulation" << endl
        //          << endl;
        //     // break;
        // }

        /* ---- reset landscape table ---- */
        world.resetCounts();

        if (timeStep % freqRfll == 0)
            world.resetResources();

        timeStep++;
    }

    /* ---- free memory ---- */
    cout << "simulation ended, free memory, call destructors" << endl
         << endl;

    /* diets table */
    int rmax = memberMatchingListsSize; // the exact dimension of the diet table
    for (int row = 0; row < rmax; row++)
    {
        delete[] dietsTable[row];
    }
    delete[] dietsTable;
    dietsTable = NULL;

    /* individuals matching lists */
    delete[] memberTypes;
    delete[] typeTags;
    delete[] indexInLandscape;
    memberTypes = NULL;
    typeTags = NULL;
    indexInLandscape = NULL;
}

/*
To compile and run:
> g++ -Wall -g chapter2modelV3.cpp -o outfile.o
> ./outfile.o
*/