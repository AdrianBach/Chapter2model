#include <iostream>  // input and output from the console
#include <string>    // manipulate char strings
#include <fstream>   // handle files
#include <time.h>    // get the time for random number generator
#include <stdlib.h>  // random generator tools
#include <algorithm> // random shuffle
#include <vector>    // needed for random shuffle

using namespace std; // not to have to write std:: in front of every call

/* ---------------------------- Global variables ---------------------------- */

/* give a name to the simulation */
// hard coded NOT GOOD
string simulationName = "test";

/* define animals and resources modelled */
// hard coded NOT GOOD
int resourceTypesNb = 2; // Number of resources available
int preyTypesNb = 2;     // Number of preys modelled
int predatorTypesNb = 1; // Number of predators modelled

/* pointers to individual types arrays */
string *resourceTypes; // resource1..n
string *preyTypes;     // prey1..n
string *predatorTypes; // types of predator simulated: predator1..n

/* pointers to initial population sizes arrays */
int *predatorsInitialDensities;
int *preysInitialDensities;

/* pointers to members' matching lists (arrays) + size variables */
int memberMatchingListsSize = resourceTypesNb + preyTypesNb + predatorTypesNb; // total number of types

string *memberTypes;   //
int *typeTags;         //
int *indexInLandscape; // column index of each type in the landscape table

/* pointer to diets' table: will indicate who eats who and who does not */
int **dietsTable;

/* world size */
int worldSize; // side of the squared lanscape in number of cells [0;+inf[

/* time variables */
int timeMax; // simulation time

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

double randomNumberGenerator(double min, double max) // generates a random number between min and max
{
    double res;
    res = min + (double)rand() * (max - min + 1) / (RAND_MAX - 1);

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

void makeDietsTable() // this table allows for each member of the system to eat and be eaten by another
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
    dietsTable[0][2] = 1; // prey1 feeds on resource1
    dietsTable[1][3] = 1; // prey2 feeds on resource2
    dietsTable[2][4] = 1; // predator1 feeds on prey1
    dietsTable[3][4] = 1; // predator1 feeds on prey2

    /* debug : OK
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
    */
}

vector<int> getDietLandscapeIndexes(int MembersMatchingListsIndex) // get the diet of a particular member of the food chain from its typeTag.
{

    vector<int> dietIndexes; // declare the vector to be returned

    int rowMax = memberMatchingListsSize; // row number of dietTable

    int typesInDiet = sumColumn(dietsTable, rowMax, MembersMatchingListsIndex);

    /* loop over lines to fill vector with diet's members column index in landscape table */
    for (int row = 0; row < rowMax; row++)
    {
        if (dietsTable[row][MembersMatchingListsIndex] == 1 && dietIndexes.size() <= typesInDiet) // control not to add to the array more than its size
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

// void doublePtrFreeMemory(int** pointerToTable){}
// Good idea but would need access to row and column number that are class protected variables

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
    int MaxResource;         // max amount of resources on a cell
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

    landscape(int size, int maxResource) // constructor: function that creates objects of the class by assigning values to or initializing the variables
    {
        Size = size;
        MaxResource = maxResource;

        /* column index where every group of info starts in the table*/
        resColumnStart = 3;                                      // before: cellCode - x - y
        preyColumnStart = resColumnStart + resourceTypesNb;      // before: x - y - resource densities
        predatorColumnStart = preyColumnStart + 2 * preyTypesNb; // before: x - y - resource densities - prey densities - prey catches

        /* table size */
        rowNb = Size * Size;
        columnNb = resColumnStart + resourceTypesNb + 2 * preyTypesNb + predatorTypesNb;
    }

    void create() // function to allocate memory to and fill the landscape table
    {

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
                landscapeTablePtr[r][resColumnStart + res] = MaxResource;

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

    void freeMemory() // free memory allocated to landscape table
    {

        /* free matching list memory */
        delete[] XYcoordinates;
        XYcoordinates = NULL;
        delete[] cellCode;
        cellCode = NULL;

        /* free landscape tabel memory */
        for (int row = 0; row < rowNb; row++) // free memory allocated to each row
        {
            delete[] landscapeTablePtr[row];
        }

        delete[] landscapeTablePtr; // free the memory allocated to the array of pointer to each row

        landscapeTablePtr = NULL; // erase the address of the array of pointers to rows
    }

    void resetLandscape() // function to reset resources to maximum and counts to 0
    {

        /* debug */
        cout << "reinitialising landscape ..." << endl
             << endl;

        int r = 0;
        while (r < rowNb)
        {

            /* refill resources */
            for (int res = 0; res < resourceTypesNb; res++)
                landscapeTablePtr[r][resColumnStart + res] = MaxResource;

            /* resetting counts to O */
            for (int col = preyColumnStart; col < columnNb; col++)
                landscapeTablePtr[r][col] = 0;

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
                snapshotTable << preyTypes[i] << "PopulationSize"
                              << ",";
            for (int i = 0; i < preyTypesNb; i++)
                snapshotTable << preyTypes[i] << "catches"
                              << ",";
            for (int i = 0; i < predatorTypesNb; i++)
            {
                snapshotTable << predatorTypes[i] << "PopulationSize";
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
        /* debug */
        cout << "appending " << name << endl
             << endl;

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
        /* debug */
        cout << "appending " << name << endl
             << endl;

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

                snapshotTable << "\n";

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
    int rowNb;         // total number of animals of all types
    int columnNb;      // number of info stored in the table
    int landscapeSize; //

    /* population level constants that might have different values according to animal types */
    int initialDensity;
    int typeTag;          // a integer tag for each animal type
    int maxMove;          // max number of cells an animal can move by in each direction
    int reproductionCost; // resources needed to reproduce
    int maintenanceCost;  // resources needed to survive a time step
    int maxOffspring;     // max number of offspring when passing reproduction trial

protected:
    /* population level variables */
    int currentPopulationSize; // current population size of the animal

    /* matching informations */
    int membersMatchingListsIndex;    //
    vector<int> dietLandscapeIndexes; // array containing their column index in the landscape table

public:
    animals(int InitialDensity, int TypeTag, int MaxMove, int ReproductionCost, int MaintenanceCost, int LandscapeSize, int MaxOffspring) // initialise the constants shared by all animal types
    {

        initialDensity = InitialDensity;
        typeTag = TypeTag;
        maxMove = MaxMove;
        reproductionCost = ReproductionCost;
        maintenanceCost = MaintenanceCost;
        maxOffspring = MaxOffspring;
        landscapeSize = LandscapeSize;

        if (maxMove >= landscapeSize)
            cout << "Warning animal's moving range is larger than the landscape size, can mess cell coordinates when moving" << endl;

        currentPopulationSize = initialDensity; // initialise current population size variable

        membersMatchingListsIndex = getMemberIndexFromTag(typeTag);
        dietLandscapeIndexes = getDietLandscapeIndexes(membersMatchingListsIndex);
    }

    int **create(string *XYcoordinates, int *CellCodes) // function to allocate memory to and fill the animal table
    {

        columnNb = 6; // x - y - cellCode - resourceConsumed - deadOrAlive - offspring

        rowNb = initialDensity;

        /* create a 2D dynamic array */
        int **pointerToTable;

        pointerToTable = new int *[rowNb]; // define the pointer to an array of int pointers for each row

        for (int row = 0; row < rowNb; row++) // for each row, create a pointer to an integer array of size columnNb
        {
            pointerToTable[row] = new int[columnNb];
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

            /* assign random coordinates between 0 and Size*/
            pointerToTable[r][0] = int(randomNumberGenerator(0, landscapeSize - 1));
            pointerToTable[r][1] = int(randomNumberGenerator(0, landscapeSize - 1));

            /* update pointerToTable with prey info - hard coded NOT GOOD */
            pointerToTable[r][2] = getCellCode(XYcoordinates, CellCodes, landscapeSize, pointerToTable[r][0], pointerToTable[r][1]);
            pointerToTable[r][3] = 0; // initialise resource consumed
            pointerToTable[r][4] = 1; // the animal is alive
            pointerToTable[r][5] = 0; // the animal has no offspring yet

            r++;
        }

        return pointerToTable;
    }

    void freeMemory(int **pointerToTable) // free memory allocated to animals table
    {

        for (int row = 0; row < rowNb; row++) // free memory allocated to each row
        {
            delete[] pointerToTable[row];
        }

        delete[] pointerToTable; // free the memory allocated to the array of pointer to each row

        pointerToTable = NULL; // erase the address of the array of pointers to rows
    }

    void randomMove(int **pointerToTable, string *XYcoordinates, int *CellCodes)
    {
        /* debug */
        cout << memberTypes[membersMatchingListsIndex] << " are moving at random" << endl
             << endl;

        /* iterate through individuals */
        int ind = 0;

        while (ind < currentPopulationSize)
        {
            /* update position with random number in moving range */
            pointerToTable[ind][0] += randomNumberGenerator(-1 * maxMove, maxMove); // xCell
            pointerToTable[ind][1] += randomNumberGenerator(-1 * maxMove, maxMove); // yCell

            /* check boundaries */
            pointerToTable[ind][0] = transmissiveBoundaries(pointerToTable[ind][0], landscapeSize); // xCell
            pointerToTable[ind][1] = transmissiveBoundaries(pointerToTable[ind][1], landscapeSize); // yCell

            pointerToTable[ind][2] = getCellCode(XYcoordinates, CellCodes, landscapeSize, pointerToTable[ind][0], pointerToTable[ind][1]); // update cell code

            ind++; // next individual
        }
    }

    void survivalTrial(int **pointerToTable)
    {
        /* debug */
        cout << memberTypes[membersMatchingListsIndex] << " survival trials" << endl
             << endl;

        /* iterate through individuals */
        int ind = 0;

        while (ind < currentPopulationSize)
        {

            /* if resonsume < maintenance cost -> change deadOrAlive to 0 else remove maintenance cost from the individual's resource pool */
            if (pointerToTable[ind][3] < maintenanceCost)
                pointerToTable[ind][4] = 0;
            else
                pointerToTable[ind][3] -= maintenanceCost;

            ind++; // next individual
        }
    }

    void reproductionTrial(int **pointerToTable)
    {
        /* debug */
        cout << memberTypes[membersMatchingListsIndex] << " reproduction trials" << endl
             << endl;

        /* iterate through individuals */
        int ind = 0;

        while (ind < currentPopulationSize)
        {

            /* if resConsumed >= reproduction cost -> change offspring and substract reproduction cost from the individual's resource pool */
            if (pointerToTable[ind][3] >= reproductionCost)
            {
                pointerToTable[ind][5] = randomNumberGenerator(0, maxOffspring); // even if it has the resource it might not reproduce (and thus not paying the cost)
                if (pointerToTable[ind][5] > 0)
                    pointerToTable[ind][3] -= reproductionCost;
            }

            ind++; // next individual
        }
    }

    void measureDensity(int **pointerToTable, int **LandscapeTable)
    {

        /* debug */
        cout << "measuring " << memberTypes[membersMatchingListsIndex] << " densities" << endl
             << endl;

        /* get animal's column index in landscapeTable */
        int colIndex = indexInLandscape[membersMatchingListsIndex];

        /* iterate through population table */
        int ind = 0;
        while (ind < currentPopulationSize)
        {
            /* increment landscape table at the right position */
            int rowIndex = pointerToTable[ind][2];
            LandscapeTable[rowIndex][colIndex] += 1;

            /* debug
            cout << "incrementing cell [" << rowIndex << "] [" << colIndex << "]" << endl
                 << endl;
            */

            ind++; // next individual
        }
    }

    void getInfo(int **pointerToTable) // function to cast out the animalsTable and check if all good
    {

        /* cast out the column names */
        cout << typeTag << "'s population table" << endl;
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
             << endl;

        /* iterate through lines and column to cast out the values */
        int r = 0;
        while (r < rowNb)
        {
            for (int column = 0; column < columnNb; column++)
                cout << pointerToTable[r][column] << " ";
            cout << endl;

            r++;
        }
        cout << endl;
    }

    /* next functions:
    - updateTable(int **pointerToTable) // has to reset offspring to zero as well
    */
};

class prey : public animals // preys are one type of animal object, they share the same charasteristics but also have specificities
{

private:
    int maxConsumption; // max resource units the prey can consume in a time step

public:
    prey(int preyInitialDensity, int preyTypeTag, int preyMaxMovingDistance, int preyReproductionCost, int preyMaintenanceCost, int LandscapeSize, int preyMaxOffspring) : animals(preyInitialDensity, preyTypeTag, preyMaxMovingDistance, preyReproductionCost, preyMaintenanceCost, LandscapeSize, preyMaxOffspring) //
    {
    }

    void assignPreyVariables(int preyMaxConsumption)
    {
        maxConsumption = preyMaxConsumption;
    }

    void feed(int **pointerToTable, int **landscapeTable)
    {

        /* debug */
        cout << memberTypes[membersMatchingListsIndex] << " are feeding" << endl
             << endl;

        /* shuffle the order of the individuals */
        vector<int> shuffledPop = shuffleOrder(currentPopulationSize);

        int ind = 0;                        // initialise individual counter
        while (ind < currentPopulationSize) // iterate through individuals
        {

            int rowIndex = shuffledPop[ind]; // get shuffled row index

            int indCellCode = pointerToTable[rowIndex][2]; // get individual's cellCode

            /* iteration through diet and feed */
            int dietSize = dietLandscapeIndexes.size();

            for (int i = 0; i < dietSize; i++)
            {
                /* get resource amount */
                int resourceColIndex = dietLandscapeIndexes[i]; // resource column index in landscape table

                int resourceAvailable = landscapeTable[indCellCode][resourceColIndex];

                /* compute maxConsumption */
                if (resourceAvailable >= maxConsumption) // if there is more than maxConsumption rate, consume max
                {
                    landscapeTable[indCellCode][resourceColIndex] -= maxConsumption; // subtract to landscape cell
                    pointerToTable[rowIndex][3] += maxConsumption;
                }
                else // else consume what is left (should also work if resourceAvailable=0)
                {
                    landscapeTable[indCellCode][resourceColIndex] -= resourceAvailable;
                    pointerToTable[rowIndex][3] += resourceAvailable;
                }
            }

            ind++; // next individual
        }
    }
};

class predator : public animals // predators are another type of animal object
{

protected:
    int maxCatches;     // number of preys a predator can catch a day
    int conversionRate; // equivalent of a catch in resources

public:
    predator(int initialDensity, int typeTag, int maxMovingDistance, int predatorReproductionCost, int predatorMaintenanceCost, int LandscapeSize, int predatorMaxOffspring) : animals(initialDensity, typeTag, maxMovingDistance, predatorReproductionCost, predatorMaintenanceCost, LandscapeSize, predatorMaxOffspring) //
    {
    }

    void assignPredatorVariables(int predatorMaxCatches, int predatorConversionRate)
    {
        maxCatches = predatorMaxCatches;
        conversionRate = predatorConversionRate;
    }

    void hunt(int **pointerToTable, int **landscapeTable)
    {

        /* debug */
        cout << memberTypes[membersMatchingListsIndex] << " are hunting" << endl
             << endl;

        /* shuffle the order of the individuals */
        vector<int> shuffledPop = shuffleOrder(currentPopulationSize);

        int ind = 0;                        // initialise individual counter
        while (ind < currentPopulationSize) // iterate through individuals
        {

            int rowIndex = shuffledPop[ind]; // get shuffled row index

            int indCellCode = pointerToTable[rowIndex][2]; // get individual's cellCode

            vector<int> shuffledDiet = dietLandscapeIndexes; // store diet landscape index before shuffling it

            random_shuffle(shuffledDiet.begin(), shuffledDiet.end()); // shuffle

            int catches = 0; // initialise a catch counter

            /* iterate through prey columns and while catches < maxCatches and
               that there are prey available, catch them */

            for (int i = 0; i < shuffledDiet.size(); i++)
            {
                /* get density */
                int dens = landscapeTable[indCellCode][shuffledDiet[i]];

                if (dens > 0 && catches < maxCatches) // if prey present and maxCatches is not attained yet
                {
                    landscapeTable[indCellCode][shuffledDiet[i] + preyTypesNb] += 1; // increment corresponding catch cell in landscape table
                    pointerToTable[rowIndex][3] += 1 * conversionRate;               // increment predator resource pool

                    catches++; // increment catches counter

                    /* debug */
                    cout << "a " << memberTypes[membersMatchingListsIndex] << " caught a prey on cell " << indCellCode << endl
                         << endl;
                }
            }

            ind++; // next individual
        }
    }
};

/* ------------------------------ Main program ------------------------------ */

int main()
{

    /* enter different types and initital densities */

    resourceTypes = new string[resourceTypesNb];
    resourceTypes[0] = "resource1"; // hard coded NOT GOOD
    resourceTypes[1] = "resource2"; // hard coded NOT GOOD

    preyTypes = new string[preyTypesNb];
    preyTypes[0] = "prey1"; // hard coded NOT GOOD
    preyTypes[1] = "prey2"; // hard coded NOT GOOD

    preysInitialDensities = new int[preyTypesNb];
    preysInitialDensities[0] = 10; // hard coded NOT GOOD
    preysInitialDensities[1] = 15; // hard coded NOT GOOD

    predatorTypes = new string[predatorTypesNb]; // assign memory to the pointer
    predatorTypes[0] = "predator1";              // hard coded NOT GOOD

    predatorsInitialDensities = new int[predatorTypesNb];
    predatorsInitialDensities[0] = 5; // hard coded NOT GOOD

    assignTagsIndexes(); // creates individuals' matching tables. NEED 3 delete[] : OK

    makeDietsTable(); // NEED 1 delete[] : OK

    worldSize = 3; // hard coded NOT GOOD

    srand(time(0)); // set random generator seed with instant time

    /* intiate world */

    /* contruct and create landscape */
    landscape world(worldSize, 10);
    world.create();

    /* check if everything is where expected: OK
    world.getInfo();
    */

    /* construct and create animals table (iteratively if possible) : OK but not iterative yet
       should work like this:
       > prey preyN() // class constructor
       > int **preyNTablePtr = preyN.create(world.XYcoordinates, world.cellCode); // allocate memory to temporary animalsTablePointer and initialise values
       .
       .
       .
       > preyN.freeMemory(prey1TablePtr)

       CAREFUL : free memory each time using animals.create()

       CAN BE OPTIMISED: https://www.youtube.com/watch?v=T8f4ajtFU9g&list=PL43pGnjiVwgTJg7uz8KUGdXRdGKE0W_jN&index=6&ab_channel=CodeBeauty
    */

    prey prey1(preysInitialDensities[0], 201, 1, 4, 5, worldSize, 1); // construct prey1 population = assigning values to the constants, intitialise some variables, compute others
    prey1.assignPreyVariables(5);
    int **prey1TablePtr = prey1.create(world.XYcoordinates, world.cellCode); // allocate memory for prey1 population table

    prey prey2(preysInitialDensities[1], 202, 1, 6, 5, worldSize, 1);
    prey2.assignPreyVariables(7);
    int **prey2TablePtr = prey2.create(world.XYcoordinates, world.cellCode);

    predator pred1(predatorsInitialDensities[0], 301, 1, 5, 1, worldSize, 1);
    pred1.assignPredatorVariables(1, 1);
    int **pred1TablePtr = pred1.create(world.XYcoordinates, world.cellCode);

    /* check create function : OK
    prey1.getInfo(prey1TablePtr);
    prey1.getInfo(prey2TablePtr);
    pred1.getInfo(pred1TablePtr);
    */

    /* create results and snapshot csv files */

    /* file names */
    string resultsTableName = simulationName + "-ResultsTable.csv";
    string snapshotTableName = simulationName + "-SnapshotTable.csv";

    world.createResultsTable(resultsTableName);
    world.createSnapshotTable(snapshotTableName);

    /* start simulation */

    timeMax = 3; // hard coded NOT GOOD

    /* debug */
    cout << "starting simulation over " << timeMax << " time steps" << endl
         << endl;

    int timeStep = 0; // time step variable

    while (timeStep < timeMax)
    {

        /* debug */
        cout << "time step " << timeStep << endl
             << endl;

        /* life events
           CAN BE OPTIMIZED with iteration using pointers to objects : see above */

        /* moving */
        prey1.randomMove(prey1TablePtr, world.XYcoordinates, world.cellCode);
        prey2.randomMove(prey2TablePtr, world.XYcoordinates, world.cellCode);
        pred1.randomMove(pred1TablePtr, world.XYcoordinates, world.cellCode);

        /* feeding */
        prey1.feed(prey1TablePtr, world.landscapeTablePtr);
        prey2.feed(prey2TablePtr, world.landscapeTablePtr);
        pred1.hunt(pred1TablePtr, world.landscapeTablePtr);

        /* reproducing */
        prey1.reproductionTrial(prey1TablePtr);
        prey2.reproductionTrial(prey2TablePtr);
        pred1.reproductionTrial(pred1TablePtr);

        /* surviving */
        prey1.survivalTrial(prey1TablePtr);
        prey2.survivalTrial(prey2TablePtr);
        pred1.survivalTrial(pred1TablePtr);

        /* update animals dynamic arrays with birth catches and deaths */

        /* measure densities */
        prey1.measureDensity(prey1TablePtr, world.landscapeTablePtr);
        prey2.measureDensity(prey2TablePtr, world.landscapeTablePtr);
        pred1.measureDensity(pred1TablePtr, world.landscapeTablePtr);

        /* save measures and snapshot in files : OK */
        world.saveMeasures(resultsTableName, timeStep);
        world.snapshot(snapshotTableName, timeStep);

        /* reset landscape table */
        world.resetLandscape();

        timeStep++;
    }

    /* free memory */

    /* animals tables : including animal's diet and dietLandscapeIndexes arrays */
    prey1.freeMemory(prey1TablePtr);
    prey2.freeMemory(prey2TablePtr);
    pred1.freeMemory(pred1TablePtr);

    /* landscape table : including landscape matching lists */
    world.freeMemory();

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

    /* parameter arrays */
    delete[] predatorTypes;
    delete[] predatorsInitialDensities;
    delete[] preyTypes;
    delete[] preysInitialDensities;
    delete[] resourceTypes;
    predatorTypes = NULL;
    predatorsInitialDensities = NULL;
    preyTypes = NULL;
    preysInitialDensities = NULL;
    resourceTypes = NULL;
}

/*
To compile :
> g++ -Wall -g chapter2modelV1.cpp -o outfile.o
> ./outfile.o
*/