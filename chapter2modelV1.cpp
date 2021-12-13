#include <iostream> // input and output from the console
#include <string>   // manipulate char strings
#include <fstream>  // handle files
#include <time.h>   // get the time for random number generator
#include <stdlib.h> // random generator tools

using namespace std; // not to have to write std:: in front of every call

/* ---------------------------- Global variables ---------------------------- */

/* give a name to the simulation */
// hard coded NOT GOOD
string simulationName = "test";

/* define animals and resources modelled */
// hard coded NOT GOOD
int predatorTypesNb = 1; // Number of predators modelled
int preyTypesNb = 2;     // Number of preys modelled
int resourceNb = 2;      // Number of resources available

string *predatorTypes;        // types of predator simulated: predator1..n prey1..n
string *preyTypes;            // prey1..n
string *primaryResourceTypes; // resource1..n

/* initial population sizes */
int *predatorsInitialDensities;
int *preysInitialDensities;

/* world size */
int worldSize; // side of the squared lanscape in number of cells [0;+inf[

/* time variables */
int timeMax;

/* ---------------------- Some useful global functions ---------------------- */

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

/* MATCHING CELLS FUNCTION */

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

protected:                   // variables that should not be modified directly, nor accessed from the main function, but accessible to the other classes
    int **landscapeTablePtr; // pointer to the landscape table

public:                                  // functions that can modify the private and protected variables and can be called in the main function
    landscape(int size, int maxResource) // constructor: function that creates objects of the class by assigning values to or initializing the variables
    {
        Size = size;
        MaxResource = maxResource;
        rowNb = Size * Size;
        columnNb = 2 + resourceNb + 2 * preyTypesNb + predatorTypesNb;

        /* column index where every group of info starts in the table*/
        resColumnStart = 2;                                     // before: x - y
        preyColumnStart = 2 + resourceNb;                       // before: x - y - resource densities
        predatorColumnStart = 2 + resourceNb + 2 * preyTypesNb; // before: x - y - resource densities - prey densities - prey catches
    }

    void create() // function to allocate memory to and fill the landscape table
    {

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

            landscapeTablePtr[r][0] = x;
            landscapeTablePtr[r][1] = y;

            y++;
            r++;
        }

        /* initialise resources to maximum */

        r = 0;
        while (r < rowNb)
        {

            for (int res = 0; res < resourceNb; res++)
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

        for (int row = 0; row < rowNb; row++) // free memory allocated to each row
        {
            delete[] landscapeTablePtr[row];
        }

        delete[] landscapeTablePtr; // free the memory allocated to the array of pointer to each row

        landscapeTablePtr = NULL; // erase the address of the array of pointers to rows
    }

    void fill() // function to reset resources to maximum
    {
        int r = 0;
        while (r < rowNb)
        {
            for (int res = 0; res < resourceNb; res++)
                // landscapeTablePtr[r][2 + res] = MaxResource;

                /* a specific number to check if it works as expected */
                landscapeTablePtr[r][2 + res] = 10 * MaxResource;

            r++;
        }
    }

    void resetCatches() // function to reset catches to zero
    {
        int r = 0;
        while (r < rowNb)
        {
            for (int prey = 0; prey < preyTypesNb; prey++)
                landscapeTablePtr[r][2 + resourceNb + preyTypesNb + prey] = 0;

            r++;
        }
    }

    void getInfo() // function to cast out the landscape table and check if all good
    {
        /* cast out the column names */
        cout << "xCell"
             << " "
             << "yCell"
             << " ";
        for (int i = 0; i < resourceNb; i++)
            cout << primaryResourceTypes[i] << " ";
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
    }

    void createResultsTable(string name) // creates a resultsTable
    {

        /* write headers */
        resultsTable.open(name, ios::out);
        if (resultsTable.is_open())
        {
            resultsTable << "timeStep"
                         << ",";
            for (int i = 0; i < resourceNb; i++)
                resultsTable << primaryResourceTypes[i] << "amount"
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

        /* write headers */
        snapshotTable.open(name, ios::out);
        if (snapshotTable.is_open())
        {
            snapshotTable << "timeStep"
                          << ","
                          << "xCell"
                          << ","
                          << "yCell"
                          << ",";
            for (int i = 0; i < resourceNb; i++)
                snapshotTable << primaryResourceTypes[i] << "amount"
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

        resultsTable.open(name, ios::app);
        if (resultsTable.is_open())
        {
            resultsTable << ts
                         << ",";

            /* write the sum of the measure columns */

            for (int i = 0; i < resourceNb; i++)
                resultsTable << sumColumn(landscapeTablePtr, rowNb, 2 + i)
                             << ",";
            for (int i = 0; i < preyTypesNb; i++)
                resultsTable << sumColumn(landscapeTablePtr, rowNb, 2 + resourceNb + i)
                             << ",";
            for (int i = 0; i < preyTypesNb; i++)
                resultsTable << sumColumn(landscapeTablePtr, rowNb, 2 + resourceNb + preyTypesNb + i)
                             << ",";
            for (int i = 0; i < predatorTypesNb; i++)
            {
                resultsTable << sumColumn(landscapeTablePtr, rowNb, 2 + resourceNb + 2 * preyTypesNb + i);
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
    int rowNb;    // total number of animals of all types
    int columnNb; // number of info stored in the table

    /* these constants might have different values according to animal types */
    int initialDensity;
    int typeTag; // a integer tag for each animal type
    // int maxMove;          // max number of cells an animal can move by in each direction
    // int reproductionCost; // resources needed to reproduce
    // int maintenanceCost;  // resources needed to survive a time step

    /* these are individual level-variables that will change during simulation */
    // int xCell, yCell;
    // int resourceConsumed; // an animal's resource pool
    // int deadOrAlive;      // 0 if dead 1 if alive

public:
    animals(int InitialDensity, int TypeTag) // , int MaxMove, int ReproductionCost, int MaintenanceCost / initialise the constants shared by all animal types
    {

        initialDensity = InitialDensity;
        typeTag = TypeTag;
        // maxMove = MaxMove;
        // reproductionCost = ReproductionCost;
        // maintenanceCost = MaintenanceCost;
    }

    int **create(int landscapeSize) // function to allocate memory to and fill the animal table
    {
        columnNb = 5; // x - y - animalType - resourceConsumed - deadOrAlive

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
            pointerToTable[r][0] = int(randomNumberGenerator(0, landscapeSize-1));
            pointerToTable[r][1] = int(randomNumberGenerator(0, landscapeSize-1));

            /* update pointerToTable with prey info - hard coded NOT GOOD */
            pointerToTable[r][2] = typeTag;
            pointerToTable[r][3] = 0; // initialise resource consumed
            pointerToTable[r][4] = 1; // the animal is alive

            r++;
        }

        return pointerToTable;

        /* delete former table pointer without freeing the memory
           WARNING HERE: be sure that the previous line works fine otherwise the address of the table will be lost FOR EVER */
        pointerToTable = NULL;
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

    void getInfo(int **pointerToTable) // function to cast out the animalsTable and check if all good
    {

        /* cast out the column names */
        cout << "xCell"
             << " "
             << "yCell"
             << " "
             << "animalTag"
             << " "
             << "resourcesConsumed"
             << " "
             << "deadOrAlive"
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
    }

    /* next functions:
    - measureDensity(int **tablePtr) -> sum of animals in each cell to update landscapeTable
    - randomMove(int worldSize)
    - reproduce()
    - die()
    */

    // virtual void feed() = 0;
};

class prey : public animals // preys are one type of animal object, they share the same charasteristics but also have specificities
{

public:
    prey(int initialDensity, int typeTag) : animals(initialDensity, typeTag) // , int maxMovingDistance, int preyReproductionCost, int preyMaintenanceCost / , maxMovingDistance, preyReproductionCost, preyMaintenanceCost
    {
    }

    /* PROBLEM HERE complicated to add objects information
       - need to know where the previous object stopped to continue assigning
       - should every object/population have its own table ? no need to create a pop object ?
       - is it going to a problem to match cells ?
       - where does the matching function should be ? */

    /* Preys feeding function */
};

/* ------------------------------ Main program ------------------------------ */

int main()
{

    /* enter different types and initital densities */

    predatorTypes = new string[predatorTypesNb]; // assign memory to the pointer
    predatorTypes[0] = "predator1";              // hard coded NOT GOOD

    predatorsInitialDensities = new int[predatorTypesNb];
    predatorsInitialDensities[0] = 5; // hard coded NOT GOOD

    preyTypes = new string[preyTypesNb];
    preyTypes[0] = "prey1"; // hard coded NOT GOOD
    preyTypes[1] = "prey2"; // hard coded NOT GOOD

    preysInitialDensities = new int[preyTypesNb];
    preysInitialDensities[0] = 10; // hard coded NOT GOOD
    preysInitialDensities[1] = 15; // hard coded NOT GOOD

    primaryResourceTypes = new string[resourceNb];
    primaryResourceTypes[0] = "resource1"; // hard coded NOT GOOD
    primaryResourceTypes[1] = "resource2"; // hard coded NOT GOOD

    worldSize = 3; // hard coded NOT GOOD

    srand(time(0)); // set random generator seed with instant time

    /* intiate world */

    /* contruct and create landscape */
    landscape world(worldSize, 10);
    world.create();

    /* check if everything is where expected: OK
    world.getInfo();
    */

    /* check fill function: OK
    world.fill();
    world.getInfo();
    */

    /* check resetCatches function: OK
    world.resetCatches();
    world.getInfo();
    */

    /* convert animal types in integer tags
       matching strings : https://www.youtube.com/watch?v=aEgG4pidcKU&t=1s&ab_channel=CodeBeauty
    */

    /* edit interaction/diet table */

    /* construct and create animals table (iteratively if possible)
       should work like that: OK
       > prey prey1() // class constructor
       > int **prey1TablePtr = prey1.create(worldSize); // allocate memory to temporary animalsTablePointer and initialise values
       .
       .
       .
       > prey1.freeMemory(prey1TablePtr)

       CAREFUL : free memory each time using animals.create()

       CAN BE OPTIMISED: https://www.youtube.com/watch?v=T8f4ajtFU9g&list=PL43pGnjiVwgTJg7uz8KUGdXRdGKE0W_jN&index=6&ab_channel=CodeBeauty 
    */

    prey prey1(preysInitialDensities[0], 201);
    int **prey1TablePtr = prey1.create(worldSize);
    
    prey prey2(preysInitialDensities[1], 202);
    int **prey2TablePtr = prey2.create(worldSize);

    /* check create function : OK
    prey1.getInfo(prey1TablePtr);
    prey1.getInfo(prey2TablePtr);
    */

    /* create results and snapshot csv files */

    /* file names */
    string resultsTableName = simulationName + "-ResultsTable.csv";
    string snapshotTableName = simulationName + "-SnapshotTable.csv";

    world.createResultsTable(resultsTableName);
    world.createSnapshotTable(snapshotTableName);

    /* start simulation */

    timeMax = 3; // hard coded NOT GOOD

    int timeStep = 0; // time step variable

    while (timeStep < timeMax)
    {
        /* life events */

        /* take measures and snapshot: OK */
        world.saveMeasures(resultsTableName, timeStep);
        world.snapshot(snapshotTableName, timeStep);

        /* refill resources */
        world.fill();

        /* reset catches counts */
        world.resetCatches();

        timeStep++;
    }

    /* free memory */

    /* landscape table */
    world.freeMemory();

    /* animals tables */
    prey1.freeMemory(prey1TablePtr);
    prey2.freeMemory(prey2TablePtr);

    /* parameter arrays */
    delete[] predatorTypes;
    delete[] predatorsInitialDensities;
    delete[] preyTypes;
    delete[] preysInitialDensities;
    delete[] primaryResourceTypes;
}

/*
To compile :
> g++ -Wall -g chapter2modelV1.cpp -o outfile.o
> ./outfile.o
*/