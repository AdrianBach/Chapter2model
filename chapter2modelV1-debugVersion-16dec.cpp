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

/* The following function allocates memory to 3 arrays of the same size each containing
   the food web members' type, associated tag and assiciated column index in the
   landscape table.
   This way, you only need to know the row index of a member type to have access
   to all these informations */

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

    /* debug title and headers */
    cout << "assignTagsIndexes() :" << endl;
    cout << "member types - tag - index in landscape table" << endl;

    int r = 0; // initialise row counter
    while (r < memberMatchingListsSize)
    {
        for (int res = 0; res < resourceTypesNb; res++)
        {
            memberTypes[r] = resourceTypes[res];
            typeTags[r] = resourceTagStart + res;
            indexInLandscape[r] = 3 + res; // before: cellCode - x - y

            /* debug : OK */
            cout << memberTypes[r] << " - " << typeTags[r] << " - " << indexInLandscape[r] << endl;

            r++;
        }

        for (int prey = 0; prey < preyTypesNb; prey++)
        {
            memberTypes[r] = preyTypes[prey];
            typeTags[r] = preyTagStart + prey;
            indexInLandscape[r] = 3 + resourceTypesNb + prey; // before: x - y - resource densities

            /* debug : OK */
            cout << memberTypes[r] << " - " << typeTags[r] << " - " << indexInLandscape[r] << endl;

            r++;
        }

        for (int pred = 0; pred < predatorTypesNb; pred++)
        {
            memberTypes[r] = predatorTypes[pred];
            typeTags[r] = predatorTagStart + pred;
            indexInLandscape[r] = 3 + resourceTypesNb + 2 * preyTypesNb + pred; // before: x - y - resource densities - prey densities - prey catches

            /* debug : OK */
            cout << memberTypes[r] << " - " << typeTags[r] << " - " << indexInLandscape[r] << endl;

            r++;
        }
    }

    /* debug convinience */
    cout << endl;
    cout << endl;
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

    /* debug : OK */
    cout << "makeDietsTable() :" << endl;
    cout << "   ";
    for (int i = 0; i < memberMatchingListsSize; i++)
    {
        cout << memberTypes[i];
        if (i == (memberMatchingListsSize - 1)) // last loop
            cout << endl;
        else
            cout << " ";
    }
    for (int row = 0; row < dietsTableSize; row++)
    {
        cout << memberTypes[row] << " ";
        for (int col = 0; col < dietsTableSize; col++)
        {
            cout << dietsTable[row][col];
            if (col == (dietsTableSize - 1)) // last loop
                cout << endl;
            else
                cout << " ";
        }
    }
    cout << endl;
    cout << endl;
}

/* the following function gets a member type's diet from its index in the members'
   matching list and return the column index in the landscape table of each type
   present in the diet (useful to have access to their density in a particular cell) */

int *getDietLandscapeIndexes(int MembersMatchingListsIndex) // get the diet of a particular member of the food chain from its typeTag.
{

    int *dietArray = NULL; // declare the array to be returned

    int rowMax = memberMatchingListsSize; // row number of dietTable

    /* sum the values of the columns to get the number of members the individual feeds on */
    int typesInDiet = sumColumn(dietsTable, rowMax, MembersMatchingListsIndex);

    /* create int array of this size to store their tags */
    dietArray = new int[typesInDiet + 1]; // +1 because the first element will be the number of types in the diet

    dietArray[0] = typesInDiet;

    /* loop over lines to fill the array with their tags */
    int i = 1; // initialise counter to iterate through diet array. value of 1 to let the first element (0) contain the size of the array
    for (int row = 0; row < rowMax; row++)
    {
        if (dietsTable[row][MembersMatchingListsIndex] == 1 && i <= typesInDiet) // control not to add to the array more than its size
        {
            dietArray[i] = indexInLandscape[row]; // the row number of this match is enough to know its index in the landscape table!
            i++;                                  // increment i for next 1 value in the column
        }
    }

    return dietArray;
    /* WARNING: be sure to delete the array returned when not need anymore
       currently in animal.freeMemory() where it is called */

    dietArray = NULL;
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
    int **landscapeTablePtr;        // pointer to the landscape table
    int landscapeMatchingListsSize; //

public:                    // can be used / called in the main function
    string *XYcoordinates; // landscape matching lists
    int *cellCode;         //

    landscape(int size, int maxResource) // constructor: function that creates objects of the class by assigning values to or initializing the variables
    {
        /* debug */
        cout << "landscape constructor: assigning values to class variables..." << endl
             << endl;

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
        /* debug */
        cout << "creating landscape ... " << endl
             << endl;

        /* create landscape matching lists */
        landscapeMatchingListsSize = rowNb;
        XYcoordinates = new string[landscapeMatchingListsSize];
        cellCode = new int[landscapeMatchingListsSize];

        /* create a dynamic 2D array
        Super clear video: https://www.youtube.com/watch?v=mGl9LO-je3o&list=PL43pGnjiVwgSSRlwfahAuIqoJ8TfDIlHq&index=6&ab_channel=CodeBeauty */

        /* debug */
        cout << "allocating memory ... " << endl
             << endl;

        landscapeTablePtr = new int *[rowNb]; // define the pointer to an array of int pointers for each row

        for (int row = 0; row < rowNb; row++) // for each row, create a pointer to an array of size columnNb
        {
            landscapeTablePtr[row] = new int[columnNb];
        }

        /* fill the table with info */

        /* debug */
        cout << "initialising values ... (here dumb values just to make sure everything is where expected) " << endl
             << endl;

        /* debug title and headers */
        cout << "creating 2 corresponding arrays for x and y cell coordinates (string type) and cell code (int) :" << endl;
        cout << "concatenation of x and y - cellCode" << endl;

        /* x/y cell coordinates */
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

            /* debug : OK */
            cout << XYcoordinates[r] << " - " << cellCode[r] << endl;

            y++;
            r++;
        }

        /* debug convenience */
        cout << endl;

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
                // landscapeTablePtr[r][preyColumnStart + prey] = 0;

                /* a specific number to check if the info is where expected: OK */
                landscapeTablePtr[r][preyColumnStart + prey] = 1 + prey;

            r++;
        }

        /* initialise predator densities to 0 */

        r = 0;
        while (r < rowNb)
        {
            for (int pred = 0; pred < predatorTypesNb; pred++)
                // landscapeTablePtr[r][predatorColumnStart + pred] = 0;

                /* a specific number to check if the info is where expected: OK */
                landscapeTablePtr[r][predatorColumnStart + pred] = 11 + pred;

            r++;
        }

        /* debug */
        cout << "done." << endl
             << endl;
    }

    void freeMemory() // free memory allocated to landscape table
    {
        /* debug */
        cout << "deleting landscape table arrays and pointers..." << endl
             << endl;

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

    void fill() // function to reset resources to maximum
    {

        /* debug */
        cout << "refilling resources to maximum ... (here more to make sure it works)" << endl
             << endl;

        int r = 0;
        while (r < rowNb)
        {
            for (int res = 0; res < resourceTypesNb; res++)
                // landscapeTablePtr[r][resColumnStart + res] = MaxResource;

                /* a specific number to check if it works as expected */
                landscapeTablePtr[r][resColumnStart + res] = 10 * MaxResource;

            r++;
        }
    }

    void resetCatches() // function to reset catches to zero
    {

        /* debug */
        cout << "reseting catches to zero ... " << endl
             << endl;

        int r = 0;
        while (r < rowNb)
        {
            for (int prey = 0; prey < preyTypesNb; prey++)
                landscapeTablePtr[r][preyColumnStart + preyTypesNb + prey] = 0;

            r++;
        }
    }

    void getInfo() // function to cast out the landscape table and check if all good
    {

        cout << "landscape table" << endl;
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
        cout << endl;
    }

    void createResultsTable(string name) // creates a resultsTable
    {

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

        /* debug */
        cout << "results table csv created in the folder." << endl
             << endl;
    }

    void createSnapshotTable(string name) // creates a snapshot file
    {

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

        cout << "snapshot table csv created in the folder." << endl
             << endl;
    }

    void saveMeasures(string name, int ts) // write sum of results columns in the results file
    {

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

        cout << "results table csv appended." << endl
             << endl;
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

        cout << "snapshot table csv appended." << endl
             << endl;
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

    /* matching informations */
    int membersMatchingListsIndex; //
    int *dietLandscapeIndexes;     // array containing their column index in the landscape table

    /* individual level variables that will change during simulation */
    // int xCell, yCell;
    // int resourceConsumed; // an animal's resource pool
    // int deadOrAlive;      // 0 if dead 1 if alive

protected:
    /* population level variables */
    int currentPopulationSize; // current population size of the animal

public:
    animals(int InitialDensity, int TypeTag, int MaxMove, int ReproductionCost, int MaintenanceCost, int LandscapeSize) // initialise the constants shared by all animal types
    {

        /* debug */
        cout << "animal's constructor: assigning values to class variables..." << endl
             << endl;

        initialDensity = InitialDensity;
        typeTag = TypeTag;
        maxMove = MaxMove;
        reproductionCost = ReproductionCost;
        maintenanceCost = MaintenanceCost;
        landscapeSize = LandscapeSize;

        currentPopulationSize = initialDensity; // initialise current population size variable

        membersMatchingListsIndex = getMemberIndexFromTag(typeTag);
        dietLandscapeIndexes = getDietLandscapeIndexes(membersMatchingListsIndex);

        /* debug : OK */
        cout << "getting the animal's diet informations..." << endl;
        cout << typeTag << "'s row index in the member matching lists is " << membersMatchingListsIndex << endl;
        cout << "information about the type(s) present in the diet of type tag " << typeTag << " is(are) landscape table's column(s) number ";
        for (int i = 0; i < dietLandscapeIndexes[0]; i++)
        {
            cout << dietLandscapeIndexes[i + 1] << " "; // idem
        }
        cout << endl
             << endl;
    }

    int **create(string *XYcoordinates, int *CellCodes) // function to allocate memory to and fill the animal table
    {

        /* debug */
        cout << "creating animal's population table..." << endl
             << endl;

        columnNb = 5; // x - y - cellCode - resourceConsumed - deadOrAlive

        rowNb = initialDensity;

        /* debug */
        cout << "allocating memory to a temporary double pointer..." << endl
             << endl;

        /* create a 2D dynamic array */
        int **pointerToTable;

        pointerToTable = new int *[rowNb]; // define the pointer to an array of int pointers for each row

        for (int row = 0; row < rowNb; row++) // for each row, create a pointer to an integer array of size columnNb
        {
            pointerToTable[row] = new int[columnNb];
        }

        /* add the animals with their info */

        /* debug */
        cout << "initialising variables..." << endl
             << endl;

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

            r++;
        }

        /* debug */
        cout << "returning double pointer and deleting the temporary variable" << endl
             << endl;

        return pointerToTable;

        /* delete former table pointer without freeing the memory
           WARNING HERE: be sure that the previous line works fine otherwise the address of the table will be lost FOR EVER */
        pointerToTable = NULL;
    }

    void freeMemory(int **pointerToTable) // free memory allocated to animals table
    {

        /* debug */
        cout << "freeing animal's table array, tables and pointers..." << endl
             << endl;

        for (int row = 0; row < rowNb; row++) // free memory allocated to each row
        {
            delete[] pointerToTable[row];
        }

        delete[] pointerToTable; // free the memory allocated to the array of pointer to each row

        pointerToTable = NULL; // erase the address of the array of pointers to rows

        delete[] dietLandscapeIndexes;
        dietLandscapeIndexes = NULL;
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
        cout << endl;
    }

    /* next functions:
    - measureDensity(int **tablePtr) -> sum of animals in each cell to update landscapeTable
    - randomMove()
    - reproduce()
    - die()
    */

    // virtual void feed() = 0;
};

class prey : public animals // preys are one type of animal object, they share the same charasteristics but also have specificities
{

private:
    // float consumptionRate; // fraction of available resources collected by prey
    // max daily consumption ???

public:
    prey(int preyInitialDensity, int preyTypeTag, int preyMaxMovingDistance, int preyReproductionCost, int preyMaintenanceCost, int LandscapeSize) : animals(preyInitialDensity, preyTypeTag, preyMaxMovingDistance, preyReproductionCost, preyMaintenanceCost, LandscapeSize) //
    {
    }

    // void assignPreyVariables(float preyConsumptionRate)
    // {
    //     consumptionRate = preyConsumptionRate;
    // }

    // void feed(int **pointerToTable, int *Diet, int *DietLandscapeIndexes, int MemberMatchingListsIndex)
    // {

    //     int ind = 0; // initialise individual counter
    //     while (ind < currentPopulationSize)
    //     {
    //         int indCellCode = pointerToTable[ind][2]; // get individual's cellCode

    //         /* for each type in the diet match tag to the associated column index in landscape table
    //            note that the first element of Diet should be the number of different types in the diet */
    //         int dietSize = DietLandscapeIndexes[0];
    //         for (int i = 0; i < dietSize; i++)
    //         {
    //             /* get resource amount */
    //             int resourceIndex = DietLandscapeIndexes[i + 1]; // resource column index in landscape table
    //             int resourceAvailable = landscapeTablePtr[indCellCode][resourceIndex];

    //             /* compute consumption */
    //             if (resourceAvailable >= consumptionRate) // if there is more than consumption rate, consume max
    //             {
    //                 landscapeTablePtr[indCellCode][resourceIndex] -= consumptionRate; // subtract to landscape cell
    //                 pointerToTable[ind][4] += consumptionRate;                        // add resource consumed to the individual resource pool
    //             }
    //             else // else consume what is left (should also work if resourceAvailable=0)
    //             {
    //                 landscapeTablePtr[indCellCode][resourceIndex] -= resourceAvailable;
    //                 pointerToTable[ind][4] += resourceAvailable;
    //             }
    //         }

    //         ind++; // next individual
    //     }
    // }
};

class predator : public animals // predators are another type of animal object
{
public:
    predator(int initialDensity, int typeTag, int maxMovingDistance, int predatorReproductionCost, int predatorMaintenanceCost, int LandscapeSize) : animals(initialDensity, typeTag, maxMovingDistance, predatorReproductionCost, predatorMaintenanceCost, LandscapeSize) //
    {
    }

    /* feeding function */
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

    makeDietsTable(); // NEED 1 -r delete[] : OK

    worldSize = 3; // hard coded NOT GOOD

    srand(time(0)); // set random generator seed with instant time

    /* intiate world */

    /* contruct and create landscape */
    landscape world(worldSize, 10);
    world.create();

    /* check if everything is where expected: OK */
    world.getInfo();

    /* check fill function: OK
    world.fill();
    world.getInfo();
    */

    /* check resetCatches function: OK
    world.resetCatches();
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

    /* debug */
    cout << "prey 1 :" << endl;

    prey prey1(preysInitialDensities[0], 201, 10, 10, 10, worldSize);        // construct prey1 population = assigning values to the constants, intitialise some variables, compute others
    int **prey1TablePtr = prey1.create(world.XYcoordinates, world.cellCode); // allocate memory for prey1 population table
    prey1.getInfo(prey1TablePtr);

    /* debug */
    cout << "prey 2 :" << endl;
    prey prey2(preysInitialDensities[1], 202, 10, 10, 10, worldSize);
    int **prey2TablePtr = prey2.create(world.XYcoordinates, world.cellCode);
    prey1.getInfo(prey2TablePtr);

    /* debug */
    cout << "predator 1 :" << endl;
    predator predator1(predatorsInitialDensities[0], 301, 10, 10, 10, worldSize);
    int **predatorTablePtr = predator1.create(world.XYcoordinates, world.cellCode);
    predator1.getInfo(predatorTablePtr);

    /* create results and snapshot csv files */

    /* file names */
    string resultsTableName = simulationName + "-ResultsTable.csv";
    string snapshotTableName = simulationName + "-SnapshotTable.csv";

    world.createResultsTable(resultsTableName);
    world.createSnapshotTable(snapshotTableName);

    /* start simulation */

    timeMax = 3; // hard coded NOT GOOD

    int timeStep = 0; // time step variable

    /* debug */
    cout << "starting simulation over " << timeMax << " time steps" << endl
         << endl;

    while (timeStep < timeMax)
    {

        /* debug */
        cout << "time step " << timeStep << endl
             << endl;

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

    /* debug */
    cout << "deleting all arrays tables and pointers in reverse order ..." << endl
         << endl;

    /* animals tables : including animal's diet and dietLandscapeIndexes arrays */
    cout << "predator1" << endl;
    predator1.freeMemory(predatorTablePtr);
    cout << "prey2" << endl;
    prey2.freeMemory(prey2TablePtr);
    cout << "prey1" << endl;
    prey1.freeMemory(prey1TablePtr);

    /* landscape table : including landscape matching lists */
    world.freeMemory();

    /* diets table */
    cout << "diets table" << endl << endl;
    int rmax = memberMatchingListsSize; // the exact dimension of the diet table
    for (int row = 0; row < rmax; row++)
    {
        delete[] dietsTable[row];
    }
    delete[] dietsTable;
    dietsTable = NULL;

    /* members matching lists */
    cout << "members matching lists" << endl << endl;
    delete[] memberTypes;
    delete[] typeTags;
    delete[] indexInLandscape;
    memberTypes = NULL;
    typeTags = NULL;
    indexInLandscape = NULL;

    /* parameter arrays */
    cout << "parameter arrays" << endl << endl;
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

    /* debug */
    cout << "end of program" << endl;
}

/*
To compile :
> g++ -Wall -g chapter2modelV1-debugVersion-16dec.cpp -o outfile.o
> ./outfile.o
*/