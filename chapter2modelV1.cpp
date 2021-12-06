#include <iostream>
#include <string>
using namespace std;
// #include <list>

/* declare global variables here */

/* define animals and resources modelled */
int predatorNb = 1;           // Number of predators modelled
int preyNb = 2;               // Number of preys modelled
int resourceNb = 2;           // Number of resources available
string *predatorTypes;        // types of predator simulated: predator1..n prey1..n
string *preyTypes;            // prey1..n
string *primaryResourceTypes; // resource1..n

/* initial population sizes */
int *predatorsInitialDensities;
int *preysInitialDensities;

/*
create a class for landscape:
- all the variables that define it
- all the functions that manipulate these variables
- make sure that the values can be used by other classes
*/

/* create class */
class landscape
{
    /* list of landscape-specific variables */
private:             // variables that should not be modified directly, nor accessed from the main function
    int Size;        // side of the squared lanscape in number of cells [0;+inf[
    int MaxResource; // max amount of resources on a cell
    int rowNb;       // row number
    int columnNb;    // column number

protected:                   // variables that should not be modified directly, nor accessed from the main function, but accessible to the other classes
    int **landscapeTablePtr; // pointer to the landscape table

public:                                  // functions that can modify the private and protected variables and can be called in the main function
    landscape(int size, int maxResource) // constructor: function that creates objects of the class by assigning values to or initializing the variables
    {
        Size = size;
        MaxResource = maxResource;
        rowNb = Size * Size;
        columnNb = 2 + resourceNb + 2 * preyNb + predatorNb;
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
                landscapeTablePtr[r][2 + res] = MaxResource;

            r++;
        }

        /* initialise prey densities and catches to 0 */
        r = 0;
        while (r < rowNb)
        {
            for (int prey = 0; prey < (2 * preyNb); prey++)
                // landscapeTablePtr[r][2 + resourceNb + prey] = 0;

                /* a specific number to check if the info is where expected */
                landscapeTablePtr[r][2 + resourceNb + prey] = 1+prey;

            r++;
        }

        /* initialise predator densities to 0 */
        r = 0;
        while (r < rowNb)
        {
            for (int pred = 0; pred < predatorNb; pred++)
                // landscapeTablePtr[r][2+ resourceNb + 2*preyNb + pred] = 0;
                
                /* a specific number to check if the info is where expected */
                landscapeTablePtr[r][2+ resourceNb + 2*preyNb + pred] = 11+pred;

            r++;
        }
    }

    void fill() // function to reset resources to maximum
    {
        int r = 0;
        while (r < rowNb)
        {
            for (int res = 0; res < resourceNb; res++)
                // landscapeTablePtr[r][2 + res] = MaxResource;
                
                /* a specific number to check if it works as expected */
                landscapeTablePtr[r][2 + res] = 10*MaxResource;

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
        for (int i = 0; i < preyNb; i++)
            cout << preyTypes[i] << " " << preyTypes[i] << "catches"
                 << " ";
        for (int i = 0; i < predatorNb; i++)
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

    void freeMemory() // free memory allocated to landscape table
    {

        for (int row = 0; row < rowNb; row++) // free memory allocated to each row
        {
            delete[] landscapeTablePtr[row];
        }

        delete[] landscapeTablePtr; // free the memory allocated to the array of pointer to each row

        landscapeTablePtr = NULL; // erase the address of the array of pointers to rows
    }
};

int main()
{

    /* enter different types and initital densities */

    predatorTypes = new string[predatorNb]; // assign memory to the pointer
    predatorTypes[0] = "predator1";         // hard coded NOT GOOD

    predatorsInitialDensities = new int[predatorNb];
    predatorsInitialDensities[0] = 5; // hard coded NOT GOOD

    preyTypes = new string[preyNb];
    preyTypes[0] = "prey1"; // hard coded NOT GOOD
    preyTypes[1] = "prey2"; // hard coded NOT GOOD

    preysInitialDensities = new int[preyNb];
    preysInitialDensities[0] = 10; // hard coded NOT GOOD
    preysInitialDensities[1] = 15; // hard coded NOT GOOD

    primaryResourceTypes = new string[resourceNb];
    primaryResourceTypes[0] = "resource1"; // hard coded NOT GOOD
    primaryResourceTypes[1] = "resource2"; // hard coded NOT GOOD

    /* contruct and create landscape */
    landscape world(3, 10);
    world.create();

    /* check if everything is where expected */
    world.getInfo();

    /* check fill function */
    world.fill();
    world.getInfo();

    /* free memory */
    world.freeMemory();
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