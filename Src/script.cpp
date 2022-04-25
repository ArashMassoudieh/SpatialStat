#include "script.h"
#include "iostream"
#include "fstream"

using namespace std;

Script::Script()
{

}

Script::Script(const string &filename)
{
    Commands.clear();
    GetFromFile(filename);
}

bool Script::GetFromFile(const string &filename)
{
    ifstream file;
    file.open(filename);
    if (!file.good())
    {
        while (!file.eof())
        {   string line;
            getline(file,line);
            Command cmd = Command(line);
        }
    }
    return true;
}
