#include "interface.h"
#include "Utilities.h"

Interface::Interface()
{

}

bool Interface::HasCommand(const string &cmd)
{
    if (aquiutils::lookup(commands(),cmd)!=-1)
        return true;
    else
        return false;
}

vector<string> Interface::commands()
{
    return vector<string>();
}
