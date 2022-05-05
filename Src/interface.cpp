#include "interface.h"
#include "Utilities.h"
#include "environment.h"

Interface::Interface(Environment *_parent)
{
    parent = _parent;
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
