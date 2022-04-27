#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <vector>
#include <map>
#include <interface.h>
#include "command.h"

using namespace std;

class Environment
{
public:
    Environment();
    map<string, Interface*> Objects;
    bool Execute(const Command &cmd);

};

#endif // ENVIRONMENT_H
